/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "itkFiberCurvatureFilter.h"

#include <mitkDerivedFiberBundle.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <boost/timer/progress_display.hpp>

namespace itk{

FiberCurvatureFilter::FiberCurvatureFilter()
    : m_AngularDeviation(30)
    , m_Distance(10.0)
    , m_RemoveFibers(false)
    , m_UseMedian(false)
{

}

void FiberCurvatureFilter::GenerateData()
{
    vtkSmartPointer<vtkPolyData> inputPoly = m_InputFiberBundle->GetFiberPolyData();
    using Piece = std::vector< vnl_vector_fixed< double, 3 > >;

    // The pieces of each fiber are collected per fiber and assembled in fiber order afterwards, so
    // that the output does not depend on the number of threads.
    const int numFibers = inputPoly->GetNumberOfCells();
    std::vector< std::vector< Piece > > fiberPieces(numFibers);

    MITK_INFO << "Applying curvature threshold";
    boost::timer::progress_display disp(numFibers);
    // GetCell shares one vtkCell between all callers; GetCellPoints with a vtkIdList of its own
    // reads the fibers in parallel once the cell lookup is built.
    if (inputPoly->NeedToBuildCells())
        inputPoly->BuildCells();
#pragma omp parallel for
    for (int i=0; i<numFibers; i++)
    {
        // The progress display is not thread-safe; only its increment is serialized.
#pragma omp critical
        ++disp;

        vtkNew<vtkIdList> idBuffer;
        vtkIdType numIds = 0;
        const vtkIdType* ids = nullptr;
        inputPoly->GetCellPoints(i, numIds, ids, idBuffer);
        std::vector< vnl_vector_fixed< double, 3 > > vertices(numIds);
        for (vtkIdType j=0; j<numIds; j++)
            inputPoly->GetPoint(ids[j], vertices[j].data_block());

        // calculate curvatures
        int numPoints = vertices.size();
        std::vector< Piece >& pieces = fiberPieces[i];
        Piece container;
        for (int j=0; j<numPoints; j++)
        {
            double dist = 0;
            int c = j;
            std::vector< vnl_vector_fixed< double, 3 > > vectors;
            vnl_vector_fixed< double, 3 > meanV; meanV.fill(0.0);
            // Segments (c-1, c) back to the first segment (0, 1).
            while(dist<m_Distance/2 && c>0)
            {
              vnl_vector_fixed< double, 3 > p1 = vertices.at(c-1);
              vnl_vector_fixed< double, 3 > p2 = vertices.at(c);

                vnl_vector_fixed< double, 3 > v = p2-p1;
                dist += v.magnitude();
                v.normalize();
                vectors.push_back(v);
                if (m_UseMedian && c==j)
                    meanV += v;
                else if (!m_UseMedian)
                  meanV += v;
                c--;
            }
            c = j;
            dist = 0;
            while(dist<m_Distance/2 && c<numPoints-1)
            {
              vnl_vector_fixed< double, 3 > p1 = vertices.at(c);
              vnl_vector_fixed< double, 3 > p2 = vertices.at(c+1);

                vnl_vector_fixed< double, 3 > v = p2-p1;
                dist += v.magnitude();
                v.normalize();
                vectors.push_back(v);
                if (m_UseMedian && c==j)
                    meanV += v;
                else if (!m_UseMedian)
                  meanV += v;
                c++;
            }
            meanV.normalize();

            double dev = 0;
            for (auto vec : vectors)
            {
                double angle = dot_product(meanV, vec);
                if (angle>1.0)
                    angle = 1.0;
                if (angle<-1.0)
                    angle = -1.0;
                dev += acos(angle)*180/itk::Math::pi;
            }
            if (vectors.size()>0)
                dev /= vectors.size();

            if (dev<m_AngularDeviation)
                container.push_back(vertices.at(j));
            else
            {
                if (m_RemoveFibers)
                {
                    pieces.clear();
                    container.clear();
                    break;
                }

                // A single point is no fiber.
                if (container.size()>1)
                    pieces.push_back(container);
                container.clear();
            }
        }

        if (container.size()>1)
            pieces.push_back(container);
    }

    // The output arrays are sized before they are filled: a vtkPolyLine per piece and arrays that
    // grew point by point made this sequential part take most of the filter's time.
    vtkIdType numOutPoints = 0;
    vtkIdType numPieces = 0;
    for (const auto& pieces : fiberPieces)
    {
        for (const auto& piece : pieces)
        {
            numOutPoints += piece.size();
            ++numPieces;
        }
    }

    vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
    vtkNewPoints->SetNumberOfPoints(numOutPoints);
    vtkNew<vtkIdTypeArray> offsets;
    offsets->SetNumberOfValues(numPieces+1);
    vtkNew<vtkIdTypeArray> connectivity;
    connectivity->SetNumberOfValues(numOutPoints);
    // Each piece keeps the weight of its fiber (e.g. SIFT2 or FitFibersToImage weights).
    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
    weights->SetNumberOfValues(numPieces);

    vtkIdType pointId = 0;
    vtkIdType pieceId = 0;
    offsets->SetValue(0, 0);
    for (int i=0; i<numFibers; i++)
    {
        for (const auto& piece : fiberPieces[i])
        {
            for (const auto& point : piece)
            {
                vtkNewPoints->SetPoint(pointId, point.data_block());
                connectivity->SetValue(pointId, pointId);
                ++pointId;
            }
            weights->SetValue(pieceId, m_InputFiberBundle->GetFiberWeight(i));
            ++pieceId;
            offsets->SetValue(pieceId, pointId);
        }
    }
    vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
    vtkNewCells->SetData(offsets, connectivity);

    vtkSmartPointer<vtkPolyData> outputPoly = vtkSmartPointer<vtkPolyData>::New();
    outputPoly->SetPoints(vtkNewPoints);
    outputPoly->SetLines(vtkNewCells);
    m_OutputFiberBundle = mitk::FiberBundle::New(outputPoly);
    m_OutputFiberBundle->SetFiberWeights(weights);
    mitk::CopyFiberBundleMetadata(m_InputFiberBundle, m_OutputFiberBundle);
}

}




