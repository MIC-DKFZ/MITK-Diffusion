/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "itkTractsToFiberEndingsImageFilter.h"

// VTK
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <boost/timer/progress_display.hpp>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkFiberImageGrid.h>

namespace itk{

  template< class OutputImageType >
  TractsToFiberEndingsImageFilter< OutputImageType >::TractsToFiberEndingsImageFilter()
    : m_UpsamplingFactor(1)
    , m_InvertImage(false)
    , m_UseImageGeometry(false)
    , m_BinaryOutput(false)
    , m_InputImage(nullptr)
  {

  }

  template< class OutputImageType >
  void TractsToFiberEndingsImageFilter< OutputImageType >::GenerateData()
  {
    typename OutputImageType::Pointer outImage = this->GetOutput();

    // The same grid as the tract density image, so that endpoints and density share voxels.
    mitk::FiberImageGrid grid;
    if (m_UseImageGeometry && !m_InputImage.IsNull())
      grid = mitk::FiberImageGrid::FromReference(m_InputImage, m_UpsamplingFactor);
    else
      grid = mitk::FiberImageGrid::FromBundle(m_FiberBundle->GetGeometry(), m_UpsamplingFactor);
    typename OutputImageType::RegionType::SizeType upsampledSize = grid.Region.GetSize();

    grid.ApplyTo(outImage);
    outImage->Allocate();

    int w = upsampledSize[0];
    int h = upsampledSize[1];
    int d = upsampledSize[2];

    // set/initialize output
    OutPixelType* outImageBufferPointer = outImage->GetBufferPointer();
    for (int i=0; i<w*h*d; i++)
      outImageBufferPointer[i] = 0;

    // resample fiber bundle
    vtkSmartPointer<vtkPolyData> fiberPolyData = m_FiberBundle->GetFiberPolyData();

    int numFibers = m_FiberBundle->GetNumFibers();
    boost::timer::progress_display disp(numFibers);
    for( int i=0; i<numFibers; i++ )
    {
        ++disp;
      vtkCell* cell = fiberPolyData->GetCell(i);
      int numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      // fill output image
      if (numPoints>0)
      {
        itk::Point<float, 3> vertex = mitk::imv::GetItkPoint(points->GetPoint(0));
        itk::Index<3> index;

        if (outImage->TransformPhysicalPointToIndex(vertex, index))
        {
          if (m_BinaryOutput)
          {
            outImage->SetPixel(index, 1);
          }
          else
          {
            outImage->SetPixel(index, outImage->GetPixel(index)+1);
          }
        }
      }

      if (numPoints>=2)
      {
        itk::Point<float, 3> vertex = mitk::imv::GetItkPoint(points->GetPoint(numPoints-1));
        itk::Index<3> index;

        if (outImage->TransformPhysicalPointToIndex(vertex, index))
        {
          if (m_BinaryOutput)
          {
              outImage->SetPixel(index, 1);
          }
          else
          {
              outImage->SetPixel(index, outImage->GetPixel(index)+1);
          }
        }
      }
    }

    if (m_InvertImage)
      for (int i=0; i<w*h*d; i++)
        outImageBufferPointer[i] = 1-outImageBufferPointer[i];
  }
}
