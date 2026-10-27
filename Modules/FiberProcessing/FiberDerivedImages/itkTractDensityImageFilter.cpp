/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "itkTractDensityImageFilter.h"

// VTK
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCell.h>

// misc
#include <algorithm>
#include <cmath>
#include <boost/timer/progress_display.hpp>
#include <vtkBox.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkFiberImageGrid.h>

namespace itk{

template< class OutputImageType, class RefImageType >
TractDensityImageFilter< OutputImageType, RefImageType >::TractDensityImageFilter()
  : m_UpsamplingFactor(1)
  , m_InvertImage(false)
  , m_Mode(TDI_MODE::DENSITY)
  , m_UseImageGeometry(false)
  , m_OutputAbsoluteValues(false)
  , m_MaxDensity(0)
{
}

template< class OutputImageType, class RefImageType >
TDI_MODE TractDensityImageFilter< OutputImageType, RefImageType >::GetMode() const
{
  return m_Mode;
}

template< class OutputImageType, class RefImageType >
void TractDensityImageFilter< OutputImageType, RefImageType >::SetMode(const TDI_MODE &Mode)
{
  m_Mode = Mode;
}

template< class OutputImageType, class RefImageType >
void TractDensityImageFilter< OutputImageType, RefImageType >::GenerateData()
{
  typename OutputImageType::Pointer outImage = this->GetOutput();

  mitk::FiberImageGrid grid;
  if (m_UseImageGeometry && !m_InputImage.IsNull())
  {
    MITK_INFO << "TractDensityImageFilter: using image geometry";
    grid = mitk::FiberImageGrid::FromReference(m_InputImage, m_UpsamplingFactor);
  }
  else
  {
    MITK_INFO << "TractDensityImageFilter: using fiber bundle geometry";
    grid = mitk::FiberImageGrid::FromBundle(m_FiberBundle->GetGeometry(), m_UpsamplingFactor);
  }
  const itk::Vector<double,3> newSpacing = grid.Spacing;
  typename OutputImageType::RegionType::SizeType upsampledSize = grid.Region.GetSize();

  grid.ApplyTo(outImage);
  outImage->Allocate();
  outImage->FillBuffer(0.0);

  int w = upsampledSize[0];
  int h = upsampledSize[1];
  int d = upsampledSize[2];

  // set/initialize output
  OutPixelType* outImageBufferPointer = outImage->GetBufferPointer();

  MITK_INFO << "TractDensityImageFilter: starting image generation";
  vtkSmartPointer<vtkPolyData> fiberPolyData = m_FiberBundle->GetFiberPolyData();

  int numFibers = m_FiberBundle->GetNumFibers();
  boost::timer::progress_display disp(numFibers);
  for( int i=0; i<numFibers; i++ )
  {
    ++disp;
    vtkCell* cell = fiberPolyData->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    float weight = m_FiberBundle->GetFiberWeight(i);

    // A streamline counts once per voxel it visits (Calamante's track density), however many of
    // its segments pass the voxel.
    std::vector< itk::OffsetValueType > visitedVoxels;

    // fill output image
    for( int j=0; j<numPoints-1; j++)
    {
      itk::Point<float, 3> startVertex = mitk::imv::GetItkPoint(points->GetPoint(j));
      itk::Index<3> startIndex;
      itk::ContinuousIndex<float, 3> startIndexCont;
      (void)outImage->TransformPhysicalPointToIndex(startVertex, startIndex);
      (void)outImage->TransformPhysicalPointToContinuousIndex(startVertex, startIndexCont);

      itk::Point<float, 3> endVertex = mitk::imv::GetItkPoint(points->GetPoint(j + 1));
      itk::Index<3> endIndex;
      itk::ContinuousIndex<float, 3> endIndexCont;
      (void)outImage->TransformPhysicalPointToIndex(endVertex, endIndex);
      (void)outImage->TransformPhysicalPointToContinuousIndex(endVertex, endIndexCont);

      std::vector< std::pair< itk::Index<3>, double > > segments = mitk::imv::IntersectImage(newSpacing, startIndex, endIndex, startIndexCont, endIndexCont);
      for (std::pair< itk::Index<3>, double > segment : segments)
      {
        if (!outImage->GetLargestPossibleRegion().IsInside(segment.first))
          continue;
        switch (m_Mode)
        {
        case TDI_MODE::BINARY:
        {
          outImage->SetPixel(segment.first, 1);
          break;
        }
        case TDI_MODE::VISITATION_COUNT:
        {
          visitedVoxels.push_back(outImage->ComputeOffset(segment.first));
          break;
        }
        case TDI_MODE::DENSITY:
        {
          outImage->SetPixel(segment.first, outImage->GetPixel(segment.first)+segment.second * weight);
          break;
        }
        default:
        {
          outImage->SetPixel(segment.first, outImage->GetPixel(segment.first)+segment.second * weight);
        }
        }
      }
    }

    std::sort(visitedVoxels.begin(), visitedVoxels.end());
    visitedVoxels.erase(std::unique(visitedVoxels.begin(), visitedVoxels.end()), visitedVoxels.end());
    for (auto offset : visitedVoxels)
      outImageBufferPointer[offset] += 1;
  }

  m_MaxDensity = 0;
  for (int i=0; i<w*h*d; i++)
    if (m_MaxDensity < outImageBufferPointer[i])
      m_MaxDensity = outImageBufferPointer[i];
  if (!m_OutputAbsoluteValues && m_Mode!=TDI_MODE::BINARY)
  {
    MITK_INFO << "TractDensityImageFilter: max-normalizing output image";
    if (m_MaxDensity>0)
      for (int i=0; i<w*h*d; i++)
      {
        outImageBufferPointer[i] /= m_MaxDensity;
      }
  }
  if (m_InvertImage)
  {
    MITK_INFO << "TractDensityImageFilter: inverting image";
    for (int i=0; i<w*h*d; i++)
      outImageBufferPointer[i] = 1-outImageBufferPointer[i];
  }
  MITK_INFO << "TractDensityImageFilter: finished processing";
}
}

