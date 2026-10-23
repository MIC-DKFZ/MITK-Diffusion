/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "itkTractsToRgbaImageFilter.h"

// VTK
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

// misc
#include <math.h>
#include <boost/timer/progress_display.hpp>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkFiberImageGrid.h>

namespace itk{

template< class OutputImageType >
TractsToRgbaImageFilter< OutputImageType >::TractsToRgbaImageFilter()
  : m_UpsamplingFactor(1)
  , m_UseImageGeometry(false)
  , m_InputImage(nullptr)
{

}

template< class OutputImageType >
TractsToRgbaImageFilter< OutputImageType >::~TractsToRgbaImageFilter()
{
}

template< class OutputImageType >
void TractsToRgbaImageFilter< OutputImageType >::GenerateData()
{
  if(&typeid(OutPixelType) != &typeid(itk::RGBAPixel<unsigned char>))
    return;

  typename OutputImageType::Pointer outImage = this->GetOutput();

  // The same grid as the tract density image, so that the derived images share voxels.
  mitk::FiberImageGrid grid;
  if (m_UseImageGeometry && !m_InputImage.IsNull())
    grid = mitk::FiberImageGrid::FromReference(m_InputImage, m_UpsamplingFactor);
  else
    grid = mitk::FiberImageGrid::FromBundle(m_FiberBundle->GetGeometry(), m_UpsamplingFactor);
  const itk::Vector<double,3> newSpacing = grid.Spacing;
  typename OutputImageType::RegionType::SizeType upsampledSize = grid.Region.GetSize();

  grid.ApplyTo(outImage);
  outImage->Allocate();
  OutPixelType nullp;
  nullp.Fill(0.0);
  outImage->FillBuffer(nullp);

  itk::Image< itk::RGBAPixel<double>, 3 >::Pointer double_out = itk::Image< itk::RGBAPixel<double>, 3 >::New();
  grid.ApplyTo(double_out);
  double_out->Allocate();
  auto p = itk::RGBAPixel<double>();
  p.Fill(0.0);
  double_out->FillBuffer(p);

  vtkSmartPointer<vtkPolyData> fiberPolyData = m_FiberBundle->GetFiberPolyData();
  int numFibers = m_FiberBundle->GetNumFibers();

  float scale = 100 * pow((float)m_UpsamplingFactor,3);
  boost::timer::progress_display disp(numFibers);
  for( int i=0; i<numFibers; ++i )
  {
    ++disp;
    vtkCell* cell = fiberPolyData->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();
    const float weight = m_FiberBundle->GetFiberWeight(i);

    for( int j=0; j<numPoints-1; ++j)
    {
      Point<float, 3> startVertex = mitk::imv::GetItkPoint(points->GetPoint(j));
      Index<3> startIndex;
      ContinuousIndex<float, 3> startIndexCont;
      (void)double_out->TransformPhysicalPointToIndex(startVertex, startIndex);
      (void)double_out->TransformPhysicalPointToContinuousIndex(startVertex, startIndexCont);

      Point<float, 3> endVertex = mitk::imv::GetItkPoint(points->GetPoint(j + 1));
      Index<3> endIndex;
      ContinuousIndex<float, 3> endIndexCont;
      (void)double_out->TransformPhysicalPointToIndex(endVertex, endIndex);
      (void)double_out->TransformPhysicalPointToContinuousIndex(endVertex, endIndexCont);

      Vector<float, 3> dir;
      dir[0] = fabs(endVertex[0]-startVertex[0]);
      dir[1] = fabs(endVertex[1]-startVertex[1]);
      dir[2] = fabs(endVertex[2]-startVertex[2]);
      dir.Normalize();

      std::vector< std::pair< itk::Index<3>, double > > segments = mitk::imv::IntersectImage(newSpacing, startIndex, endIndex, startIndexCont, endIndexCont);
      for (std::pair< itk::Index<3>, double > segment : segments)
      {
        if (!double_out->GetLargestPossibleRegion().IsInside(segment.first))
          continue;

        itk::RGBAPixel<double> pix = double_out->GetPixel(segment.first);
        // Direction-encoded track density: |direction| and opacity weighted by the in-voxel
        // length and the fiber weight, like the density image, not by the number of segments.
        const double amount = segment.second * weight * scale;
        pix[0] += dir[0] * amount;
        pix[1] += dir[1] * amount;
        pix[2] += dir[2] * amount;
        pix[3] += amount;

        double_out->SetPixel(segment.first, pix);
      }
    }

    float maxRgb = 0.000000001;
    float maxInt = 0.000000001;
    int w = upsampledSize[0];
    int h = upsampledSize[1];
    int d = upsampledSize[2];
    int numPix = w*h*d*4;

    double* buffer = (double*)double_out->GetBufferPointer();
    // calc maxima
    for(int i=0; i<numPix; i++)
    {
      if((i-3)%4 != 0)
      {
        if(buffer[i] > maxRgb)
          maxRgb = buffer[i];
      }
      else
      {
        if(buffer[i] > maxInt)
          maxInt = buffer[i];
      }
    }

    // write output, normalized uchar 0..255
    unsigned char* outImageBufferPointer = (unsigned char*)outImage->GetBufferPointer();
    for(int i=0; i<numPix; i++)
    {
      if((i-3)%4 != 0)
        outImageBufferPointer[i] = (unsigned char) (255.0 * buffer[i] / maxRgb);
      else
        outImageBufferPointer[i] = (unsigned char) (255.0 * buffer[i] / maxInt);
    }
  }
}
}
