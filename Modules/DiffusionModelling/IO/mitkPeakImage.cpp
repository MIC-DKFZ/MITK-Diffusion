/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkPeakImage.h"
#include "mitkImageDataItem.h"
#include "mitkImageCast.h"
#include <mitkImageToItk.h>
#include <vtkPointData.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkImagePixelReadAccessor.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>

namespace
{
  /** Sets the "FIBER_COLORS" array of the peak lines to one colour. */
  void ColorUniformly(vtkPolyData* polyData, const std::array<float, 3>& rgb)
  {
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(4);
    colors->SetName("FIBER_COLORS");
    colors->SetNumberOfTuples(polyData->GetNumberOfPoints());
    const unsigned char rgba[4] = {static_cast<unsigned char>(rgb[0]), static_cast<unsigned char>(rgb[1]), static_cast<unsigned char>(rgb[2]), 255};
    for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); ++i)
      colors->SetTypedTuple(i, rgba);
    polyData->GetPointData()->AddArray(colors);
  }

  /** Sets the "FIBER_COLORS" array of the peak lines to |direction| per channel. */
  void ColorByDirection(vtkPolyData* polyData)
  {
    vtkPoints* extrPoints = polyData->GetPoints();
    int numOfPoints = 0;
    if (extrPoints!=nullptr)
      numOfPoints = static_cast<int>(extrPoints->GetNumberOfPoints());

    //colors and alpha value for each single point, RGBA = 4 components
    unsigned char rgba[4] = {0,0,0,0};
    int componentSize = 4;

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->ReserveValues(numOfPoints * componentSize);
    colors->SetNumberOfComponents(componentSize);
    colors->SetName("FIBER_COLORS");

    int num_peaks = static_cast<int>(polyData->GetNumberOfLines());
    if (num_peaks < 1)
      return;

    vtkCellArray* fiberList = polyData->GetLines();
    fiberList->InitTraversal();
    for (int fi=0; fi<num_peaks; ++fi)
    {
      vtkIdType const* idList; // contains the point id's of the line
      vtkIdType num_points; // number of points for current line
      fiberList->GetNextCell(num_points, idList);

      if (num_points == 2)
      {
        vnl_vector_fixed< double, 3 > currentPntvtk(extrPoints->GetPoint(idList[0])[0], extrPoints->GetPoint(idList[0])[1],extrPoints->GetPoint(idList[0])[2]);
        vnl_vector_fixed< double, 3 > nextPntvtk(extrPoints->GetPoint(idList[1])[0], extrPoints->GetPoint(idList[1])[1], extrPoints->GetPoint(idList[1])[2]);

        vnl_vector_fixed< double, 3 > diff1;
        diff1 = currentPntvtk - nextPntvtk;
        diff1.normalize();

        rgba[0] = static_cast<unsigned char>(255.0 * std::fabs(diff1[0]));
        rgba[1] = static_cast<unsigned char>(255.0 * std::fabs(diff1[1]));
        rgba[2] = static_cast<unsigned char>(255.0 * std::fabs(diff1[2]));
        rgba[3] = 255;
        colors->InsertTypedTuple(idList[0], rgba);
        colors->InsertTypedTuple(idList[1], rgba);
      }
      else
      {
        MITK_WARN << "Peak with point number!=2 detected!" ;
        continue;
      }
    }

    polyData->GetPointData()->AddArray(colors);
  }
}

mitk::PeakImage::PeakImage() = default;

mitk::PeakImage::PeakImage(const PeakImage& other) : Image(other)
  , m_CustomColor(other.m_CustomColor)
{
}

mitk::PeakImage::~PeakImage() = default;

vtkSmartPointer< vtkPolyData > mitk::PeakImage::GetPolyData() const
{
  std::lock_guard<std::mutex> lock(m_PolyDataMutex);
  // With the geometry's time (BaseData::GetMTime): the lines are in world coordinates.
  const itk::ModifiedTimeType imageTime = this->GetMTime();
  if (m_PolyData == nullptr || imageTime > m_PolyDataTime)
  {
    this->ConstructPolydata();
    m_PolyDataTime = imageTime;
  }
  return m_PolyData;
}

void mitk::PeakImage::KeepPolyData()
{
  std::lock_guard<std::mutex> lock(m_PolyDataMutex);
  m_PolyDataTime = this->GetMTime();
}

void mitk::PeakImage::ConstructPolydata() const
{
  if (this->GetDimensions()[3]%3!=0)
    mitkThrow() << "Fourth dimension needs to be a multiple of 3";

  typedef mitk::ImageToItk< ItkPeakImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(this);
  caster->Update();
  ItkPeakImageType::ConstPointer itkImg = caster->GetOutput();

  int sz_x = static_cast<int>(this->GetDimensions()[0]);
  int sz_y = static_cast<int>(this->GetDimensions()[1]);
  int sz_z = static_cast<int>(this->GetDimensions()[2]);
  unsigned int num_dirs = this->GetDimensions()[3]/3;

  ItkPeakImageType::SpacingType spacing = itkImg->GetSpacing();
  const double minSpacing = std::min({spacing[0], spacing[1], spacing[2]});

  vtkSmartPointer<vtkCellArray> m_VtkCellArray = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints>    m_VtkPoints = vtkSmartPointer<vtkPoints>::New();

  // Serial, so that the lines come in the same order on every run.
  for (int x=0; x<sz_x; ++x)
    for (int y=0; y<sz_y; ++y)
      for (int z=0; z<sz_z; ++z)
      {
        itk::Index<4> idx4; idx4[0] = x; idx4[1] = y; idx4[2] = z;

        itk::ContinuousIndex<double, 4> center;
        center[0] = idx4[0];
        center[1] = idx4[1];
        center[2] = idx4[2];
        center[3] = 0;
        itk::Point<double, 4> worldCenter;
        itkImg->TransformContinuousIndexToPhysicalPoint( center, worldCenter );

        for (unsigned int i=0; i<num_dirs; i++)
        {
          vnl_vector_fixed< float, 3 > dir;
          idx4[3] = i*3;
          dir[0] = itkImg->GetPixel(idx4);
          idx4[3] += 1;
          dir[1] = itkImg->GetPixel(idx4);
          idx4[3] += 1;
          dir[2] = itkImg->GetPixel(idx4);

          // Also skips NaN directions: niftilib sets them to 0, NRRD peak files may contain them.
          if (!(dir.magnitude()>=0.001f))
            continue;

          // add direction to vector field (with spacing compensation)
          itk::Point<double> worldStart;
          worldStart[0] = worldCenter[0]-static_cast<double>(dir[0])/2*minSpacing;
          worldStart[1] = worldCenter[1]-static_cast<double>(dir[1])/2*minSpacing;
          worldStart[2] = worldCenter[2]-static_cast<double>(dir[2])/2*minSpacing;
          itk::Point<double> worldEnd;
          worldEnd[0] = worldCenter[0]+static_cast<double>(dir[0])/2*minSpacing;
          worldEnd[1] = worldCenter[1]+static_cast<double>(dir[1])/2*minSpacing;
          worldEnd[2] = worldCenter[2]+static_cast<double>(dir[2])/2*minSpacing;
          const vtkIdType line[2] = {m_VtkPoints->InsertNextPoint(worldStart.GetDataPointer()),
                                     m_VtkPoints->InsertNextPoint(worldEnd.GetDataPointer())};
          m_VtkCellArray->InsertNextCell(2, line);
        }
      }

  auto polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(m_VtkPoints);
  polyData->SetLines(m_VtkCellArray);

  // Built on demand, also from inside rendering, so the image is not marked as modified here.
  if (m_CustomColor)
    ColorUniformly(polyData, *m_CustomColor);
  else
    ColorByDirection(polyData);
  m_PolyData = polyData;
}


mitk::PeakImage::ItkPeakImageType::Pointer mitk::PeakImage::GetItkPeakImage()
{
  typedef mitk::ImageToItk< mitk::PeakImage::ItkPeakImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(this);
  caster->SetCopyMemFlag(true);
  caster->Update();
  return caster->GetOutput();
}

void mitk::PeakImage::SetCustomColor(float r, float g, float b)
{
  m_CustomColor = std::array<float, 3>{r, g, b};
  // The polydata is built on demand; an image that was never rendered has none yet.
  ColorUniformly(this->GetPolyData(), *m_CustomColor);
  this->Modified();
  this->KeepPolyData();
}

void mitk::PeakImage::ColorByOrientation()
{
  m_CustomColor.reset();
  ColorByDirection(this->GetPolyData());
  this->Modified();
  this->KeepPolyData();
}
