/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundle.h"

#include <mitkGeometry3D.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <cstring>

mitk::FiberBundle::FiberBundle( vtkPolyData* fiberPolyData )
  : m_NumFibers(0)
  , m_IsRAS(false)
{
  // Zeroed, because the TrackVis writer writes all 1000 bytes (unset fields would leak heap
  // contents into the file); hdr_size 0 marks "not set yet" for GetTrackVisHeader.
  std::memset(&m_TrackVisHeader, 0, sizeof(m_TrackVisHeader));
  m_FiberWeights = vtkSmartPointer<vtkFloatArray>::New();
  m_FiberWeights->SetName("FIBER_WEIGHTS");

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  if (fiberPolyData != nullptr)
    m_FiberPolyData = fiberPolyData;
  else
  {
    this->m_FiberPolyData->SetPoints(vtkSmartPointer<vtkPoints>::New());
    this->m_FiberPolyData->SetLines(vtkSmartPointer<vtkCellArray>::New());
  }

  this->UpdateFiberGeometry();
  this->ColorFibersByOrientation();
}

/** Used by Clone(): an independent copy of fibers, colours, weights, statistics and metadata,
 *  without the cleaning and colouring of the polydata constructor. */
mitk::FiberBundle::FiberBundle(const FiberBundle& other)
  : BaseData(other)
  , m_FiberPolyData(vtkSmartPointer<vtkPolyData>::New())
  , m_NumFibers(other.m_NumFibers)
  , m_FiberColors(vtkSmartPointer<vtkUnsignedCharArray>::New())
  , m_FiberWeights(vtkSmartPointer<vtkFloatArray>::New())
  , m_FiberLengths(other.m_FiberLengths)
  , m_MinFiberLength(other.m_MinFiberLength)
  , m_MaxFiberLength(other.m_MaxFiberLength)
  , m_MeanFiberLength(other.m_MeanFiberLength)
  , m_MedianFiberLength(other.m_MedianFiberLength)
  , m_LengthStDev(other.m_LengthStDev)
  , m_TrackVisHeader(other.m_TrackVisHeader)
  , m_IsRAS(other.m_IsRAS)
{
  m_FiberPolyData->DeepCopy(other.m_FiberPolyData);
  m_FiberColors->DeepCopy(other.m_FiberColors);
  m_FiberWeights->DeepCopy(other.m_FiberWeights);
  this->RequestUpdate();
}

mitk::FiberBundle::~FiberBundle()
{

}

mitk::FiberBundle::Pointer mitk::FiberBundle::GetDeepCopy() const
{
  return this->Clone();
}

void mitk::FiberBundle::SetFiberPolyData(vtkSmartPointer<vtkPolyData> fiberPD)
{
  // A new object instead of the old one changed in place: the old one may be in use elsewhere
  // (mappers, callers holding GetFiberPolyData()), and its point data, sized for the old points,
  // must not be carried over to the new ones.
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  if (fiberPD != nullptr)
    polyData->CopyStructure(fiberPD);
  this->m_FiberPolyData = polyData;

  m_NumFibers = static_cast<unsigned int>(m_FiberPolyData->GetNumberOfLines());

  UpdateFiberGeometry();
  ColorFibersByOrientation();
}

/*
 * return vtkPolyData
 */
vtkSmartPointer<vtkPolyData> mitk::FiberBundle::GetFiberPolyData() const
{
  return m_FiberPolyData;
}

void mitk::FiberBundle::UpdateFiberGeometry()
{
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(m_FiberPolyData);
  cleaner->PointMergingOff();
  cleaner->Update();
  m_FiberPolyData = cleaner->GetOutput();

  m_FiberLengths.clear();
  m_MeanFiberLength = 0;
  m_MedianFiberLength = 0;
  m_LengthStDev = 0;
  m_NumFibers = static_cast<unsigned int>(m_FiberPolyData->GetNumberOfCells());

  if (m_FiberColors==nullptr || m_FiberColors->GetNumberOfTuples()!=m_FiberPolyData->GetNumberOfPoints())
    this->ColorFibersByOrientation();

  if (m_FiberWeights->GetNumberOfValues()!=m_NumFibers)
  {
    m_FiberWeights = vtkSmartPointer<vtkFloatArray>::New();
    m_FiberWeights->SetName("FIBER_WEIGHTS");
    m_FiberWeights->SetNumberOfValues(m_NumFibers);
    this->SetFiberWeights(1);
  }

  if (m_NumFibers<=0) // no fibers present; apply default geometry
  {
    m_MinFiberLength = 0;
    m_MaxFiberLength = 0;
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetImageGeometry(false);
    float b[] = {0, 1, 0, 1, 0, 1};
    geometry->SetFloatBounds(b);
    SetGeometry(geometry);
    return;
  }
  double b[6];
  m_FiberPolyData->GetBounds(b);

  // calculate statistics
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto p = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();
    float length = 0;
    for (int j=0; j<p-1; j++)
    {
      double p1[3];
      points->GetPoint(j, p1);
      double p2[3];
      points->GetPoint(j+1, p2);

      double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
      length += static_cast<float>(dist);
    }
    m_FiberLengths.push_back(length);
    m_MeanFiberLength += length;
    if (i==0)
    {
      m_MinFiberLength = length;
      m_MaxFiberLength = length;
    }
    else
    {
      if (length<m_MinFiberLength)
        m_MinFiberLength = length;
      if (length>m_MaxFiberLength)
        m_MaxFiberLength = length;
    }
  }
  m_MeanFiberLength /= m_NumFibers;

  std::vector< float > sortedLengths = m_FiberLengths;
  std::sort(sortedLengths.begin(), sortedLengths.end());
  for (unsigned int i=0; i<m_NumFibers; i++)
    m_LengthStDev += (m_MeanFiberLength-sortedLengths.at(i))*(m_MeanFiberLength-sortedLengths.at(i));
  if (m_NumFibers>1)
    m_LengthStDev /= (m_NumFibers-1);
  else
    m_LengthStDev = 0;
  m_LengthStDev = std::sqrt(m_LengthStDev);
  m_MedianFiberLength = sortedLengths.at(m_NumFibers/2);

  mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
  geometry->SetFloatBounds(b);
  this->SetGeometry(geometry);

  // The bounding box describes the reference grid until a header is set.
  if (m_TrackVisHeader.hdr_size==0)
    this->SetTrackVisHeader(dynamic_cast<const mitk::Geometry3D*>(this->GetGeometry()));

  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

float mitk::FiberBundle::GetFiberWeight(unsigned int fiber) const
{
  // vtkFloatArray does not check the index
  if (fiber >= static_cast<unsigned int>(m_FiberWeights->GetNumberOfValues()))
    mitkThrow() << "Fiber index " << fiber << " out of range (" << m_FiberWeights->GetNumberOfValues() << " weights).";
  return m_FiberWeights->GetValue(fiber);
}

void mitk::FiberBundle::SetFiberWeights(float newWeight)
{
  for (int i=0; i<m_FiberWeights->GetNumberOfValues(); i++)
    m_FiberWeights->SetValue(i, newWeight);
}

void mitk::FiberBundle::SetFiberWeights(vtkSmartPointer<vtkFloatArray> weights)
{
  // One weight per fiber; a mismatch used to be logged and ignored, leaving the old weights.
  if (weights==nullptr || m_NumFibers!=weights->GetNumberOfValues())
    mitkThrow() << "Fiber weights need one value per fiber (" << m_NumFibers << " fibers, "
                << (weights==nullptr ? 0 : weights->GetNumberOfValues()) << " weights).";

  for (int i=0; i<weights->GetNumberOfValues(); i++)
    m_FiberWeights->SetValue(i, weights->GetValue(i));

  m_FiberWeights->SetName("FIBER_WEIGHTS");
}

void mitk::FiberBundle::SetFiberWeight(unsigned int fiber, float weight)
{
  if (fiber >= static_cast<unsigned int>(m_FiberWeights->GetNumberOfValues()))
    mitkThrow() << "Fiber index " << fiber << " out of range (" << m_FiberWeights->GetNumberOfValues() << " weights).";
  m_FiberWeights->SetValue(fiber, weight);
}

unsigned int mitk::FiberBundle::GetNumberOfPoints() const
{
  unsigned int points = 0;
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    points += cell->GetNumberOfPoints();
  }
  return points;
}

// reapply selected colorcoding in case PolyData structure has changed
bool mitk::FiberBundle::Equals(const mitk::FiberBundle* fib, double eps) const
{
  if (fib==nullptr)
  {
    MITK_INFO << "Reference bundle is nullptr!";
    return false;
  }

  if (m_NumFibers!=fib->GetNumFibers())
  {
    MITK_INFO << "Unequal number of fibers!";
    MITK_INFO << m_NumFibers << " vs. " << fib->GetNumFibers();
    return false;
  }

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkCell* cell2 = fib->GetFiberPolyData()->GetCell(i);
    auto numPoints2 = cell2->GetNumberOfPoints();
    vtkPoints* points2 = cell2->GetPoints();

    if (numPoints2!=numPoints)
    {
      MITK_INFO << "Unequal number of points in fiber " << i << "!";
      MITK_INFO << numPoints2 << " vs. " << numPoints;
      return false;
    }

    for (int j=0; j<numPoints; j++)
    {
      double* p1 = points->GetPoint(j);
      double* p2 = points2->GetPoint(j);
      if (fabs(p1[0]-p2[0])>eps || fabs(p1[1]-p2[1])>eps || fabs(p1[2]-p2[2])>eps)
      {
        MITK_INFO << "Unequal points in fiber " << i << " at position " << j << "!";
        MITK_INFO << "p1: " << p1[0] << ", " << p1[1] << ", " << p1[2];
        MITK_INFO << "p2: " << p2[0] << ", " << p2[1] << ", " << p2[2];
        return false;
      }
    }
  }

  return true;
}

void mitk::FiberBundle::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  os << indent << "Number of fibers: " << this->GetNumFibers() << std::endl;
  os << indent << "Min. fiber length: " << this->GetMinFiberLength() << std::endl;
  os << indent << "Max. fiber length: " << this->GetMaxFiberLength() << std::endl;
  os << indent << "Mean fiber length: " << this->GetMeanFiberLength() << std::endl;
  os << indent << "Median fiber length: " << this->GetMedianFiberLength() << std::endl;
  os << indent << "STDEV fiber length: " << this->GetLengthStDev() << std::endl;
  os << indent << "Number of points: " << this->GetNumberOfPoints() << std::endl;

  os << indent << "Extent x: " << this->GetGeometry()->GetExtentInMM(0) << "mm" << std::endl;
  os << indent << "Extent y: " << this->GetGeometry()->GetExtentInMM(1) << "mm" << std::endl;
  os << indent << "Extent z: " << this->GetGeometry()->GetExtentInMM(2) << "mm" << std::endl;
  os << indent << "Diagonal: " << this->GetGeometry()->GetDiagonalLength()  << "mm" << std::endl;

  // the header the TrackVis writer would write
  const TrackVis_header header = this->GetTrackVisHeader();
  os << "\nReference geometry:" << std::endl;
  os << indent << "Size: [" << std::defaultfloat << header.dim[0] << " " << header.dim[1] << " " << header.dim[2] << "]" << std::endl;
  os << indent << "Voxel size: [" << header.voxel_size[0] << " " << header.voxel_size[1] << " " << header.voxel_size[2] << "]" << std::endl;
  os << indent << "Origin: [" << header.origin[0] << " " << header.origin[1] << " " << header.origin[2] << "]" << std::endl;
  os << indent << "Matrix: " << std::scientific << std::endl;
  os << indent << "[[" << header.vox_to_ras[0][0] << ", " << header.vox_to_ras[0][1] << ", " << header.vox_to_ras[0][2] << ", " << header.vox_to_ras[0][3] << "]" << std::endl;
  os << indent << " [" << header.vox_to_ras[1][0] << ", " << header.vox_to_ras[1][1] << ", " << header.vox_to_ras[1][2] << ", " << header.vox_to_ras[1][3] << "]" << std::endl;
  os << indent << " [" << header.vox_to_ras[2][0] << ", " << header.vox_to_ras[2][1] << ", " << header.vox_to_ras[2][2] << ", " << header.vox_to_ras[2][3] << "]" << std::endl;
  os << indent << " [" << header.vox_to_ras[3][0] << ", " << header.vox_to_ras[3][1] << ", " << header.vox_to_ras[3][2] << ", " << header.vox_to_ras[3][3] << "]]" << std::defaultfloat << std::endl;

  if (m_FiberWeights!=nullptr && m_FiberWeights->GetNumberOfValues()>0)
  {
    std::vector< float > weights;
    for (int i=0; i<m_FiberWeights->GetNumberOfValues(); i++)
      weights.push_back(m_FiberWeights->GetValue(i));

    std::sort(weights.begin(), weights.end());

    os << "\nFiber weight statistics" << std::endl;
    os << indent << "Min: " << weights.front() << std::endl;
    os << indent << "1% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.01)) << std::endl;
    os << indent << "5% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.05)) << std::endl;
    os << indent << "25% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.25)) << std::endl;
    os << indent << "Median: " << weights.at(static_cast<unsigned long>(weights.size()*0.5)) << std::endl;
    os << indent << "75% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.75)) << std::endl;
    os << indent << "95% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.95)) << std::endl;
    os << indent << "99% quantile: " << weights.at(static_cast<unsigned long>(weights.size()*0.99)) << std::endl;
    os << indent << "Max: " << weights.back() << std::endl;
  }
  else
    os << indent << "\n\nNo fiber weights." << std::endl;

  Superclass::PrintSelf(os, indent);
}

namespace
{
  /** Header fields that describe the reference grid of the given geometry. */
  void FillTrackVisHeader(mitk::FiberBundle::TrackVis_header& header, const mitk::BaseGeometry* geometry)
  {
    // vox_to_ras holds the index-to-RAS affine of the reference grid, as the TrackVis format defines
    // it (MITK's world is LPS: x and y negated). A geometry that is no image grid (the bounding box
    // of a bundle without reference image) gets a grid whose index 0 lies at the box corner, so that
    // all points lie inside the declared volume.
    const vtkMatrix4x4* indexToLps = geometry->GetVtkMatrix();
    double shift[3] = {0, 0, 0};
    double dims[3];
    for (int i=0; i<3; ++i)
    {
      dims[i] = geometry->GetExtent(i);
      if (!geometry->GetImageGeometry())
      {
        shift[i] = geometry->GetBounds()[2*i];
        dims[i] = std::ceil(geometry->GetExtent(i)) + 1;
      }
    }
    double indexToRas[3][4];
    for (int r=0; r<3; ++r)
    {
      const double flip = r<2 ? -1.0 : 1.0;
      indexToRas[r][3] = indexToLps->GetElement(r, 3);
      for (int c=0; c<3; ++c)
      {
        indexToRas[r][c] = flip * indexToLps->GetElement(r, c);
        indexToRas[r][3] += indexToLps->GetElement(r, c) * shift[c];
      }
      indexToRas[r][3] *= flip;
    }

    // Voxel axes are reordered and flipped to point to R, A and S (voxel order "RAS"), the one order
    // that TrackVis readers, MITK's of 2023 to 2026 included, all map the same way.
    int perm[3] = {0, 1, 2};
    double sign[3] = {1, 1, 1};
    bool used[3] = {false, false, false};
    for (int k=0; k<3; ++k)
    {
      int best = -1;
      for (int j=0; j<3; ++j)
        if (!used[j] && (best<0 || std::fabs(indexToRas[k][j])>std::fabs(indexToRas[k][best])))
          best = j;
      perm[k] = best;
      used[best] = true;
      sign[k] = indexToRas[k][best]<0 ? -1.0 : 1.0;
    }

    for (int r=0; r<3; ++r)
    {
      double translation = indexToRas[r][3];
      for (int k=0; k<3; ++k)
      {
        header.vox_to_ras[r][k] = static_cast<float>(sign[k] * indexToRas[r][perm[k]]);
        if (sign[k]<0)
          translation += indexToRas[r][perm[k]] * (dims[perm[k]] - 1);
      }
      header.vox_to_ras[r][3] = static_cast<float>(translation);
    }
    for (int c=0; c<3; ++c)
      header.vox_to_ras[3][c] = 0;
    header.vox_to_ras[3][3] = 1;

    for(int k=0; k<3 ;k++)
    {
      header.dim[k]            = static_cast<short>(dims[perm[k]]);
      header.voxel_size[k]     = static_cast<float>(geometry->GetSpacing()[perm[k]]);
      header.origin[k]         = 0;
    }

    header.n_scalars = 0;
    header.n_properties = 0;
    std::memcpy(header.voxel_order, "RAS", 4);
    header.image_orientation_patient[0] = 1.0;
    header.image_orientation_patient[1] = 0.0;
    header.image_orientation_patient[2] = 0.0;
    header.image_orientation_patient[3] = 0.0;
    header.image_orientation_patient[4] = 1.0;
    header.image_orientation_patient[5] = 0.0;
    header.pad1[0] = 0;
    header.pad1[1] = 0;
    header.pad2[0] = 0;
    header.pad2[1] = 0;
    header.invert_x = 0;
    header.invert_y = 0;
    header.invert_z = 0;
    header.swap_xy = 0;
    header.swap_yz = 0;
    header.swap_zx = 0;
    header.n_count = 0;
    header.version = 2;
    header.hdr_size = 1000;
    std::memcpy(header.id_string, "TRACK", 6);
  }
}

mitk::FiberBundle::TrackVis_header mitk::FiberBundle::GetTrackVisHeader() const
{
  if (m_TrackVisHeader.hdr_size!=0)
    return m_TrackVisHeader;
  // No header set yet (a bundle constructed without fibers): that of the bounding box.
  TrackVis_header header = m_TrackVisHeader;
  if (auto geometry = dynamic_cast<const mitk::Geometry3D*>(this->GetGeometry()))
    FillTrackVisHeader(header, geometry);
  return header;
}

void mitk::FiberBundle::SetTrackVisHeader(const mitk::FiberBundle::TrackVis_header &TrackVisHeader)
{
  m_TrackVisHeader = TrackVisHeader;
}

void mitk::FiberBundle::SetTrackVisHeader(const mitk::BaseGeometry* geometry)
{
  if (geometry!=nullptr)
    FillTrackVisHeader(m_TrackVisHeader, geometry);
}

/* ESSENTIAL IMPLEMENTATION OF SUPERCLASS METHODS */
void mitk::FiberBundle::UpdateOutputInformation()
{

}
void mitk::FiberBundle::SetRequestedRegionToLargestPossibleRegion()
{

}
bool mitk::FiberBundle::RequestedRegionIsOutsideOfTheBufferedRegion()
{
  return false;
}
bool mitk::FiberBundle::VerifyRequestedRegion()
{
  return true;
}
void mitk::FiberBundle::SetRequestedRegion(const itk::DataObject* )
{

}


namespace mitk {

  bool FiberBundle::IsRAS() const
  {
    return m_IsRAS;
  }


  void FiberBundle::setIsRAS(bool newIsRAS)
  {
    m_IsRAS = newIsRAS;
  }

}
