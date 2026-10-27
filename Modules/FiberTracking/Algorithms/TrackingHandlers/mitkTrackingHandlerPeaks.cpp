/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTrackingHandlerPeaks.h"
#include <algorithm>
#include <cmath>

namespace mitk
{

TrackingHandlerPeaks::TrackingHandlerPeaks()
{

}

TrackingHandlerPeaks::~TrackingHandlerPeaks()
{
}

bool TrackingHandlerPeaks::WorldToIndex(itk::Point<float, 3>& pos, itk::Index<3>& index)
{
  return m_DummyImage->TransformPhysicalPointToIndex(pos, index);
}

void TrackingHandlerPeaks::InitForTracking()
{
  MITK_INFO << "Initializing peak tracker.";

  if (m_NeedsDataInit)
  {
    itk::Vector<double, 4> spacing4 = m_PeakImage->GetSpacing();
    itk::Point<float, 4> origin4 = m_PeakImage->GetOrigin();
    itk::Matrix<double, 4, 4> direction4 = m_PeakImage->GetDirection();
    itk::ImageRegion<4> imageRegion4 = m_PeakImage->GetLargestPossibleRegion();

    spacing3[0] = spacing4[0]; spacing3[1] = spacing4[1]; spacing3[2] = spacing4[2];
    origin3[0] = origin4[0]; origin3[1] = origin4[1]; origin3[2] = origin4[2];
    for (int r=0; r<3; r++)
      for (int c=0; c<3; c++)
      {
        direction3[r][c] = direction4[r][c];
        m_FloatImageRotation[r][c] = direction4[r][c];
      }
    imageRegion3.SetSize(0, imageRegion4.GetSize()[0]);
    imageRegion3.SetSize(1, imageRegion4.GetSize()[1]);
    imageRegion3.SetSize(2, imageRegion4.GetSize()[2]);

    m_DummyImage = ItkUcharImgType::New();
    m_DummyImage->SetSpacing( spacing3 );
    m_DummyImage->SetOrigin( origin3 );
    m_DummyImage->SetDirection( direction3 );
    m_DummyImage->SetRegions( imageRegion3 );
    m_DummyImage->Allocate();
    m_DummyImage->FillBuffer(0.0);

    if (imageRegion4.GetSize(3)==0 || imageRegion4.GetSize(3)%3!=0)
      mitkThrow() << "The 4th dimension of a peak image holds three components per peak, but its size is " << imageRegion4.GetSize(3) << ".";
    m_NumDirs = imageRegion4.GetSize(3)/3;

    m_NeedsDataInit = false;
  }

  this->CalculateMinVoxelSize();
  std::cout << "TrackingHandlerPeaks - Peak threshold: " << m_Parameters->m_Cutoff << std::endl;
  if (m_Parameters->m_Mode == MODE::PROBABILISTIC)
    std::cout << "TrackingHandlerPeaks - Peak jitter: " << m_Parameters->m_PeakJitter << std::endl;
}

vnl_vector_fixed<float,3> TrackingHandlerPeaks::GetMatchingDirection(itk::Index<3> idx3, vnl_vector_fixed<float,3>& oldDir)
{
  vnl_vector_fixed<float,3> out_dir; out_dir.fill(0);
  float angle = 0;
  float mag = oldDir.magnitude();
  if (mag<mitk::eps)
  {
    bool found = false;

    if (!m_Parameters->m_FixRandomSeed)
    {
      // try m_NumDirs times to get a non-zero random direction
      for (int j=0; j<m_NumDirs; j++)
      {
        int i = static_cast<int>(GetThreadRngItk()->GetIntegerVariate(m_NumDirs-1));
        out_dir = GetDirection(idx3, i);

        if (out_dir.magnitude()>mitk::eps)
        {
          oldDir[0] = out_dir[0];
          oldDir[1] = out_dir[1];
          oldDir[2] = out_dir[2];
          found = true;
          break;
        }
      }
    }

    if (!found)
    {
      // if you didn't find a non-zero random direction, take first non-zero direction you find
      for (int i=0; i<m_NumDirs; i++)
      {
        out_dir = GetDirection(idx3, i);
        if (out_dir.magnitude()>mitk::eps)
        {
          oldDir[0] = out_dir[0];
          oldDir[1] = out_dir[1];
          oldDir[2] = out_dir[2];
          found = true;
          break;
        }
      }
    }

    // not the last candidate, which may be NaN (e.g. MRtrix peak files)
    if (!found)
      out_dir.fill(0);
  }
  else
  {
    for (int i=0; i<m_NumDirs; i++)
    {
      vnl_vector_fixed<float,3> dir = GetDirection(idx3, i);
      mag = dir.magnitude();
      if (mag>mitk::eps)
        dir.normalize();
      float a = dot_product(dir, oldDir);
      if (fabs(a)>angle)
      {
        angle = fabs(a);
        if (a<0)
          out_dir = -dir;
        else
          out_dir = dir;
        out_dir *= mag;
        out_dir *= angle; // shrink contribution of direction if is less parallel to previous direction
      }
    }
  }

  return out_dir;
}

vnl_vector_fixed<float,3> TrackingHandlerPeaks::GetDirection(itk::Index<3> idx3, int dirIdx)
{
  vnl_vector_fixed<float,3> dir; dir.fill(0.0);
  if ( !m_DummyImage->GetLargestPossibleRegion().IsInside(idx3) )
    return dir;

  PeakImgType::IndexType idx4;
  idx4.SetElement(0,idx3[0]);
  idx4.SetElement(1,idx3[1]);
  idx4.SetElement(2,idx3[2]);

  for (int k=0; k<3; k++)
  {
    idx4.SetElement(3, dirIdx*3 + k);
    dir[k] = m_PeakImage->GetPixel(idx4);
  }

  if (m_Parameters->m_FlipX)
    dir[0] *= -1;
  if (m_Parameters->m_FlipY)
    dir[1] *= -1;
  if (m_Parameters->m_FlipZ)
    dir[2] *= -1;
  if (m_Parameters->m_ApplyDirectionMatrix)
    dir = m_FloatImageRotation*dir;

  return dir;
}

vnl_vector_fixed<float,3> TrackingHandlerPeaks::GetDirection(itk::Point<float, 3> itkP, bool interpolate, vnl_vector_fixed<float,3> oldDir){
  // transform physical point to index coordinates
  itk::Index<3> idx3;
  itk::ContinuousIndex< float, 3> cIdx;

  vnl_vector_fixed<float,3> dir; dir.fill(0.0);
  if ( !m_DummyImage->TransformPhysicalPointToIndex(itkP, idx3) || !m_DummyImage->TransformPhysicalPointToContinuousIndex(itkP, cIdx) )
    return dir;

  if (interpolate)
  {
    // Trilinear interpolation with the corner indices clamped to the image, as ITK's interpolators
    // do for the FA/GFA, ODF and mask images, so that the outermost half voxel is tracked too.
    itk::Index<3> lo, hi;
    float frac[3];
    for (int d=0; d<3; ++d)
    {
      auto size = static_cast<itk::IndexValueType>(m_DummyImage->GetLargestPossibleRegion().GetSize(d));
      auto base = static_cast<itk::IndexValueType>(std::floor(cIdx[d]));
      frac[d] = 1 - (cIdx[d] - base);   // weight of the lower corner
      lo[d] = std::clamp<itk::IndexValueType>(base, 0, size-1);
      hi[d] = std::clamp<itk::IndexValueType>(base+1, 0, size-1);
    }
    const float frac_x = frac[0], frac_y = frac[1], frac_z = frac[2];
    auto corner = [&lo, &hi](bool x, bool y, bool z) { itk::Index<3> c = {{x ? hi[0] : lo[0], y ? hi[1] : lo[1], z ? hi[2] : lo[2]}}; return c; };

    vnl_vector_fixed<float, 8> interpWeights;
    interpWeights[0] = (  frac_x)*(  frac_y)*(  frac_z);
    interpWeights[1] = (1-frac_x)*(  frac_y)*(  frac_z);
    interpWeights[2] = (  frac_x)*(1-frac_y)*(  frac_z);
    interpWeights[3] = (  frac_x)*(  frac_y)*(1-frac_z);
    interpWeights[4] = (1-frac_x)*(1-frac_y)*(  frac_z);
    interpWeights[5] = (  frac_x)*(1-frac_y)*(1-frac_z);
    interpWeights[6] = (1-frac_x)*(  frac_y)*(1-frac_z);
    interpWeights[7] = (1-frac_x)*(1-frac_y)*(1-frac_z);

    dir = GetMatchingDirection(corner(false, false, false), oldDir) * interpWeights[0];
    dir +=  GetMatchingDirection(corner(true, false, false), oldDir) * interpWeights[1];
    dir +=  GetMatchingDirection(corner(false, true, false), oldDir) * interpWeights[2];
    dir +=  GetMatchingDirection(corner(false, false, true), oldDir) * interpWeights[3];
    dir +=  GetMatchingDirection(corner(true, true, false), oldDir) * interpWeights[4];
    dir +=  GetMatchingDirection(corner(false, true, true), oldDir) * interpWeights[5];
    dir +=  GetMatchingDirection(corner(true, false, true), oldDir) * interpWeights[6];
    dir +=  GetMatchingDirection(corner(true, true, true), oldDir) * interpWeights[7];
  }
  else
    dir = GetMatchingDirection(idx3, oldDir);

  return dir;
}

vnl_vector_fixed<float,3> TrackingHandlerPeaks::ProposeDirection(const itk::Point<float, 3>& pos, std::deque<vnl_vector_fixed<float, 3> >& olddirs, itk::Index<3>& oldIndex)
{
  // CHECK: wann wird wo normalisiert
  vnl_vector_fixed<float,3> output_direction; output_direction.fill(0);

  itk::Index<3> index;
  if (!m_DummyImage->TransformPhysicalPointToIndex(pos, index))
    return output_direction;

  vnl_vector_fixed<float,3> oldDir; oldDir.fill(0.0);
  if (!olddirs.empty())
    oldDir = olddirs.back();
  float old_mag = oldDir.magnitude();

  if (!m_Parameters->m_InterpolateTractographyData && oldIndex==index)
    return oldDir;

  output_direction = GetDirection(pos, m_Parameters->m_InterpolateTractographyData, oldDir);
  float mag = output_direction.magnitude();

  if (mag>=m_Parameters->m_Cutoff)
  {
    if (m_Parameters->m_Mode == MODE::PROBABILISTIC)
    {
      // Each component of the unit peak direction gets normal noise with standard deviation
      // m_PeakJitter (ITK's GetNormalVariate takes the variance), so that the spread depends
      // neither on the peak's orientation nor on its length.
      output_direction.normalize();
      auto rng = GetThreadRngItk();
      double variance = static_cast<double>(m_Parameters->m_PeakJitter)*m_Parameters->m_PeakJitter;
      for (int c=0; c<3; ++c)
        output_direction[c] += rng->GetNormalVariate(0, variance);
    }

    output_direction.normalize();

    float a = 1;
    if (old_mag>0.5)
      a = dot_product(output_direction, oldDir);
    if (a>=m_Parameters->GetAngularThresholdDot())
      output_direction *= mag;
    else
      output_direction.fill(0);
  }
  else
    output_direction.fill(0);

  return output_direction;
}

}

