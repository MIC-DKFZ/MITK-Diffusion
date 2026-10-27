/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTrackingHandlerTensor.h"
#include <algorithm>
#include <cmath>

namespace mitk
{

TrackingHandlerTensor::TrackingHandlerTensor()
  : m_InterpolateTensors(true)
  , m_NumberOfInputs(0)
{
  m_FaInterpolator = itk::LinearInterpolateImageFunction< itk::Image< float, 3 >, float >::New();
}

TrackingHandlerTensor::~TrackingHandlerTensor()
{
}

void TrackingHandlerTensor::InitForTracking()
{
  MITK_INFO << "Initializing tensor tracker.";

  if (m_NeedsDataInit)
  {
    m_NumberOfInputs = m_TensorImages.size();
    // all tensor images are read at the voxel index of the first one
    for (const auto& tensorImage : m_TensorImages)
      if (tensorImage->GetLargestPossibleRegion()!=m_TensorImages.at(0)->GetLargestPossibleRegion()
          || !tensorImage->GetSpacing().GetVnlVector().is_equal(m_TensorImages.at(0)->GetSpacing().GetVnlVector(), 1e-6)
          || tensorImage->GetOrigin().EuclideanDistanceTo(m_TensorImages.at(0)->GetOrigin())>1e-6
          || !tensorImage->GetDirection().GetVnlMatrix().is_equal(m_TensorImages.at(0)->GetDirection().GetVnlMatrix(), 1e-6))
        mitkThrow() << "All tensor images used for multi-tensor tracking need the same geometry.";
    for (int i=0; i<m_NumberOfInputs; i++)
    {
      ItkPDImgType::Pointer pdImage = ItkPDImgType::New();
      pdImage->SetSpacing( m_TensorImages.at(0)->GetSpacing() );
      pdImage->SetOrigin( m_TensorImages.at(0)->GetOrigin() );
      pdImage->SetDirection( m_TensorImages.at(0)->GetDirection() );
      pdImage->SetRegions( m_TensorImages.at(0)->GetLargestPossibleRegion() );
      pdImage->Allocate();
      m_PdImage.push_back(pdImage);

      ItkDoubleImgType::Pointer emaxImage = ItkDoubleImgType::New();
      emaxImage->SetSpacing( m_TensorImages.at(0)->GetSpacing() );
      emaxImage->SetOrigin( m_TensorImages.at(0)->GetOrigin() );
      emaxImage->SetDirection( m_TensorImages.at(0)->GetDirection() );
      emaxImage->SetRegions( m_TensorImages.at(0)->GetLargestPossibleRegion() );
      emaxImage->Allocate();
      emaxImage->FillBuffer(1.0);
      m_EmaxImage.push_back(emaxImage);
    }

    bool useUserFaImage = true;
    if (m_FaImage.IsNull())
    {
      m_FaImage = ItkFloatImgType::New();
      m_FaImage->SetSpacing( m_TensorImages.at(0)->GetSpacing() );
      m_FaImage->SetOrigin( m_TensorImages.at(0)->GetOrigin() );
      m_FaImage->SetDirection( m_TensorImages.at(0)->GetDirection() );
      m_FaImage->SetRegions( m_TensorImages.at(0)->GetLargestPossibleRegion() );
      m_FaImage->Allocate();
      m_FaImage->FillBuffer(0.0);
      useUserFaImage = false;
    }

    typedef itk::DiffusionTensor3D<float>    TensorType;

    for (int x=0; x<(int)m_TensorImages.at(0)->GetLargestPossibleRegion().GetSize()[0]; x++)
      for (int y=0; y<(int)m_TensorImages.at(0)->GetLargestPossibleRegion().GetSize()[1]; y++)
        for (int z=0; z<(int)m_TensorImages.at(0)->GetLargestPossibleRegion().GetSize()[2]; z++)
        {
          ItkTensorImageType::IndexType index;
          index[0] = x; index[1] = y; index[2] = z;
          for (int i=0; i<m_NumberOfInputs; i++)
          {
            TensorType::EigenValuesArrayType eigenvalues;
            TensorType::EigenVectorsMatrixType eigenvectors;

            ItkTensorImageType::PixelType tensor;
            vnl_vector_fixed<float,3> dir;
            tensor = m_TensorImages.at(i)->GetPixel(index);
            tensor.ComputeEigenAnalysis(eigenvalues, eigenvectors);

            dir[0] = eigenvectors(2, 0);
            dir[1] = eigenvectors(2, 1);
            dir[2] = eigenvectors(2, 2);
            if (dir.magnitude()>mitk::eps)
              dir.normalize();
            else
              dir.fill(0.0);

            m_PdImage.at(i)->SetPixel(index, dir);
            if (!useUserFaImage)
              m_FaImage->SetPixel(index, m_FaImage->GetPixel(index)+tensor.GetFractionalAnisotropy());
            // 0 for zero tensors: inf would make the interpolated TEND tensor NaN (0*inf) next to them
            m_EmaxImage.at(i)->SetPixel(index, eigenvalues[2]>0 ? 2/eigenvalues[2] : 0);
          }
          if (!useUserFaImage)
            m_FaImage->SetPixel(index, m_FaImage->GetPixel(index)/m_NumberOfInputs);
        }

    auto double_dir = m_TensorImages.at(0)->GetDirection().GetVnlMatrix();
    for (int r=0; r<3; r++)
      for (int c=0; c<3; c++)
      {
        m_FloatImageRotation[r][c] = double_dir[r][c];
      }

    m_NeedsDataInit = false;
  }


  if (m_Parameters->m_F+m_Parameters->m_G>1.0)
  {
    float temp = m_Parameters->m_F+m_Parameters->m_G;
    m_Parameters->m_F /= temp;
    m_Parameters->m_G /= temp;
  }

  m_FaInterpolator->SetInputImage(m_FaImage);
  this->CalculateMinVoxelSize();

  std::cout << "TrackingHandlerTensor - FA threshold: " << m_Parameters->m_Cutoff << std::endl;
  std::cout << "TrackingHandlerTensor - f: " << m_Parameters->m_F << std::endl;
  std::cout << "TrackingHandlerTensor - g: " << m_Parameters->m_G << std::endl;
}

vnl_vector_fixed<float,3> TrackingHandlerTensor::GetMatchingDirection(itk::Index<3> idx, vnl_vector_fixed<float,3>& oldDir, int& image_num)
{
  vnl_vector_fixed<float,3> out_dir; out_dir.fill(0);
  float angle = 0;
  float mag = oldDir.magnitude();
  if (mag<mitk::eps)
  {
    for (unsigned int i=0; i<m_PdImage.size(); i++)
    {
      out_dir = m_PdImage.at(i)->GetPixel(idx);

      if (out_dir.magnitude()>0.5)
      {
        image_num = i;
        oldDir[0] = out_dir[0];
        oldDir[1] = out_dir[1];
        oldDir[2] = out_dir[2];
        break;
      }
    }
  }
  else
  {
    for (unsigned int i=0; i<m_PdImage.size(); i++)
    {
      vnl_vector_fixed<float,3> dir = m_PdImage.at(i)->GetPixel(idx);

      float a = dot_product(dir, oldDir);
      if (fabs(a)>angle)
      {
        image_num = i;
        angle = fabs(a);
        if (a<0)
          out_dir = -dir;
        else
          out_dir = dir;
        out_dir *= angle; // shrink contribution of direction if is less parallel to previous direction
      }
    }
  }

  return out_dir;
}

bool TrackingHandlerTensor::WorldToIndex(itk::Point<float, 3>& pos, itk::Index<3>& index)
{
  return m_TensorImages.at(0)->TransformPhysicalPointToIndex(pos, index);
}

vnl_vector_fixed<float,3> TrackingHandlerTensor::GetDirection(itk::Point<float, 3> itkP, vnl_vector_fixed<float,3> oldDir, TensorType& tensor)
{
  // transform physical point to index coordinates
  itk::Index<3> idx;
  itk::ContinuousIndex< float, 3> cIdx;

  vnl_vector_fixed<float,3> dir; dir.fill(0.0);
  if ( !m_TensorImages.at(0)->TransformPhysicalPointToIndex(itkP, idx) || !m_TensorImages.at(0)->TransformPhysicalPointToContinuousIndex(itkP, cIdx) )
    return dir;

  int image_num = -1;
  if (!m_Parameters->m_InterpolateTractographyData)
  {
    dir = GetMatchingDirection(idx, oldDir, image_num);
    if (image_num>=0)
      tensor = m_TensorImages[image_num]->GetPixel(idx) * m_EmaxImage[image_num]->GetPixel(idx);
  }
  else
  {
    // Trilinear interpolation with the corner indices clamped to the image, as ITK's interpolators
    // do for the FA image and the masks, so that the outermost half voxel is tracked too.
    itk::Index<3> lo, hi;
    float frac[3];
    for (int d=0; d<3; ++d)
    {
      auto size = static_cast<itk::IndexValueType>(m_TensorImages.at(0)->GetLargestPossibleRegion().GetSize(d));
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

    const bool corners[8][3] = {{false, false, false}, {true, false, false}, {false, true, false}, {false, false, true},
                                {true, true, false}, {false, true, true}, {true, false, true}, {true, true, true}};
    for (int k=0; k<8; ++k)
    {
      itk::Index<3> c = corner(corners[k][0], corners[k][1], corners[k][2]);
      if (k==0)
        dir = GetMatchingDirection(c, oldDir, image_num) * interpWeights[k];
      else
        dir += GetMatchingDirection(c, oldDir, image_num) * interpWeights[k];
      if (image_num>=0)
        tensor += m_TensorImages[image_num]->GetPixel(c) * m_EmaxImage[image_num]->GetPixel(c) * interpWeights[k];
    }
  }

  return dir;
}

vnl_vector_fixed<float,3> TrackingHandlerTensor::GetLargestEigenvector(TensorType& tensor)
{
  vnl_vector_fixed<float,3> dir;
  TensorType::EigenValuesArrayType eigenvalues;
  TensorType::EigenVectorsMatrixType eigenvectors;
  tensor.ComputeEigenAnalysis(eigenvalues, eigenvectors);
  dir[0] = eigenvectors(2, 0);
  dir[1] = eigenvectors(2, 1);
  dir[2] = eigenvectors(2, 2);
  if (dir.magnitude()<mitk::eps)
    dir.fill(0.0);

  return dir;
}

vnl_vector_fixed<float,3> TrackingHandlerTensor::ProposeDirection(const itk::Point<float, 3>& pos, std::deque<vnl_vector_fixed<float, 3> >& olddirs, itk::Index<3>& oldIndex)
{
  vnl_vector_fixed<float,3> output_direction; output_direction.fill(0);
  TensorType tensor; tensor.Fill(0);

  try
  {
    itk::Index<3> index;
    if (!m_TensorImages.at(0)->TransformPhysicalPointToIndex(pos, index))
      return output_direction;

    float fa = mitk::imv::GetImageValue<float>(pos, m_Parameters->m_InterpolateTractographyData, m_FaInterpolator);
    if (fa<m_Parameters->m_Cutoff)
      return output_direction;

    vnl_vector_fixed<float,3> oldDir; oldDir.fill(0.0);
    if (!olddirs.empty())
      oldDir = olddirs.back();

    // same voxel without interpolation: keep the previous direction, still in the tracking frame
    if (!m_Parameters->m_InterpolateTractographyData && oldIndex==index)
      return oldDir;

    FlipDirection(oldDir);  // to the data frame
    float old_mag = oldDir.magnitude();

    output_direction = GetDirection(pos, oldDir, tensor);
    float mag = output_direction.magnitude();

    if (mag>=mitk::eps)
    {
      output_direction.normalize();

      if (old_mag>0.5 && m_Parameters->m_G>mitk::eps)  // TEND tracking
      {

        output_direction[0] = m_Parameters->m_F*output_direction[0] + (1-m_Parameters->m_F)*( (1-m_Parameters->m_G)*oldDir[0] + m_Parameters->m_G*(tensor[0]*oldDir[0] + tensor[1]*oldDir[1] + tensor[2]*oldDir[2]));
        output_direction[1] = m_Parameters->m_F*output_direction[1] + (1-m_Parameters->m_F)*( (1-m_Parameters->m_G)*oldDir[1] + m_Parameters->m_G*(tensor[1]*oldDir[0] + tensor[3]*oldDir[1] + tensor[4]*oldDir[2]));
        output_direction[2] = m_Parameters->m_F*output_direction[2] + (1-m_Parameters->m_F)*( (1-m_Parameters->m_G)*oldDir[2] + m_Parameters->m_G*(tensor[2]*oldDir[0] + tensor[4]*oldDir[1] + tensor[5]*oldDir[2]));
        output_direction.normalize();
      }

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

  }
  catch(...)
  {

  }

  if (m_Parameters->m_FlipX)
    output_direction[0] *= -1;
  if (m_Parameters->m_FlipY)
    output_direction[1] *= -1;
  if (m_Parameters->m_FlipZ)
    output_direction[2] *= -1;
  if (m_Parameters->m_ApplyDirectionMatrix)
    output_direction = m_FloatImageRotation*output_direction;

  return output_direction;
}

}

