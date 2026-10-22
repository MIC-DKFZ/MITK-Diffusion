/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageHelper_h_
#define __mitkDiffusionImageHelper_h_


#include <MitkDiffusionImageExports.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkVectorContainer.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkLinearInterpolateImageFunction.h>
#include <mitkImage.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>

namespace mitk{

class MITKDIFFUSIONIMAGE_EXPORT imv
{
public:

  /** Length of the segment from sf to ef (continuous indices) inside each voxel it crosses; si and ei are the voxels of sf and ef. */
  static std::vector< std::pair< itk::Index<3>, double > > IntersectImage(const itk::Vector<double,3>& spacing, itk::Index<3> si, itk::Index<3> ei, const itk::ContinuousIndex<float, 3>& sf, const itk::ContinuousIndex<float, 3>& ef);

  static itk::Point<float, 3> GetItkPoint(double point[3])
  {
    itk::Point<float, 3> itkPoint;
    itkPoint[0] = static_cast<float>(point[0]);
    itkPoint[1] = static_cast<float>(point[1]);
    itkPoint[2] = static_cast<float>(point[2]);
    return itkPoint;
  }

  template< class TPixelType, class TOutPixelType=TPixelType >
  static TOutPixelType GetImageValue(const itk::Point<float, 3>& itkP, bool interpolate, typename itk::LinearInterpolateImageFunction< itk::Image< TPixelType, 3 >, float >::Pointer interpolator)
  {
    TPixelType nullval = TPixelType(0.0);
    if (interpolator==nullptr)
      return nullval;

    itk::ContinuousIndex< float, 3> cIdx;
    interpolator->ConvertPointToContinuousIndex(itkP, cIdx);

    if (interpolator->IsInsideBuffer(cIdx))
    {
      if (interpolate)
        return interpolator->EvaluateAtContinuousIndex(cIdx);
      else
      {
        itk::Index<3> idx;
        interpolator->ConvertContinuousIndexToNearestIndex(cIdx, idx);
        return interpolator->EvaluateAtIndex(idx);
      }
    }
    else
      return nullval;
  }

  template< class TPixelType=unsigned char >
  static bool IsInsideMask(const itk::Point<float, 3>& itkP, bool interpolate, typename itk::LinearInterpolateImageFunction< itk::Image< TPixelType, 3 >, float >::Pointer interpolator, float threshold=0.5)
  {
    if (interpolator==nullptr)
      return false;

    itk::ContinuousIndex< float, 3> cIdx;
    interpolator->ConvertPointToContinuousIndex(itkP, cIdx);

    if (interpolator->IsInsideBuffer(cIdx))
    {
      float value = 0.0;
      if (interpolate)
        value = interpolator->EvaluateAtContinuousIndex(cIdx);
      else
      {
        itk::Index<3> idx;
        interpolator->ConvertContinuousIndexToNearestIndex(cIdx, idx);
        value = interpolator->EvaluateAtIndex(idx);
      }

      if (value>=threshold)
        return true;
    }
    return false;
  }

  template< class TType=float >
  static vnl_matrix_fixed< TType, 3, 3 > GetRotationMatrixVnl(TType rx, TType ry, TType rz)
  {
    rx = rx*itk::Math::pi/180;
    ry = ry*itk::Math::pi/180;
    rz = rz*itk::Math::pi/180;

    vnl_matrix_fixed< TType, 3, 3 > rotX; rotX.set_identity();
    rotX[1][1] = cos(rx);
    rotX[2][2] = rotX[1][1];
    rotX[1][2] = -sin(rx);
    rotX[2][1] = -rotX[1][2];

    vnl_matrix_fixed< TType, 3, 3 > rotY; rotY.set_identity();
    rotY[0][0] = cos(ry);
    rotY[2][2] = rotY[0][0];
    rotY[0][2] = sin(ry);
    rotY[2][0] = -rotY[0][2];

    vnl_matrix_fixed< TType, 3, 3 > rotZ; rotZ.set_identity();
    rotZ[0][0] = cos(rz);
    rotZ[1][1] = rotZ[0][0];
    rotZ[0][1] = -sin(rz);
    rotZ[1][0] = -rotZ[0][1];

    vnl_matrix_fixed< TType, 3, 3 > rot = rotZ*rotY*rotX;
    return rot;
  }

  template< class TType=float >
  static itk::Matrix< TType, 3, 3 > GetRotationMatrixItk(TType rx, TType ry, TType rz)
  {
    rx = rx*itk::Math::pi/180;
    ry = ry*itk::Math::pi/180;
    rz = rz*itk::Math::pi/180;

    itk::Matrix< TType, 3, 3 > rotX; rotX.SetIdentity();
    rotX[1][1] = cos(rx);
    rotX[2][2] = rotX[1][1];
    rotX[1][2] = -sin(rx);
    rotX[2][1] = -rotX[1][2];

    itk::Matrix< TType, 3, 3 > rotY; rotY.SetIdentity();
    rotY[0][0] = cos(ry);
    rotY[2][2] = rotY[0][0];
    rotY[0][2] = sin(ry);
    rotY[2][0] = -rotY[0][2];

    itk::Matrix< TType, 3, 3 > rotZ; rotZ.SetIdentity();
    rotZ[0][0] = cos(rz);
    rotZ[1][1] = rotZ[0][0];
    rotZ[0][1] = -sin(rz);
    rotZ[1][0] = -rotZ[0][1];

    itk::Matrix< TType, 3, 3 > rot = rotZ*rotY*rotX;
    return rot;
  }

  /**
   * DWI volumes stored as one 4-D image (the NIfTI layout, volume t is the slowest index) as MITK's
   * vector image (volume t is component t of each voxel), with the geometry of the first three
   * dimensions.
   */
  template< class TPixelType >
  static typename itk::VectorImage< TPixelType, 3 >::Pointer VolumesToVectorImage(const itk::Image< TPixelType, 4 >* volumes)
  {
    using VectorImageType = itk::VectorImage< TPixelType, 3 >;
    const auto& region4 = volumes->GetBufferedRegion();
    typename VectorImageType::RegionType region;
    typename VectorImageType::SpacingType spacing;
    typename VectorImageType::PointType origin;
    typename VectorImageType::DirectionType direction;
    for (unsigned int i = 0; i < 3; ++i)
    {
      region.SetIndex(i, region4.GetIndex(i));
      region.SetSize(i, region4.GetSize(i));
      spacing[i] = volumes->GetSpacing()[i];
      origin[i] = volumes->GetOrigin()[i];
      for (unsigned int j = 0; j < 3; ++j)
        direction[i][j] = volumes->GetDirection()[i][j];
    }

    auto image = VectorImageType::New();
    image->SetRegions(region);
    image->SetSpacing(spacing);
    image->SetOrigin(origin);
    image->SetDirection(direction);
    image->SetVectorLength(region4.GetSize(3));
    image->Allocate();

    // The buffers are transposed directly: the 4-D image stores one volume after the other, the
    // vector image all volumes of a voxel next to each other. GetPixel/SetPixel per value takes
    // about ten times as long as the NIfTI file IO itself.
    const std::size_t numberOfVolumes = region4.GetSize(3);
    const std::size_t numberOfVoxels = region.GetNumberOfPixels();
    const TPixelType* in = volumes->GetBufferPointer();
    TPixelType* out = image->GetBufferPointer();
    for (std::size_t v = 0; v < numberOfVoxels; ++v)
      for (std::size_t t = 0; t < numberOfVolumes; ++t)
        out[v * numberOfVolumes + t] = in[t * numberOfVoxels + v];
    return image;
  }

  /** Inverse of VolumesToVectorImage; the fourth dimension gets spacing 1, origin 0 and no rotation. */
  template< class TPixelType >
  static typename itk::Image< TPixelType, 4 >::Pointer VectorImageToVolumes(const itk::VectorImage< TPixelType, 3 >* image)
  {
    using VolumesType = itk::Image< TPixelType, 4 >;
    const auto& region = image->GetBufferedRegion();
    typename VolumesType::RegionType region4;
    typename VolumesType::SpacingType spacing4;
    typename VolumesType::PointType origin4;
    typename VolumesType::DirectionType direction4;
    spacing4.Fill(1);
    origin4.Fill(0);
    direction4.SetIdentity();
    for (unsigned int i = 0; i < 3; ++i)
    {
      region4.SetIndex(i, region.GetIndex(i));
      region4.SetSize(i, region.GetSize(i));
      spacing4[i] = image->GetSpacing()[i];
      origin4[i] = image->GetOrigin()[i];
      for (unsigned int j = 0; j < 3; ++j)
        direction4[i][j] = image->GetDirection()[i][j];
    }
    region4.SetIndex(3, 0);
    region4.SetSize(3, image->GetVectorLength());

    auto volumes = VolumesType::New();
    volumes->SetRegions(region4);
    volumes->SetSpacing(spacing4);
    volumes->SetOrigin(origin4);
    volumes->SetDirection(direction4);
    volumes->Allocate();

    // Voxel by voxel as in VolumesToVectorImage: the vector image is read sequentially, each volume
    // is written sequentially.
    const std::size_t numberOfVolumes = image->GetVectorLength();
    const std::size_t numberOfVoxels = region.GetNumberOfPixels();
    const TPixelType* in = image->GetBufferPointer();
    TPixelType* out = volumes->GetBufferPointer();
    for (std::size_t v = 0; v < numberOfVoxels; ++v)
      for (std::size_t t = 0; t < numberOfVolumes; ++t)
        out[t * numberOfVoxels + v] = in[v * numberOfVolumes + t];
    return volumes;
  }
};

class MITKDIFFUSIONIMAGE_EXPORT gradients
{
private:
  using IndiciesVector = std::vector<unsigned int>;
  using BValueMap = mitk::BValueMapProperty::BValueMap;

  using GradientDirectionContainerType = DiffusionPropertyHelper::GradientDirectionsContainerType;
  using GradientDirectionType = DiffusionPropertyHelper::GradientDirectionType;

public:

  static void Cart2Sph(double x, double y, double z, double* spherical);
  static vnl_vector_fixed<double, 3> Sph2Cart(const double& theta, const double& phi, const double& rad);
  static GradientDirectionContainerType::ConstPointer ReadBvalsBvecs(std::string bvals_file, std::string bvecs_file, double& reference_bval);
  static void WriteBvalsBvecs(std::string bvals_file, std::string bvecs_file, GradientDirectionContainerType::ConstPointer gradients, double reference_bval);
  /**
   * Measurement frame of FSL bvecs for an image with the given direction matrix: applied to the bvecs
   * it gives world (LPS) gradients. ReadBvalsBvecs and WriteBvalsBvecs do not apply it.
   */
  static DiffusionPropertyHelper::MeasurementFrameType FslBvecsMeasurementFrame(const itk::Matrix<double, 3, 3>& direction);
  static std::vector<unsigned int> GetAllUniqueDirections(const BValueMap &bValueMap, GradientDirectionContainerType::ConstPointer refGradientsContainer );

  static vnl_matrix<double> ComputeSphericalFromCartesian(const IndiciesVector  & refShell, const GradientDirectionContainerType * refGradientsContainer);
  static GradientDirectionContainerType::Pointer CreateNormalizedUniqueGradientDirectionContainer(const BValueMap &bValueMap, const GradientDirectionContainerType * origninalGradentcontainer);
  };

}

#endif //__mitkDiffusionImageHelper_h_

