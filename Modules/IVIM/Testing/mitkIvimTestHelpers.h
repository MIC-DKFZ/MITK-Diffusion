/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkIvimTestHelpers_h
#define mitkIvimTestHelpers_h

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <cppunit/TestAssert.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

/** Synthetic DWIs and mitk-diffusion-test-data input shared by the IVIM and kurtosis fit tests. */
namespace mitk::IvimTesting
{
  using VectorImageType = DiffusionPropertyHelper::ImageType;
  using GradientType = DiffusionPropertyHelper::GradientDirectionType;
  using GradientContainerType = DiffusionPropertyHelper::GradientDirectionsContainerType;

  /** The fits use only the b-value, ref. b-value * |g|^2, so the directions just cycle through
      the axes. */
  inline GradientContainerType::Pointer MakeGradients(const std::vector<double>& bValues, double referenceBValue)
  {
    auto gradients = GradientContainerType::New();
    for (unsigned int i = 0; i < bValues.size(); ++i)
    {
      GradientType g(0.0);
      g[i % 3] = std::sqrt(bValues[i] / referenceBValue);
      gradients->InsertElement(i, g);
    }
    return gradients;
  }

  /** One voxel per signal vector, along x; values are rounded to the DWI pixel type. */
  inline VectorImageType::Pointer MakeImage(const std::vector<std::vector<double>>& signals)
  {
    auto image = VectorImageType::New();
    VectorImageType::SizeType size = {{signals.size(), 1, 1}};
    image->SetRegions(size);
    image->SetVectorLength(signals.front().size());
    image->Allocate();
    for (unsigned int x = 0; x < signals.size(); ++x)
    {
      VectorImageType::PixelType pixel(signals[x].size());
      for (unsigned int i = 0; i < signals[x].size(); ++i)
        pixel[i] = static_cast<short>(std::lround(signals[x][i]));
      image->SetPixel({{x, 0, 0}}, pixel);
    }
    return image;
  }

  /** The model of the IVIM filter: D* is added to D in the perfusion compartment. */
  inline std::vector<double> IvimSignal(const std::vector<double>& bValues, double s0, double f, double d, double dStar)
  {
    std::vector<double> signal;
    for (double b : bValues)
      signal.push_back(s0 * ((1 - f) * std::exp(-b * d) + f * std::exp(-b * (d + dStar))));
    return signal;
  }

  inline std::vector<double> KurtosisSignal(const std::vector<double>& bValues, double s0, double d, double k)
  {
    std::vector<double> signal;
    for (double b : bValues)
      signal.push_back(s0 * std::exp(-b * d + b * b * d * d * k / 6.0));
    return signal;
  }

  struct LoadedDwi
  {
    Image::Pointer dwi;
    VectorImageType::Pointer image;
    GradientContainerType::ConstPointer gradients;
    double referenceBValue;
  };

  inline LoadedDwi LoadDwi(const std::string& path)
  {
    PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    auto dwi = IOUtil::Load<Image>(path, &functor);

    LoadedDwi result;
    result.dwi = dwi;
    result.image = VectorImageType::New();
    CastToItkImage(dwi, result.image);
    result.gradients = DiffusionPropertyHelper::GetGradientContainer(dwi);
    result.referenceBValue = DiffusionPropertyHelper::GetReferenceBValue(dwi);
    return result;
  }

  template <class TImage>
  std::vector<double> Values(const TImage* image)
  {
    std::vector<double> values;
    const auto* buffer = image->GetBufferPointer();
    for (std::size_t i = 0; i < image->GetLargestPossibleRegion().GetNumberOfPixels(); ++i)
      values.push_back(buffer[i]);
    return values;
  }

  inline std::vector<double> Values(const Image* map)
  {
    itk::Image<double, 3>::Pointer itkMap;
    CastToItkImage(map, itkMap);
    return Values(itkMap.GetPointer());
  }

  /** A mask with the geometry of the image, 1 at the given voxels and 0 elsewhere. */
  inline Image::Pointer MakeMask(const VectorImageType* image, const std::vector<itk::Index<3>>& voxels)
  {
    auto mask = itk::Image<short, 3>::New();
    mask->CopyInformation(image);
    mask->SetRegions(image->GetLargestPossibleRegion());
    mask->Allocate();
    mask->FillBuffer(0);
    for (const auto& voxel : voxels)
      mask->SetPixel(voxel, 1);
    return GrabItkImageMemory(mask);
  }

  /** A double mask with the geometry of the image, the given values at the given voxels and 0
      elsewhere. */
  inline Image::Pointer MakeMaskWithValues(const VectorImageType* image,
                                           const std::vector<std::pair<itk::Index<3>, double>>& values)
  {
    auto mask = itk::Image<double, 3>::New();
    mask->CopyInformation(image);
    mask->SetRegions(image->GetLargestPossibleRegion());
    mask->Allocate();
    mask->FillBuffer(0);
    for (const auto& [voxel, value] : values)
      mask->SetPixel(voxel, value);
    return GrabItkImageMemory(mask);
  }

  /** The buffer offset of a voxel of the 11 x 11 x n FWE DWI. */
  inline std::size_t FweOffset(const itk::Index<3>& index)
  {
    return index[0] + 11 * index[1] + 121 * index[2];
  }

  /** Moves the origin of an image and widens its spacing by the given amounts in mm. */
  inline void MoveVoxelGrid(itk::ImageBase<3>* image, double originShift, double spacingChange)
  {
    auto origin = image->GetOrigin();
    auto spacing = image->GetSpacing();
    for (unsigned int i = 0; i < 3; ++i)
    {
      origin[i] += originShift;
      spacing[i] += spacingChange;
    }
    image->SetOrigin(origin);
    image->SetSpacing(spacing);
  }

  /** A mask of the FWE DWI, 1 at the given voxels, with the origin of
      FWE/ReferenceFWImage.nrrd of the test data: that image was stored with fewer digits than the
      DWI, so its origin differs from the DWI's by rounding (up to 4.5e-5 mm). */
  inline Image::Pointer MakeMaskWithRoundedOrigin(const VectorImageType* image,
                                                  const std::vector<itk::Index<3>>& voxels,
                                                  const std::string& roundedImagePath)
  {
    auto mask = MakeMask(image, voxels);
    const auto roundedOrigin = IOUtil::Load<Image>(roundedImagePath)->GetGeometry()->GetOrigin();
    const auto origin = mask->GetGeometry()->GetOrigin();
    CPPUNIT_ASSERT(!Equal(origin, roundedOrigin, 1e-5));
    CPPUNIT_ASSERT(Equal(origin, roundedOrigin, 1e-4));
    mask->GetGeometry()->SetOrigin(roundedOrigin);
    return mask;
  }

  /** Recorded map values are compared relative to their size; the maps are float images. */
  inline void AssertClose(double expected, double actual, double relativeTolerance = 1e-5)
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, actual, relativeTolerance * std::abs(expected));
  }

  inline double Mean(const std::vector<double>& values)
  {
    double sum = 0;
    for (double v : values)
      sum += v;
    return sum / values.size();
  }
}

#endif
