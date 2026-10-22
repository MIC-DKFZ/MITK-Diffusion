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
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <cppunit/TestAssert.h>

#include <cmath>
#include <string>
#include <vector>

/** Synthetic DWIs and MITK-Data input shared by the IVIM and kurtosis fit tests. */
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
    VectorImageType::Pointer image;
    GradientContainerType::ConstPointer gradients;
    double referenceBValue;
  };

  inline LoadedDwi LoadDwi(const std::string& path)
  {
    PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    auto dwi = IOUtil::Load<Image>(path, &functor);

    LoadedDwi result;
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
