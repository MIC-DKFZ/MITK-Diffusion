/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkTotalVariationDenoisingImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <mitkDiffusionImageDenoising.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionImageDenoisingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageDenoisingTestSuite);
  MITK_TEST(Denoise_Gaussian_FiltersEachVolumeOnItsOwn);
  MITK_TEST(Denoise_TotalVariation_FiltersEachVolumeOnItsOwn);
  MITK_TEST(Denoise_NonLocalMeans_KeepsVolumesApart);
  MITK_TEST(Denoise_NotADwi_Throws);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using VolumeType = itk::Image<short, 3>;
  using Parameters = mitk::DiffusionImageDenoising::Parameters;
  using Method = mitk::DiffusionImageDenoising::Method;

  mitk::Image::Pointer m_Dwi;

  static VolumeType::Pointer Volume(const DwiType* dwi, unsigned int index)
  {
    auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiType, VolumeType>::New();
    extractor->SetInput(dwi);
    extractor->SetIndex(index);
    extractor->Update();
    return extractor->GetOutput();
  }

  static double Mean(const VolumeType* volume)
  {
    double sum = 0;
    itk::ImageRegionConstIterator<VolumeType> it(volume, volume->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      sum += it.Get();
    return sum / volume->GetLargestPossibleRegion().GetNumberOfPixels();
  }

  static void AssertEqualVolumes(const VolumeType* expected, const VolumeType* actual)
  {
    itk::ImageRegionConstIterator<VolumeType> e(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<VolumeType> a(actual, actual->GetLargestPossibleRegion());
    for (; !e.IsAtEnd(); ++e, ++a)
      CPPUNIT_ASSERT_EQUAL(e.Get(), a.Get());
  }

  /** The first (b0) and the last volume of the result equal the filter applied to that volume alone. */
  template <class TFilter>
  void CheckVolumeWise(const Parameters& parameters, typename TFilter::Pointer reference)
  {
    auto result = mitk::DiffusionImageDenoising::Denoise(m_Dwi, parameters);
    auto input = PropHelper::GetItkVectorImage(m_Dwi);
    auto output = PropHelper::GetItkVectorImage(result);
    CPPUNIT_ASSERT_EQUAL(input->GetVectorLength(), output->GetVectorLength());

    for (unsigned int index : {0u, input->GetVectorLength() - 1})
    {
      reference->SetInput(Volume(input, index));
      reference->Update();
      AssertEqualVolumes(reference->GetOutput(), Volume(output, index));
    }

    auto gradients = PropHelper::GetGradientContainer(m_Dwi);
    auto resultGradients = PropHelper::GetGradientContainer(result);
    CPPUNIT_ASSERT_EQUAL(gradients->Size(), resultGradients->Size());
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      CPPUNIT_ASSERT(gradients->ElementAt(i) == resultGradients->ElementAt(i));
  }

public:
  void setUp() override
  {
    // test.dwi: volume 0 is the b0 volume, much brighter than the weighted ones.
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void Denoise_Gaussian_FiltersEachVolumeOnItsOwn()
  {
    Parameters parameters;
    parameters.method = Method::Gaussian;
    parameters.gaussVariance = 2;
    auto reference = itk::DiscreteGaussianImageFilter<VolumeType, VolumeType>::New();
    reference->SetVariance(2);
    CheckVolumeWise<itk::DiscreteGaussianImageFilter<VolumeType, VolumeType>>(parameters, reference);
  }

  void Denoise_TotalVariation_FiltersEachVolumeOnItsOwn()
  {
    Parameters parameters;
    parameters.method = Method::TotalVariation;
    parameters.tvLambda = 0.1;
    parameters.tvIterations = 2;
    auto reference = itk::TotalVariationDenoisingImageFilter<VolumeType, VolumeType>::New();
    reference->SetLambda(0.1);
    reference->SetNumberIterations(2);
    CheckVolumeWise<itk::TotalVariationDenoisingImageFilter<VolumeType, VolumeType>>(parameters, reference);
  }

  void Denoise_NonLocalMeans_KeepsVolumesApart()
  {
    // The patch sampling is random, so only the separation of the volumes is checked: the b0
    // volume stays much brighter than the last weighted one.
    Parameters parameters;
    parameters.method = Method::NonLocalMeans;
    parameters.nlmSearchRadius = 1;
    parameters.nlmNumberOfPatches = 4;
    parameters.nlmIterations = 1;
    auto result = mitk::DiffusionImageDenoising::Denoise(m_Dwi, parameters);

    auto input = PropHelper::GetItkVectorImage(m_Dwi);
    auto output = PropHelper::GetItkVectorImage(result);
    const unsigned int last = input->GetVectorLength() - 1;
    const double inputB0 = Mean(Volume(input, 0));
    const double inputLast = Mean(Volume(input, last));
    CPPUNIT_ASSERT(inputB0 > 2 * inputLast);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(inputB0, Mean(Volume(output, 0)), 0.2 * inputB0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(inputLast, Mean(Volume(output, last)), 0.2 * inputB0);
  }

  void Denoise_NotADwi_Throws()
  {
    auto image = mitk::Image::New();
    image->Initialize(mitk::MakeScalarPixelType<short>(), *m_Dwi->GetGeometry());
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageDenoising::Denoise(image, Parameters()), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageDenoising)
