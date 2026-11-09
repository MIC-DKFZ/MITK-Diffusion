/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkBSplineInterpolateImageFunction.h>
#include <itkCastImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleDwiImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

#include <algorithm>
#include <cmath>

class mitkDiffusionResampleDwiImageFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionResampleDwiImageFilterTestSuite);
  MITK_TEST(Resample_IdentityDirection_KeepsCorner);
  MITK_TEST(Resample_FlippedAxis_KeepsCorner);
  MITK_TEST(Resample_EachSetterAndInterpolation_ResamplesEveryChannelAlike);
  MITK_TEST(Resample_SetterBeforeInputAndUpstreamFilter_RunsPipeline);
  MITK_TEST(ResampledImageGrid_RoundingErrorBelowInteger_KeepsVoxel);
  MITK_TEST(Resample_Linear_RoundsValues);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::ResampleDwiImageFilter<short>;
  using ChannelType = itk::Image<short, 3>;

  /**
   * The filter's grid and interpolators applied to one channel with itk::ResampleImageFilter, in
   * double; the filter rounds the values into its short output.
   */
  static itk::Image<double, 3>::Pointer ResampleChannel(const DwiType* input, unsigned int channel, const DwiType* output, FilterType::Interpolation interpolation)
  {
    auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiType, ChannelType>::New();
    extractor->SetInput(input);
    extractor->SetIndex(channel);
    auto resampler = itk::ResampleImageFilter<ChannelType, itk::Image<double, 3>>::New();
    resampler->SetInput(extractor->GetOutput());
    resampler->SetSize(output->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(output->GetSpacing());
    resampler->SetOutputOrigin(itk::ResampledImageOrigin(input->GetOrigin(), input->GetDirection(), input->GetSpacing(), output->GetSpacing()));
    resampler->SetOutputDirection(input->GetDirection());
    switch (interpolation)
    {
      case FilterType::Interpolate_NearestNeighbour:
        resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ChannelType>::New());
        break;
      case FilterType::Interpolate_BSpline:
        resampler->SetInterpolator(itk::BSplineInterpolateImageFunction<ChannelType>::New());
        break;
      case FilterType::Interpolate_WindowedSinc:
        resampler->SetInterpolator(itk::WindowedSincInterpolateImageFunction<ChannelType, 3>::New());
        break;
      default:
        resampler->SetInterpolator(itk::LinearInterpolateImageFunction<ChannelType>::New());
    }
    resampler->Update();
    return resampler->GetOutput();
  }

  static void CheckChannels(const DwiType* input, const DwiType* output, FilterType::Interpolation interpolation)
  {
    CPPUNIT_ASSERT_EQUAL(input->GetVectorLength(), output->GetVectorLength());
    CPPUNIT_ASSERT(output->GetDirection() == input->GetDirection());
    for (unsigned int c = 0; c < input->GetVectorLength(); c += 7)
    {
      auto expected = ResampleChannel(input, c, output, interpolation);
      itk::ImageRegionConstIteratorWithIndex<DwiType> it(output, output->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
        CPPUNIT_ASSERT_EQUAL(static_cast<short>(std::clamp(std::round(expected->GetPixel(it.GetIndex())), -32768.0, 32767.0)), it.Get()[c]);
    }
  }

  DwiType::Pointer m_Dwi;

  /** World position of the outer corner of voxel 0. */
  static itk::Point<double, 3> Corner(const DwiType* image)
  {
    itk::Point<double, 3> corner = image->GetOrigin();
    corner -= image->GetDirection() * (image->GetSpacing() * 0.5);
    return corner;
  }

  /**
   * Upsampling by 2 with nearest-neighbour interpolation copies input voxel i/2 into output voxel
   * i, if the outer corners of both grids coincide.
   */
  void CheckUpsampledByTwo(const DwiType* input)
  {
    auto filter = FilterType::New();
    filter->SetInput(input);
    itk::Vector<double, 3> factor;
    factor.Fill(2.0);
    filter->SetSamplingFactor(factor);
    filter->SetInterpolation(FilterType::Interpolate_NearestNeighbour);
    filter->Update();
    auto output = filter->GetOutput();

    const auto expectedCorner = Corner(input);
    const auto corner = Corner(output);
    for (unsigned int d = 0; d < 3; ++d)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedCorner[d], corner[d], 1e-9);
      CPPUNIT_ASSERT_EQUAL(2 * input->GetLargestPossibleRegion().GetSize(d), output->GetLargestPossibleRegion().GetSize(d));
    }
    CPPUNIT_ASSERT(output->GetDirection() == input->GetDirection());

    itk::ImageRegionConstIteratorWithIndex<DwiType> it(output, output->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      DwiType::IndexType inputIndex;
      for (unsigned int d = 0; d < 3; ++d)
        inputIndex[d] = it.GetIndex()[d] / 2;
      CPPUNIT_ASSERT(it.Get() == input->GetPixel(inputIndex));
    }
  }

public:
  void setUp() override
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
    m_Dwi = PropHelper::GetItkVectorImage(image);
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void Resample_IdentityDirection_KeepsCorner()
  {
    CheckUpsampledByTwo(m_Dwi);
  }

  void Resample_FlippedAxis_KeepsCorner()
  {
    // A flipped y axis as ITK reads it from most clinical NIfTI files.
    DwiType::DirectionType direction;
    direction.SetIdentity();
    direction[1][1] = -1;
    m_Dwi->SetDirection(direction);
    CheckUpsampledByTwo(m_Dwi);
  }

  void Resample_EachSetterAndInterpolation_ResamplesEveryChannelAlike()
  {
    // test.dwi: 11 x 11 x 3 voxels of 2 mm, flipped y axis.
    DwiType::DirectionType direction;
    direction.SetIdentity();
    direction[1][1] = -1;
    m_Dwi->SetDirection(direction);

    const std::vector<FilterType::Interpolation> interpolations = { FilterType::Interpolate_NearestNeighbour, FilterType::Interpolate_Linear, FilterType::Interpolate_BSpline, FilterType::Interpolate_WindowedSinc };
    for (auto interpolation : interpolations)
    {
      // Sampling factor: sizes 11 * 1.5 and 3 * 1.5 are truncated to 16 and 4.
      auto filter = FilterType::New();
      filter->SetInput(m_Dwi);
      itk::Vector<double, 3> factor;
      factor[0] = 1.5;
      factor[1] = 0.7;
      factor[2] = 1.5;
      filter->SetSamplingFactor(factor);
      filter->SetInterpolation(interpolation);
      filter->Update();
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(16), filter->GetOutput()->GetLargestPossibleRegion().GetSize(0));
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(7), filter->GetOutput()->GetLargestPossibleRegion().GetSize(1));
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(4), filter->GetOutput()->GetLargestPossibleRegion().GetSize(2));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0 / 1.5, filter->GetOutput()->GetSpacing()[0], 1e-12);
      CheckChannels(m_Dwi, filter->GetOutput(), interpolation);

      // New spacing: 22 mm / 3 mm = 7.33 voxels are truncated to 7.
      filter = FilterType::New();
      filter->SetInput(m_Dwi);
      itk::Vector<double, 3> spacing;
      spacing.Fill(3.0);
      filter->SetNewSpacing(spacing);
      filter->SetInterpolation(interpolation);
      filter->Update();
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(7), filter->GetOutput()->GetLargestPossibleRegion().GetSize(0));
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(2), filter->GetOutput()->GetLargestPossibleRegion().GetSize(2));
      CPPUNIT_ASSERT(filter->GetOutput()->GetSpacing() == spacing);
      CheckChannels(m_Dwi, filter->GetOutput(), interpolation);

      // New size: the spacing follows from it.
      filter = FilterType::New();
      filter->SetInput(m_Dwi);
      itk::ImageRegion<3> region;
      region.SetSize(0, 5);
      region.SetSize(1, 20);
      region.SetSize(2, 3);
      filter->SetNewImageSize(region);
      filter->SetInterpolation(interpolation);
      filter->Update();
      CPPUNIT_ASSERT(filter->GetOutput()->GetLargestPossibleRegion().GetSize() == region.GetSize());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0 / 5, filter->GetOutput()->GetSpacing()[0], 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0 / 20, filter->GetOutput()->GetSpacing()[1], 1e-12);
      CheckChannels(m_Dwi, filter->GetOutput(), interpolation);
    }
  }

  // Sizes that are integers in exact arithmetic came out just below them and were truncated: an
  // oblique 2 mm image (spacing read as 1.9999999999999998) resampled to 1 mm had 21 instead of 22
  // voxels along x, 1.2 mm -> 0.4 mm gave 32 instead of 33.
  void ResampledImageGrid_RoundingErrorBelowInteger_KeepsVoxel()
  {
    itk::ImageRegion<3> region;
    region.SetSize({{11, 11, 11}});
    itk::Vector<double, 3> spacing;
    spacing.Fill(1.2);
    itk::Vector<double, 3> value;
    value.Fill(0.4);
    itk::Vector<double, 3> newSpacing;
    itk::ImageRegion<3> newRegion;
    itk::ResampledImageGrid(spacing, region, itk::ResampleGridMode::Spacing, value, newSpacing, newRegion);
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(33), newRegion.GetSize(0));

    spacing.Fill(2);
    spacing[0] = 1.9999999999999998;
    value.Fill(1);
    itk::ResampledImageGrid(spacing, region, itk::ResampleGridMode::Spacing, value, newSpacing, newRegion);
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(22), newRegion.GetSize(i));

    spacing.Fill(1.2);
    value.Fill(3);
    itk::ResampledImageGrid(spacing, region, itk::ResampleGridMode::SamplingFactor, value, newSpacing, newRegion);
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(33), newRegion.GetSize(0));

    // Sizes that are not integers are still truncated: 11 * 2 / 3 = 7.33.
    spacing.Fill(2);
    value.Fill(3);
    itk::ResampledImageGrid(spacing, region, itk::ResampleGridMode::Spacing, value, newSpacing, newRegion);
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(7), newRegion.GetSize(0));
  }

  // The resampled values were truncated, a bias of about -0.5 in every value.
  void Resample_Linear_RoundsValues()
  {
    auto dwi = DwiType::New();
    dwi->SetRegions(DwiType::RegionType({{0, 0, 0}}, {{2, 1, 1}}));
    dwi->SetVectorLength(2);
    dwi->Allocate();
    DwiType::PixelType pixel(2);
    pixel[0] = 10;
    pixel[1] = 7;
    dwi->SetPixel({{0, 0, 0}}, pixel);
    pixel[0] = 11;
    dwi->SetPixel({{1, 0, 0}}, pixel);

    auto filter = FilterType::New();
    filter->SetInput(dwi);
    itk::Vector<double, 3> factor;
    factor.Fill(1);
    factor[0] = 2;
    filter->SetSamplingFactor(factor);
    filter->SetInterpolation(FilterType::Interpolate_Linear);
    filter->Update();

    // Output voxel x samples the input at (x + 0.5) / 2 - 0.5: 10 (clamped at the border), 10.25,
    // 10.75 and 11; channel 1 is constant.
    const short expected[4] = {10, 10, 11, 11};
    for (unsigned int x = 0; x < 4; ++x)
    {
      CPPUNIT_ASSERT_EQUAL(expected[x], filter->GetOutput()->GetPixel({{x, 0, 0}})[0]);
      CPPUNIT_ASSERT_EQUAL(short(7), filter->GetOutput()->GetPixel({{x, 0, 0}})[1]);
    }
  }

  void Resample_SetterBeforeInputAndUpstreamFilter_RunsPipeline()
  {
    // The grid may be requested before the input is known, and an upstream filter that has not run
    // yet runs with the resampler.
    auto filter = FilterType::New();
    itk::Vector<double, 3> factor;
    factor.Fill(2.0);
    filter->SetSamplingFactor(factor);
    filter->SetInterpolation(FilterType::Interpolate_Linear);
    auto upstream = itk::CastImageFilter<DwiType, DwiType>::New();
    upstream->SetInput(m_Dwi);
    filter->SetInput(upstream->GetOutput());
    filter->Update();

    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(22), filter->GetOutput()->GetLargestPossibleRegion().GetSize(0));
    CheckChannels(m_Dwi, filter->GetOutput(), FilterType::Interpolate_Linear);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionResampleDwiImageFilter)
