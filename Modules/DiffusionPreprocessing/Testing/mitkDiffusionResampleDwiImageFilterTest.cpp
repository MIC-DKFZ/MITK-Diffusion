/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkBSplineInterpolateImageFunction.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleDwiImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionResampleDwiImageFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionResampleDwiImageFilterTestSuite);
  MITK_TEST(Resample_IdentityDirection_KeepsCorner);
  MITK_TEST(Resample_FlippedAxis_KeepsCorner);
  MITK_TEST(Resample_EachSetterAndInterpolation_ResamplesEveryChannelAlike);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::ResampleDwiImageFilter<short>;
  using ChannelType = itk::Image<short, 3>;

  /** The filter's grid and interpolators applied to one channel with itk::ResampleImageFilter. */
  static ChannelType::Pointer ResampleChannel(const DwiType* input, unsigned int channel, const DwiType* output, FilterType::Interpolation interpolation)
  {
    auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiType, ChannelType>::New();
    extractor->SetInput(input);
    extractor->SetIndex(channel);
    auto resampler = itk::ResampleImageFilter<ChannelType, ChannelType>::New();
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
        CPPUNIT_ASSERT_EQUAL(expected->GetPixel(it.GetIndex()), it.Get()[c]);
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
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionResampleDwiImageFilter)
