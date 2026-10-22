/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkResampleDwiImageFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionResampleDwiImageFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionResampleDwiImageFilterTestSuite);
  MITK_TEST(Resample_IdentityDirection_KeepsCorner);
  MITK_TEST(Resample_FlippedAxis_KeepsCorner);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::ResampleDwiImageFilter<short>;

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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionResampleDwiImageFilter)
