/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionImageCorrectionFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionImageCorrectionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageCorrectionFilterTestSuite);
  MITK_TEST(CorrectDirections_OneAffine_RotatesAllWeightedGradients);
  MITK_TEST(CorrectDirections_OneMatrixNoB0_RotatesFirstVolumeToo);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using FilterType = mitk::DiffusionImageCorrectionFilter;

  mitk::Image::Pointer m_Image;
  FilterType::TransformMatrixType m_Rotation;

  void CheckRotated(PropHelper::GradientDirectionsContainerType::ConstPointer input, const FilterType::TransformMatrixType& transform)
  {
    auto filter = FilterType::New();
    filter->SetImage(m_Image);
    filter->CorrectDirections(transform);

    auto output = PropHelper::GetGradientContainer(m_Image);
    CPPUNIT_ASSERT_EQUAL(input->Size(), output->Size());
    for (unsigned int i = 0; i < input->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (output->ElementAt(i) - m_Rotation * input->ElementAt(i)).two_norm(), 1e-9);
  }

public:
  void setUp() override
  {
    m_Image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
    // 90 degrees about z
    m_Rotation.fill(0.0);
    m_Rotation[0][1] = -1;
    m_Rotation[1][0] = 1;
    m_Rotation[2][2] = 1;
  }

  void tearDown() override
  {
    m_Image = nullptr;
  }

  void CorrectDirections_OneAffine_RotatesAllWeightedGradients()
  {
    // The rotation part of rotation * scaling is the rotation.
    FilterType::TransformMatrixType scaling;
    scaling.set_identity();
    scaling[0][0] = 2;
    CheckRotated(PropHelper::GetGradientContainer(m_Image), m_Rotation * scaling);
  }

  void CorrectDirections_OneMatrixNoB0_RotatesFirstVolumeToo()
  {
    // Volume 0, the only b0 volume of test.dwi, becomes a weighted one: the lowest shell is then
    // b = 3000 and starts with volume 0.
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : PropHelper::GetGradientContainer(m_Image)->CastToSTLConstContainer())
      gradients->push_back(g);
    gradients->ElementAt(0).fill(0.0);
    gradients->ElementAt(0)[0] = 1.0;
    m_Image->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    PropHelper::SetGradientContainer(m_Image, gradients);
    PropHelper::InitializeImage(m_Image);
    CPPUNIT_ASSERT_EQUAL(size_t(1), PropHelper::GetBValueMap(m_Image).size());

    CheckRotated(gradients.GetPointer(), m_Rotation);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageCorrectionFilter)
