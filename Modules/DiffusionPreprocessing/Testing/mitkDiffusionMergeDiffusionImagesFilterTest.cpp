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

#include <itkImageDuplicator.h>
#include <itkImageRegionConstIterator.h>
#include <itkMergeDiffusionImagesFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionMergeDiffusionImagesFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionMergeDiffusionImagesFilterTestSuite);
  MITK_TEST(Merge_SameGeometry_ConcatenatesVolumes);
  MITK_TEST(Merge_DifferentSpacing_Throws);
  MITK_TEST(Merge_DifferentOrigin_Throws);
  MITK_TEST(Merge_DifferentDirection_Throws);
  MITK_TEST(Merge_MoreGradientsThanVolumes_Throws);
  MITK_TEST(Merge_ZeroReferenceBValue_Throws);
  MITK_TEST(Merge_BValueRoundingToZero_IsB0);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::MergeDiffusionImagesFilter<short>;

  DwiType::Pointer m_Dwi;
  DwiType::Pointer m_Copy;
  PropHelper::GradientDirectionsContainerType::ConstPointer m_Gradients;
  double m_BValue = 0;

  FilterType::Pointer CreateFilter(PropHelper::GradientDirectionsContainerType::ConstPointer secondGradients, double secondBValue)
  {
    auto filter = FilterType::New();
    filter->SetImageVolumes({m_Dwi, m_Copy});
    filter->SetGradientLists({m_Gradients, secondGradients});
    filter->SetBValues({m_BValue, secondBValue});
    return filter;
  }

public:
  void setUp() override
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
    m_Dwi = PropHelper::GetItkVectorImage(image);
    m_Gradients = PropHelper::GetGradientContainer(image);
    m_BValue = PropHelper::GetReferenceBValue(image);

    auto duplicator = itk::ImageDuplicator<DwiType>::New();
    duplicator->SetInputImage(m_Dwi);
    duplicator->Update();
    m_Copy = duplicator->GetOutput();
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
    m_Copy = nullptr;
    m_Gradients = nullptr;
  }

  void Merge_SameGeometry_ConcatenatesVolumes()
  {
    auto filter = CreateFilter(m_Gradients, m_BValue);
    filter->Update();

    const unsigned int n = m_Dwi->GetVectorLength();
    auto output = filter->GetOutput();
    CPPUNIT_ASSERT_EQUAL(2 * n, output->GetVectorLength());
    CPPUNIT_ASSERT_EQUAL(m_BValue, filter->GetB_Value());
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(2 * n), static_cast<size_t>(filter->GetOutputGradients()->Size()));
    for (unsigned int c = 0; c < n; ++c)
    {
      // Equal reference b-values keep the gradients, b-value included in their length.
      const auto expected = m_Gradients->ElementAt(c);
      for (unsigned int copy = 0; copy < 2; ++copy)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (filter->GetOutputGradients()->ElementAt(c + copy * n) - expected).two_norm(), 1e-9);
    }

    itk::ImageRegionConstIterator<DwiType> in(m_Dwi, m_Dwi->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> out(output, output->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
      for (unsigned int c = 0; c < n; ++c)
      {
        CPPUNIT_ASSERT_EQUAL(in.Get()[c], out.Get()[c]);
        CPPUNIT_ASSERT_EQUAL(in.Get()[c], out.Get()[c + n]);
      }
  }

  void Merge_DifferentSpacing_Throws()
  {
    auto spacing = m_Copy->GetSpacing();
    spacing[2] *= 2;
    m_Copy->SetSpacing(spacing);
    CPPUNIT_ASSERT_THROW(CreateFilter(m_Gradients, m_BValue)->Update(), itk::ExceptionObject);
  }

  void Merge_DifferentOrigin_Throws()
  {
    auto origin = m_Copy->GetOrigin();
    origin[0] += 1;
    m_Copy->SetOrigin(origin);
    CPPUNIT_ASSERT_THROW(CreateFilter(m_Gradients, m_BValue)->Update(), itk::ExceptionObject);
  }

  void Merge_DifferentDirection_Throws()
  {
    DwiType::DirectionType direction;
    direction.SetIdentity();
    direction[1][1] = -1;
    m_Copy->SetDirection(direction);
    CPPUNIT_ASSERT_THROW(CreateFilter(m_Gradients, m_BValue)->Update(), itk::ExceptionObject);
  }

  void Merge_MoreGradientsThanVolumes_Throws()
  {
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : m_Gradients->CastToSTLConstContainer())
      gradients->push_back(g);
    gradients->push_back(m_Gradients->ElementAt(1));
    CPPUNIT_ASSERT_THROW(CreateFilter(gradients.GetPointer(), m_BValue)->Update(), itk::ExceptionObject);
  }

  void Merge_ZeroReferenceBValue_Throws()
  {
    auto filter = FilterType::New();
    filter->SetImageVolumes({m_Dwi, m_Copy});
    filter->SetGradientLists({m_Gradients, m_Gradients});
    filter->SetBValues({0, 0});
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Merge_BValueRoundingToZero_IsB0()
  {
    // Volume 0 of the second image gets b = 3000 * 0.01^2 = 0.3, which rounds to 0 as in the
    // b-value map: its merged gradient is the zero vector, not a rescaled short one.
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : m_Gradients->CastToSTLConstContainer())
      gradients->push_back(g);
    gradients->ElementAt(0)[0] = 0.01;
    auto filter = CreateFilter(gradients.GetPointer(), m_BValue);
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(0.0, filter->GetOutputGradients()->ElementAt(m_Dwi->GetVectorLength()).two_norm());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionMergeDiffusionImagesFilter)
