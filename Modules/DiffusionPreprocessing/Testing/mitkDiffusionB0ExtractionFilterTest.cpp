/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkB0ImageExtractionImageFilter.h>
#include <itkB0ImageExtractionToSeparateImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionB0ExtractionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionB0ExtractionFilterTestSuite);
  MITK_TEST(Average_OneB0_CopiesIt);
  MITK_TEST(Average_NoB0_Throws);
  MITK_TEST(Separate_OneB0_CopiesIt);
  MITK_TEST(Separate_NoB0_Throws);
  MITK_TEST(Average_GradientCountMismatch_Throws);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using AverageFilterType = itk::B0ImageExtractionImageFilter<short, short>;
  using SeparateFilterType = itk::B0ImageExtractionToSeparateImageFilter<short, short>;

  DwiType::Pointer m_Dwi;
  PropHelper::GradientDirectionsContainerType::ConstPointer m_Gradients;
  PropHelper::GradientDirectionsContainerType::Pointer m_GradientsWithoutB0;

public:
  void setUp() override
  {
    // test.dwi: volume 0 is the only b0 volume.
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
    m_Dwi = PropHelper::GetItkVectorImage(image);
    m_Gradients = PropHelper::GetGradientContainer(image);
    m_GradientsWithoutB0 = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : m_Gradients->CastToSTLConstContainer())
      m_GradientsWithoutB0->push_back(g);
    m_GradientsWithoutB0->ElementAt(0)[0] = 1.0;
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
    m_Gradients = nullptr;
    m_GradientsWithoutB0 = nullptr;
  }

  void Average_OneB0_CopiesIt()
  {
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_Gradients);
    filter->Update();

    itk::ImageRegionConstIteratorWithIndex<AverageFilterType::OutputImageType> it(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      CPPUNIT_ASSERT_EQUAL(m_Dwi->GetPixel(it.GetIndex())[0], it.Get());
  }

  void Average_NoB0_Throws()
  {
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_GradientsWithoutB0.GetPointer());
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Separate_OneB0_CopiesIt()
  {
    auto filter = SeparateFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_Gradients);
    filter->Update();

    auto output = filter->GetOutput();
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(1), output->GetLargestPossibleRegion().GetSize(3));
    itk::ImageRegionConstIteratorWithIndex<DwiType> it(m_Dwi, m_Dwi->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      SeparateFilterType::OutputImageType::IndexType index;
      for (unsigned int d = 0; d < 3; ++d)
        index[d] = it.GetIndex()[d];
      index[3] = 0;
      CPPUNIT_ASSERT_EQUAL(it.Get()[0], output->GetPixel(index));
    }
  }

  void Separate_NoB0_Throws()
  {
    auto filter = SeparateFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_GradientsWithoutB0.GetPointer());
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Average_GradientCountMismatch_Throws()
  {
    m_GradientsWithoutB0->ElementAt(0).fill(0.0);
    m_GradientsWithoutB0->push_back(m_Gradients->ElementAt(1));
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_GradientsWithoutB0.GetPointer());
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionB0ExtractionFilter)
