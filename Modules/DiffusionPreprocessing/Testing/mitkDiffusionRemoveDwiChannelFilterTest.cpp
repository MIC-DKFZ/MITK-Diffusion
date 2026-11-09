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

#include <itkImageRegionConstIterator.h>
#include <itkRemoveDwiChannelFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionRemoveDwiChannelFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRemoveDwiChannelFilterTestSuite);
  MITK_TEST(Remove_UnsortedAndDuplicateIndices_RemovesEachOnce);
  MITK_TEST(Remove_IndexOutOfRange_Throws);
  MITK_TEST(Remove_AllChannels_Throws);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::RemoveDwiChannelFilter<short>;

  DwiType::Pointer m_Dwi;
  PropHelper::GradientDirectionsContainerType::ConstPointer m_Gradients;

  FilterType::Pointer CreateFilter(const std::vector<unsigned int>& indices)
  {
    auto filter = FilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_Gradients);
    filter->SetChannelIndices(indices);
    return filter;
  }

public:
  void setUp() override
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
    m_Dwi = PropHelper::GetItkVectorImage(image);
    m_Gradients = PropHelper::GetGradientContainer(image);
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
    m_Gradients = nullptr;
  }

  void Remove_UnsortedAndDuplicateIndices_RemovesEachOnce()
  {
    auto filter = CreateFilter({5, 3, 5});
    filter->Update();

    std::vector<unsigned int> kept;
    for (unsigned int i = 0; i < m_Dwi->GetVectorLength(); ++i)
      if (i != 3 && i != 5)
        kept.push_back(i);

    auto output = filter->GetOutput();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(kept.size()), output->GetVectorLength());
    CPPUNIT_ASSERT_EQUAL(kept.size(), static_cast<size_t>(filter->GetNewDirections()->Size()));
    for (unsigned int c = 0; c < kept.size(); ++c)
      CPPUNIT_ASSERT(filter->GetNewDirections()->ElementAt(c) == m_Gradients->ElementAt(kept[c]));

    itk::ImageRegionConstIterator<DwiType> in(m_Dwi, m_Dwi->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> out(output, output->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
      for (unsigned int c = 0; c < kept.size(); ++c)
        CPPUNIT_ASSERT_EQUAL(in.Get()[kept[c]], out.Get()[c]);
  }

  void Remove_IndexOutOfRange_Throws()
  {
    auto filter = CreateFilter({m_Dwi->GetVectorLength()});
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Remove_AllChannels_Throws()
  {
    std::vector<unsigned int> all;
    for (unsigned int i = 0; i < m_Dwi->GetVectorLength(); ++i)
      all.push_back(i);
    auto filter = CreateFilter(all);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRemoveDwiChannelFilter)
