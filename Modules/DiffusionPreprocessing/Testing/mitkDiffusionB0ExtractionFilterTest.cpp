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

#include <itkB0ImageExtractionImageFilter.h>
#include <itkB0ImageExtractionToSeparateImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

#include <cmath>

class mitkDiffusionB0ExtractionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionB0ExtractionFilterTestSuite);
  MITK_TEST(Average_OneB0_CopiesIt);
  MITK_TEST(Average_NoB0_Throws);
  MITK_TEST(Separate_OneB0_CopiesIt);
  MITK_TEST(Separate_NoB0_Throws);
  MITK_TEST(Average_GradientCountMismatch_Throws);
  MITK_TEST(Average_SeveralB0_IsTheirMean);
  MITK_TEST(Separate_SeveralB0_KeepsEachInOrder);
  MITK_TEST(Average_BValueRoundingToZero_IsB0);
  MITK_TEST(Separate_BValueRoundingToZero_IsB0);
  MITK_TEST(Average_NoReferenceBValue_Throws);
  MITK_TEST(Average_IntegerMean_Exact);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using AverageFilterType = itk::B0ImageExtractionImageFilter<short, short>;
  using SeparateFilterType = itk::B0ImageExtractionToSeparateImageFilter<short, short>;

  DwiType::Pointer m_Dwi;
  PropHelper::GradientDirectionsContainerType::ConstPointer m_Gradients;
  PropHelper::GradientDirectionsContainerType::Pointer m_GradientsWithoutB0;
  double m_BValue = 0;

  /**
   * Files with several b0 volumes: test_multi.dwi (9 b0 among 99 volumes, identity frame) and
   * InputDiffusionWeightedImage.dwi (b0 at 0, 10, 31 and 52, non-zero origin, NRRD frame field).
   */
  static std::vector<std::string> MultiB0Files()
  {
    return { "ImageReconstruction/test_multi.dwi", "FWE/InputDiffusionWeightedImage.dwi" };
  }

  /** Gradients with volume 0 (b0) shortened to b = 3000 * 0.01^2 = 0.3, which rounds to 0. */
  PropHelper::GradientDirectionsContainerType::Pointer GradientsWithNearlyZeroB0() const
  {
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : m_Gradients->CastToSTLConstContainer())
      gradients->push_back(g);
    gradients->ElementAt(0)[0] = 0.01;
    return gradients;
  }

  static std::vector<unsigned int> ZeroGradientIndices(const PropHelper::GradientDirectionsContainerType* gradients)
  {
    std::vector<unsigned int> indices;
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      if (gradients->ElementAt(i).two_norm() == 0)
        indices.push_back(i);
    return indices;
  }

public:
  void setUp() override
  {
    // test.dwi: volume 0 is the only b0 volume.
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
    m_Dwi = PropHelper::GetItkVectorImage(image);
    m_Gradients = PropHelper::GetGradientContainer(image);
    m_BValue = PropHelper::GetReferenceBValue(image);
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
    filter->SetReferenceBValue(m_BValue);
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
    filter->SetReferenceBValue(m_BValue);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Separate_OneB0_CopiesIt()
  {
    auto filter = SeparateFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_Gradients);
    filter->SetReferenceBValue(m_BValue);
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
    filter->SetReferenceBValue(m_BValue);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Average_GradientCountMismatch_Throws()
  {
    m_GradientsWithoutB0->ElementAt(0).fill(0.0);
    m_GradientsWithoutB0->push_back(m_Gradients->ElementAt(1));
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_GradientsWithoutB0.GetPointer());
    filter->SetReferenceBValue(m_BValue);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Average_SeveralB0_IsTheirMean()
  {
    for (const auto& file : MultiB0Files())
    {
      auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(file));
      auto dwi = PropHelper::GetItkVectorImage(image);
      auto gradients = PropHelper::GetGradientContainer(image);
      const auto b0 = ZeroGradientIndices(gradients);
      CPPUNIT_ASSERT(b0.size() > 1);

      auto filter = AverageFilterType::New();
      filter->SetInput(dwi);
      filter->SetDirections(gradients);
      filter->SetReferenceBValue(PropHelper::GetReferenceBValue(image));
      filter->Update();
      auto output = filter->GetOutput();
      CPPUNIT_ASSERT(output->GetOrigin() == dwi->GetOrigin());
      CPPUNIT_ASSERT(output->GetSpacing() == dwi->GetSpacing());
      CPPUNIT_ASSERT(output->GetDirection() == dwi->GetDirection());

      // The mean, rounded to short.
      itk::ImageRegionConstIteratorWithIndex<DwiType> it(dwi, dwi->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
      {
        double sum = 0;
        for (auto i : b0)
          sum += it.Get()[i];
        CPPUNIT_ASSERT_EQUAL(static_cast<short>(std::round(sum / b0.size())), output->GetPixel(it.GetIndex()));
      }
    }
  }

  void Separate_SeveralB0_KeepsEachInOrder()
  {
    for (const auto& file : MultiB0Files())
    {
      auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(file));
      auto dwi = PropHelper::GetItkVectorImage(image);
      auto gradients = PropHelper::GetGradientContainer(image);
      const auto b0 = ZeroGradientIndices(gradients);

      auto filter = SeparateFilterType::New();
      filter->SetInput(dwi);
      filter->SetDirections(gradients);
      filter->SetReferenceBValue(PropHelper::GetReferenceBValue(image));
      filter->Update();
      auto output = filter->GetOutput();
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(b0.size()), output->GetLargestPossibleRegion().GetSize(3));
      for (unsigned int d = 0; d < 3; ++d)
      {
        CPPUNIT_ASSERT_EQUAL(dwi->GetOrigin()[d], output->GetOrigin()[d]);
        CPPUNIT_ASSERT_EQUAL(dwi->GetSpacing()[d], output->GetSpacing()[d]);
      }

      itk::ImageRegionConstIteratorWithIndex<DwiType> it(dwi, dwi->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
      {
        SeparateFilterType::OutputImageType::IndexType index;
        for (unsigned int d = 0; d < 3; ++d)
          index[d] = it.GetIndex()[d];
        for (unsigned int t = 0; t < b0.size(); ++t)
        {
          index[3] = t;
          CPPUNIT_ASSERT_EQUAL(it.Get()[b0[t]], output->GetPixel(index));
        }
      }
    }
  }

  void Average_BValueRoundingToZero_IsB0()
  {
    // As in the b-value map: a b-value that rounds to 0 marks a b0 volume, not only a zero vector.
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(GradientsWithNearlyZeroB0().GetPointer());
    filter->SetReferenceBValue(m_BValue);
    filter->Update();

    itk::ImageRegionConstIteratorWithIndex<AverageFilterType::OutputImageType> it(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      CPPUNIT_ASSERT_EQUAL(m_Dwi->GetPixel(it.GetIndex())[0], it.Get());
  }

  void Separate_BValueRoundingToZero_IsB0()
  {
    auto filter = SeparateFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(GradientsWithNearlyZeroB0().GetPointer());
    filter->SetReferenceBValue(m_BValue);
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(1), filter->GetOutput()->GetLargestPossibleRegion().GetSize(3));
  }

  // The mean was accumulated as x/9 in float and truncated: nine times 96 gave 95, and a mean of
  // 100.56 gave 100.
  void Average_IntegerMean_Exact()
  {
    auto dwi = DwiType::New();
    dwi->SetRegions(DwiType::RegionType({{0, 0, 0}}, {{2, 1, 1}}));
    dwi->SetVectorLength(10);
    dwi->Allocate();
    DwiType::PixelType pixel(10);
    pixel.Fill(96);
    pixel[9] = 500;
    dwi->SetPixel({{0, 0, 0}}, pixel);
    pixel.Fill(100);
    pixel[0] = 101; // mean 100.11 -> 100
    pixel[1] = 104; // mean 100.56 -> 101
    dwi->SetPixel({{1, 0, 0}}, pixel);

    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < 10; ++i)
      gradients->push_back(PropHelper::GradientDirectionType(i < 9 ? 0.0 : 1.0, 0, 0));

    auto filter = AverageFilterType::New();
    filter->SetInput(dwi);
    filter->SetDirections(gradients);
    filter->SetReferenceBValue(1000);
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(short(96), filter->GetOutput()->GetPixel({{0, 0, 0}}));
    CPPUNIT_ASSERT_EQUAL(short(101), filter->GetOutput()->GetPixel({{1, 0, 0}}));
  }

  void Average_NoReferenceBValue_Throws()
  {
    auto filter = AverageFilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetDirections(m_Gradients);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionB0ExtractionFilter)
