/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkElectrostaticRepulsionDiffusionGradientReductionFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

class mitkDiffusionGradientReductionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionGradientReductionFilterTestSuite);
  MITK_TEST(Reduce_SameSeed_SameSelection);
  MITK_TEST(Reduce_UseFirstN_KeepsFirstVolumesOfEachShell);
  MITK_TEST(Reduce_DefaultSeed_KnownSelection);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::ElectrostaticRepulsionDiffusionGradientReductionFilter<short, short>;

  mitk::Image::Pointer m_Image;

  /** Keeps the b0 volume and 15 of the 30 directions of the b = 3000 shell. */
  FilterType::IndicesVector Reduce(bool useFirstN)
  {
    FilterType::BValueMap shells;
    for (const auto& shell : PropHelper::GetBValueMap(m_Image))
      shells[shell.first] = shell.second;
    CPPUNIT_ASSERT_EQUAL(size_t(2), shells.size());

    auto filter = FilterType::New();
    filter->SetInput(PropHelper::GetItkVectorImage(m_Image));
    filter->SetOriginalGradientDirections(PropHelper::GetGradientContainer(m_Image));
    filter->SetShellSelectionBValueMap(shells);
    filter->SetNumGradientDirections({1, 15});
    filter->SetUseFirstN(useFirstN);
    filter->Update();

    FilterType::IndicesVector selected;
    for (const auto& g : filter->GetGradientDirections()->CastToSTLConstContainer())
    {
      // Map each output gradient back to its input volume; the directions of the shell differ.
      const auto input = PropHelper::GetGradientContainer(m_Image);
      for (unsigned int i = 0; i < input->Size(); ++i)
        if (input->ElementAt(i) == g)
        {
          selected.push_back(i);
          break;
        }
    }
    CPPUNIT_ASSERT_EQUAL(size_t(16), selected.size());
    CPPUNIT_ASSERT_EQUAL(16u, filter->GetOutput()->GetVectorLength());
    return selected;
  }

public:
  void setUp() override
  {
    m_Image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
  }

  void tearDown() override
  {
    m_Image = nullptr;
  }

  void Reduce_SameSeed_SameSelection()
  {
    const auto first = Reduce(false);
    const auto second = Reduce(false);
    CPPUNIT_ASSERT(first == second);
  }

  void Reduce_UseFirstN_KeepsFirstVolumesOfEachShell()
  {
    FilterType::IndicesVector expected;
    for (unsigned int i = 0; i < 16; ++i)
      expected.push_back(i);
    CPPUNIT_ASSERT(Reduce(true) == expected);
  }

  void Reduce_DefaultSeed_KnownSelection()
  {
    // The selection of the random search with seed 0, as computed before the filter was
    // restructured; a faster cost computation must keep it.
    const FilterType::IndicesVector expected = {0, 2, 5, 7, 10, 13, 15, 17, 20, 23, 24, 25, 26, 27, 28, 29};
    CPPUNIT_ASSERT(Reduce(false) == expected);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionGradientReductionFilter)
