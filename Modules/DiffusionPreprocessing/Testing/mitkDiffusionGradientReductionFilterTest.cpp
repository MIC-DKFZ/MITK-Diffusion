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

#include <itkElectrostaticRepulsionDiffusionGradientReductionFilter.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>

#include <algorithm>
#include <map>
#include <vector>

class mitkDiffusionGradientReductionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionGradientReductionFilterTestSuite);
  MITK_TEST(Reduce_SameSeed_SameSelection);
  MITK_TEST(Reduce_UseFirstN_KeepsFirstVolumesOfEachShell);
  MITK_TEST(Reduce_DefaultSeed_KnownSelection);
  MITK_TEST(ReduceMultiShell_CountsAndVolumesPerShell);
  MITK_TEST(ReduceMultiShell_UseFirstN_FirstOfEachShell);
  MITK_TEST(ReduceMultiShell_SameSeed_SameSelection);
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
    m_Image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
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
    // The recorded selection of the random search with seed 0; a faster cost computation must
    // keep it.
    const FilterType::IndicesVector expected = {0, 2, 5, 7, 10, 13, 15, 17, 20, 23, 24, 25, 26, 27, 28, 29};
    CPPUNIT_ASSERT(Reduce(false) == expected);
  }

  /** The multi-shell voxel with b-values rounded to multiples of 1000: 27 b0, 81 per shell 1000, 2000, 3000. */
  static mitk::Image::Pointer LoadMultiShell()
  {
    auto dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"));
    auto rounded = mitk::DiffusionImageEditing::RoundBValues(dwi, 1000);
    const auto &shells = PropHelper::GetBValueMap(rounded);
    CPPUNIT_ASSERT_EQUAL(std::size_t(4), shells.size());
    for (const unsigned int b : {0u, 1000u, 2000u, 3000u})
      CPPUNIT_ASSERT_EQUAL(b == 0 ? std::size_t(27) : std::size_t(81), shells.at(b).size());
    return rounded;
  }

  /**
   * For each output volume the input volume it came from, found by the same gradient within the
   * shell and the same voxel value; b0 volumes, whose gradients are all zero, by their value and
   * order. Every output volume must have one.
   */
  static std::vector<unsigned int> Sources(mitk::Image *input, mitk::Image *output)
  {
    const auto inGradients = PropHelper::GetGradientContainer(input);
    const auto outGradients = PropHelper::GetGradientContainer(output);
    const auto inVoxel = PropHelper::GetItkVectorImage(input)->GetPixel({{0, 0, 0}});
    const auto outVoxel = PropHelper::GetItkVectorImage(output)->GetPixel({{0, 0, 0}});
    std::vector<unsigned int> sources;
    unsigned int nextB0 = 0;
    for (unsigned int j = 0; j < outGradients->Size(); ++j)
    {
      bool found = false;
      const bool b0 = outGradients->ElementAt(j).two_norm() == 0;
      for (unsigned int i = b0 ? nextB0 : 0; i < inGradients->Size() && !found; ++i)
      {
        if (inGradients->ElementAt(i) == outGradients->ElementAt(j) && inVoxel[i] == outVoxel[j])
        {
          sources.push_back(i);
          found = true;
          if (b0)
            nextB0 = i + 1;
        }
      }
      CPPUNIT_ASSERT_MESSAGE("output volume without source", found);
    }
    return sources;
  }

  // Per shell the requested number of volumes, each with its own gradient and value; shells that
  // are not listed are dropped, and the output's b-value map lists what it holds.
  void ReduceMultiShell_CountsAndVolumesPerShell()
  {
    auto dwi = LoadMultiShell();
    const std::map<unsigned int, unsigned int> wanted = {{0, 5}, {1000, 30}, {3000, 81}};
    auto reduced = mitk::DiffusionImageEditing::ReduceGradients(dwi, wanted, false);

    const auto &shells = PropHelper::GetBValueMap(reduced);
    CPPUNIT_ASSERT_EQUAL(wanted.size(), shells.size());
    for (const auto &[b, number] : wanted)
      CPPUNIT_ASSERT_EQUAL(std::size_t(number), shells.at(b).size());
    CPPUNIT_ASSERT_EQUAL(116u, PropHelper::GetItkVectorImage(reduced)->GetVectorLength());

    const auto sources = Sources(dwi, reduced);
    const auto &inputShells = PropHelper::GetBValueMap(dwi);
    for (const auto &[b, indices] : shells)
    {
      std::vector<unsigned int> fromShell;
      for (const auto j : indices)
        fromShell.push_back(sources[j]);
      // Each kept volume is a distinct volume of the same shell of the input.
      std::sort(fromShell.begin(), fromShell.end());
      CPPUNIT_ASSERT(std::adjacent_find(fromShell.begin(), fromShell.end()) == fromShell.end());
      const auto &inputIndices = inputShells.at(b);
      for (const auto i : fromShell)
        CPPUNIT_ASSERT(std::find(inputIndices.begin(), inputIndices.end(), i) != inputIndices.end());
    }
  }

  void ReduceMultiShell_UseFirstN_FirstOfEachShell()
  {
    auto dwi = LoadMultiShell();
    const std::map<unsigned int, unsigned int> wanted = {{0, 3}, {1000, 10}, {2000, 20}, {3000, 5}};
    auto reduced = mitk::DiffusionImageEditing::ReduceGradients(dwi, wanted, true);
    const auto sources = Sources(dwi, reduced);
    const auto &inputShells = PropHelper::GetBValueMap(dwi);
    const auto &shells = PropHelper::GetBValueMap(reduced);
    for (const auto &[b, number] : wanted)
    {
      std::vector<unsigned int> expected(inputShells.at(b).begin(), inputShells.at(b).begin() + number);
      std::vector<unsigned int> actual;
      for (const auto j : shells.at(b))
        actual.push_back(sources[j]);
      std::sort(actual.begin(), actual.end());
      CPPUNIT_ASSERT(expected == actual);
    }
  }

  void ReduceMultiShell_SameSeed_SameSelection()
  {
    auto dwi = LoadMultiShell();
    const std::map<unsigned int, unsigned int> wanted = {{0, 27}, {1000, 40}, {2000, 40}, {3000, 40}};
    auto first = mitk::DiffusionImageEditing::ReduceGradients(dwi, wanted, false);
    auto second = mitk::DiffusionImageEditing::ReduceGradients(dwi, wanted, false);
    CPPUNIT_ASSERT(Sources(dwi, first) == Sources(dwi, second));
    CPPUNIT_ASSERT_EQUAL(147u, PropHelper::GetItkVectorImage(first)->GetVectorLength());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionGradientReductionFilter)
