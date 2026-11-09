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
#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>

#include <cmath>

/**
 * The volume and gradient edits of the preprocessing view, on test_multi.dwi (99 volumes, 9 b0
 * volumes at 0, 1, 2, 33, ..., shells 1000, 2000 and 3000 of 30 directions each).
 */
class mitkDiffusionVolumeEditingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionVolumeEditingTestSuite);
  MITK_TEST(ExtractVolume_EachIndex_CopiesTheChannel);
  MITK_TEST(RemoveVolumes_TwoVolumes_KeepsTheOthersInOrder);
  MITK_TEST(ExtractB0_SeveralB0_AverageAndSeparateVolumes);
  MITK_TEST(FlipGradients_XAndZ_NegatesComponents);
  MITK_TEST(PermuteGradientComponents_ZXY_ReordersComponents);
  MITK_TEST(PermuteGradientComponents_RepeatedAxis_Throws);
  MITK_TEST(MirrorGradientsToHalfSphere_NegativeX_Negated);
  MITK_TEST(ReduceGradients_UseFirstN_KeepsFirstDirectionsOfSelectedShells);
  MITK_TEST(ReduceGradients_UnknownShell_Throws);
  MITK_TEST(AverageRepeatedGradients_DuplicatedVolumes_AveragedBack);
  MITK_TEST(Edits_NotADwi_Throw);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;
  using GradientsType = PropHelper::GradientDirectionsContainerType;

  mitk::Image::Pointer m_Dwi;

  static VolumeType::Pointer ToItk(const mitk::Image* image)
  {
    VolumeType::Pointer itkImage;
    mitk::CastToItkImage(image, itkImage);
    return itkImage;
  }

  /** Channel c of the output equals channel inputChannel of the input in every voxel. */
  static void CheckChannel(const DwiType* input, unsigned int inputChannel, const DwiType* output, unsigned int c)
  {
    itk::ImageRegionConstIterator<DwiType> in(input, input->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> out(output, output->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
      CPPUNIT_ASSERT_EQUAL(in.Get()[inputChannel], out.Get()[c]);
  }

  static void CheckVolume(const DwiType* dwi, unsigned int channel, const VolumeType* volume)
  {
    itk::ImageRegionConstIterator<DwiType> d(dwi, dwi->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<VolumeType> v(volume, volume->GetLargestPossibleRegion());
    for (; !d.IsAtEnd(); ++d, ++v)
      CPPUNIT_ASSERT_EQUAL(d.Get()[channel], v.Get());
  }

  /** Each modified gradient of the result is the changed input gradient; data and b-value stay. */
  void CheckGradientEdit(mitk::Image* result, const std::function<PropHelper::GradientDirectionType(const PropHelper::GradientDirectionType&)>& change)
  {
    auto input = PropHelper::GetGradientContainer(m_Dwi);
    auto output = PropHelper::GetGradientContainer(result);
    CPPUNIT_ASSERT_EQUAL(input->Size(), output->Size());
    for (unsigned int i = 0; i < input->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (change(input->ElementAt(i)) - output->ElementAt(i)).two_norm(), 1e-12);
    CPPUNIT_ASSERT_EQUAL(PropHelper::GetReferenceBValue(m_Dwi), PropHelper::GetReferenceBValue(result));

    // The originals are derived again, so a later InitializeImage keeps the gradients.
    auto copy = result->Clone();
    PropHelper::InitializeImage(copy);
    for (unsigned int i = 0; i < output->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (PropHelper::GetGradientContainer(copy)->ElementAt(i) - output->ElementAt(i)).two_norm(), 1e-12);

    auto in = PropHelper::GetItkVectorImage(m_Dwi);
    auto out = PropHelper::GetItkVectorImage(result);
    itk::ImageRegionConstIterator<DwiType> itIn(in, in->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> itOut(out, out->GetLargestPossibleRegion());
    for (; !itIn.IsAtEnd(); ++itIn, ++itOut)
      CPPUNIT_ASSERT(itIn.Get() == itOut.Get());
  }

public:
  void setUp() override
  {
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_multi.dwi"));
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void ExtractVolume_EachIndex_CopiesTheChannel()
  {
    // The volumes are kept while the DWI is accessed again: they must not hold an access to it.
    const unsigned int n = PropHelper::GetGradientContainer(m_Dwi)->Size();
    const std::vector<unsigned int> indices = {0, 1, n - 1};
    std::vector<mitk::Image::Pointer> volumes;
    for (auto i : indices)
      volumes.push_back(mitk::DiffusionImageEditing::ExtractVolume(m_Dwi, i));
    auto dwi = PropHelper::GetItkVectorImage(m_Dwi);
    for (unsigned int k = 0; k < indices.size(); ++k)
    {
      CPPUNIT_ASSERT(mitk::Equal(*m_Dwi->GetGeometry(), *volumes[k]->GetGeometry(), mitk::eps, mitk::eps, false));
      CheckVolume(dwi, indices[k], ToItk(volumes[k]));
    }
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::ExtractVolume(m_Dwi, n), itk::ExceptionObject);
  }

  void RemoveVolumes_TwoVolumes_KeepsTheOthersInOrder()
  {
    auto result = mitk::DiffusionImageEditing::RemoveVolumes(m_Dwi, {5, 1});
    auto input = PropHelper::GetGradientContainer(m_Dwi);
    auto output = PropHelper::GetGradientContainer(result);
    CPPUNIT_ASSERT_EQUAL(input->Size() - 2, output->Size());

    auto in = PropHelper::GetItkVectorImage(m_Dwi);
    auto out = PropHelper::GetItkVectorImage(result);
    unsigned int c = 0;
    for (unsigned int i = 0; i < input->Size(); ++i)
    {
      if (i == 1 || i == 5)
        continue;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (input->ElementAt(i) - output->ElementAt(c)).two_norm(), 1e-12);
      if (i % 10 == 0)
        CheckChannel(in, i, out, c);
      ++c;
    }
    // Volume 1 is one of the nine b0 volumes, volume 5 a weighted one.
    CPPUNIT_ASSERT_EQUAL(size_t(8), PropHelper::GetBValueMap(result).at(0).size());
  }

  void ExtractB0_SeveralB0_AverageAndSeparateVolumes()
  {
    const auto b0 = PropHelper::GetBValueMap(m_Dwi).at(0);
    auto dwi = PropHelper::GetItkVectorImage(m_Dwi);

    auto average = ToItk(mitk::DiffusionImageEditing::ExtractB0(m_Dwi));
    itk::ImageRegionConstIteratorWithIndex<DwiType> it(dwi, dwi->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      // The mean, rounded to short.
      double sum = 0;
      for (auto i : b0)
        sum += it.Get()[i];
      CPPUNIT_ASSERT_EQUAL(static_cast<short>(std::round(sum / b0.size())), average->GetPixel(it.GetIndex()));
    }

    auto separate = mitk::DiffusionImageEditing::ExtractB0Volumes(m_Dwi);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(b0.size()), separate->GetTimeSteps());
    mitk::ImagePixelReadAccessor<short, 4> volumes(separate);
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      itk::Index<4> index;
      for (unsigned int d = 0; d < 3; ++d)
        index[d] = it.GetIndex()[d];
      for (unsigned int t = 0; t < b0.size(); ++t)
      {
        index[3] = t;
        CPPUNIT_ASSERT_EQUAL(it.Get()[b0[t]], volumes.GetPixelByIndex(index));
      }
    }
  }

  void FlipGradients_XAndZ_NegatesComponents()
  {
    auto result = mitk::DiffusionImageEditing::FlipGradients(m_Dwi, {true, false, true});
    CheckGradientEdit(result, [](const PropHelper::GradientDirectionType& g)
    {
      PropHelper::GradientDirectionType e = g;
      e[0] *= -1;
      e[2] *= -1;
      return e;
    });
  }

  void PermuteGradientComponents_ZXY_ReordersComponents()
  {
    auto result = mitk::DiffusionImageEditing::PermuteGradientComponents(m_Dwi, {2, 0, 1});
    CheckGradientEdit(result, [](const PropHelper::GradientDirectionType& g)
    {
      PropHelper::GradientDirectionType e;
      e[0] = g[2];
      e[1] = g[0];
      e[2] = g[1];
      return e;
    });
  }

  void PermuteGradientComponents_RepeatedAxis_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::PermuteGradientComponents(m_Dwi, {0, 0, 2}), mitk::Exception);
  }

  void MirrorGradientsToHalfSphere_NegativeX_Negated()
  {
    auto result = mitk::DiffusionImageEditing::MirrorGradientsToHalfSphere(m_Dwi);
    CheckGradientEdit(result, [](const PropHelper::GradientDirectionType& g)
    {
      return g[0] < 0 ? PropHelper::GradientDirectionType(-g) : g;
    });
  }

  void ReduceGradients_UseFirstN_KeepsFirstDirectionsOfSelectedShells()
  {
    const auto shells = PropHelper::GetBValueMap(m_Dwi);
    auto result = mitk::DiffusionImageEditing::ReduceGradients(m_Dwi, {{0, 9}, {2000, 10}}, true);

    // The b0 volumes and the first 10 volumes of shell 2000, in input order; shells 1000 and 3000
    // are dropped.
    std::vector<unsigned int> kept(shells.at(0).begin(), shells.at(0).end());
    kept.insert(kept.end(), shells.at(2000).begin(), shells.at(2000).begin() + 10);
    std::sort(kept.begin(), kept.end());

    auto input = PropHelper::GetGradientContainer(m_Dwi);
    auto output = PropHelper::GetGradientContainer(result);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(kept.size()), output->Size());
    auto in = PropHelper::GetItkVectorImage(m_Dwi);
    auto out = PropHelper::GetItkVectorImage(result);
    for (unsigned int c = 0; c < kept.size(); ++c)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (input->ElementAt(kept[c]) - output->ElementAt(c)).two_norm(), 1e-12);
      CheckChannel(in, kept[c], out, c);
    }
    CPPUNIT_ASSERT_EQUAL(size_t(2), PropHelper::GetBValueMap(result).size());
  }

  void ReduceGradients_UnknownShell_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::ReduceGradients(m_Dwi, {{1500, 10}}, true), mitk::Exception);
  }

  void AverageRepeatedGradients_DuplicatedVolumes_AveragedBack()
  {
    // A merge of the DWI with itself repeats every direction; averaging the repetitions gives the
    // input back, its nine b0 volumes averaged into one.
    auto doubled = mitk::DiffusionImageEditing::MergeDwis({m_Dwi, m_Dwi});
    auto result = mitk::DiffusionImageEditing::AverageRepeatedGradients(doubled, 0.001);
    const auto shells = PropHelper::GetBValueMap(result);
    CPPUNIT_ASSERT_EQUAL(size_t(1), shells.at(0).size());
    CPPUNIT_ASSERT_EQUAL(size_t(30), shells.at(1000).size());
    CPPUNIT_ASSERT_EQUAL(size_t(30), shells.at(3000).size());
  }

  void Edits_NotADwi_Throw()
  {
    auto volume = mitk::DiffusionImageEditing::ExtractVolume(m_Dwi, 0);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::RemoveVolumes(volume, {0}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::ExtractB0(volume), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::FlipGradients(volume, {true, false, false}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::AverageRepeatedGradients(volume, 0.001), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionVolumeEditing)
