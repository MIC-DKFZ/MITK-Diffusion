/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkADCAverageFunctor.h>
#include <itkBiExpFitFunctor.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkKurtosisFitFunctor.h>
#include <itkRadialMultishellToSingleshellImageFilter.h>

#include <cmath>
#include <limits>

/**
 * The radial multi-shell to single-shell filter with its three functors, as MitkMultishellMethods
 * runs it ("-b mean"), against the references in MITK-Data.
 */
class mitkDiffusionReconstructionMultishellTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionMultishellTestSuite);
  MITK_TEST(AdcAverage_MatchesReference);
  MITK_TEST(KurtosisFit_MatchesReference);
  MITK_TEST(BiExpFit_MatchesReference);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FilterType = itk::RadialMultishellToSingleshellImageFilter<short, short>;

  DwiType::Pointer m_Dwi;
  PropHelper::GradientDirectionsContainerType::Pointer m_Gradients;
  double m_ReferenceBValue = 0;
  PropHelper::BValueMapType m_BValueMap;
  vnl_vector<double> m_BValues;

  mitk::Image::Pointer LoadDwi(const std::string &path)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath(path), &functor);
  }

public:
  void setUp() override
  {
    auto image = LoadDwi("DiffusionImaging/ImageReconstruction/test_multi.dwi");
    m_Dwi = DwiType::New();
    mitk::CastToItkImage(image, m_Dwi);

    // b-values rounded to multiples of 1000 as MitkMultishellMethods does it
    // (mitk::DiffusionImageEditing::RoundBValues, which this module cannot link): the smallest
    // shell becomes the reference b-value, gradients are rescaled to the rounded b-values.
    const auto gradients = PropHelper::GetGradientContainer(image);
    const double bRef = PropHelper::GetReferenceBValue(image);
    auto rounded = [&](const PropHelper::GradientDirectionType &g) {
      const double b = bRef * g.two_norm() * g.two_norm();
      return static_cast<double>(int((b + 500) / 1000) * 1000);
    };
    m_ReferenceBValue = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
    {
      const double b = rounded(gradients->ElementAt(i));
      if (b > 0 && b < m_ReferenceBValue)
        m_ReferenceBValue = b;
    }
    m_Gradients = PropHelper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
    {
      PropHelper::GradientDirectionType g(0.0);
      const double b = rounded(gradients->ElementAt(i));
      if (b > 0)
      {
        g = gradients->ElementAt(i);
        g.normalize();
        g *= std::sqrt(b / m_ReferenceBValue);
      }
      m_Gradients->push_back(g);
    }
    m_BValueMap = mitk::BValueMapProperty::CreateBValueMap(m_Gradients, m_ReferenceBValue);
    CPPUNIT_ASSERT_EQUAL(std::size_t(4), m_BValueMap.size());

    m_BValues.set_size(m_BValueMap.size() - 1);
    unsigned int i = 0;
    for (const auto &shell : m_BValueMap)
    {
      if (shell.first != 0)
        m_BValues[i++] = shell.first;
    }
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
    m_Gradients = nullptr;
  }

  void Check(itk::DWIVoxelFunctor *functor, const std::string &reference)
  {
    auto filter = FilterType::New();
    filter->SetInput(m_Dwi);
    filter->SetOriginalGradientDirections(m_Gradients.GetPointer());
    filter->SetOriginalBValueMap(m_BValueMap);
    filter->SetOriginalBValue(m_ReferenceBValue);
    filter->SetFunctor(functor);
    filter->Update();

    auto expectedImage = LoadDwi("DiffusionImaging/ImageReconstruction/" + reference);
    auto expected = DwiType::New();
    mitk::CastToItkImage(expectedImage, expected);
    auto output = filter->GetOutput();
    CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfComponentsPerPixel(), output->GetNumberOfComponentsPerPixel());
    CPPUNIT_ASSERT(expected->GetLargestPossibleRegion() == output->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> itExpected(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> itOutput(output, expected->GetLargestPossibleRegion());
    for (; !itExpected.IsAtEnd(); ++itExpected, ++itOutput)
      CPPUNIT_ASSERT(itExpected.Get() == itOutput.Get());

    const auto expectedGradients = PropHelper::GetGradientContainer(expectedImage);
    const auto gradients = filter->GetTargetGradientDirections();
    CPPUNIT_ASSERT_EQUAL(expectedGradients->Size(), gradients->Size());
    for (unsigned int i = 0; i < gradients->Size(); ++i)
    {
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedGradients->ElementAt(i)[c], gradients->ElementAt(i)[c], 1e-6);
    }
  }

  void AdcAverage_MatchesReference()
  {
    auto functor = itk::ADCAverageFunctor::New();
    functor->setListOfBValues(m_BValues);
    functor->setTargetBValue(m_BValues.mean());
    Check(functor, "test_multi_ADC.dwi");
  }

  void KurtosisFit_MatchesReference()
  {
    auto functor = itk::KurtosisFitFunctor::New();
    functor->setListOfBValues(m_BValues);
    functor->setTargetBValue(m_BValues.mean());
    Check(functor, "test_multi_AKC.dwi");
  }

  void BiExpFit_MatchesReference()
  {
    auto functor = itk::BiExpFitFunctor::New();
    functor->setListOfBValues(m_BValues);
    functor->setTargetBValue(m_BValues.mean());
    Check(functor, "test_multi_BiExp.dwi");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionMultishell)
