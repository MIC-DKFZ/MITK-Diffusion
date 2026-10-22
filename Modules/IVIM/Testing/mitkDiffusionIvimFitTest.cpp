/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.h>

#include "mitkIvimTestHelpers.h"

#include <iomanip>

using namespace mitk::IvimTesting;

class mitkDiffusionIvimFitTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionIvimFitTestSuite);
  MITK_TEST(FitAll_SyntheticVoxels_RecoversParameters);
  MITK_TEST(DStarFix_SyntheticVoxels_RecoversFAndD);
  MITK_TEST(DThenDStar_SyntheticVoxels_RecoversParameters);
  MITK_TEST(LinearDThenF_SyntheticVoxels_RecoversParameters);
  MITK_TEST(Fwe_AllMethods_MatchRecordedValues);
  MITK_TEST(Defaults_ThresholdsAndDStarNotSet_SameAsViewDefaults);
  MITK_TEST(Update_ReferenceBValueNotSet_Throws);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<short, float>;

  struct Parameters
  {
    double f;
    double d;
    double dStar;
  };

  struct Maps
  {
    std::vector<double> f;
    std::vector<double> d;
    std::vector<double> dStar;
  };

  const std::vector<double> m_BValues = {0, 0, 10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 500, 600, 800, 1000};
  const std::vector<Parameters> m_Truth = {{0.1, 0.001, 0.02}, {0.3, 0.0007, 0.05}};

  FilterType::Pointer MakeFilter(const VectorImageType* image,
                                 const GradientContainerType* gradients,
                                 double referenceBValue,
                                 FilterType::IVIM_Method method)
  {
    auto filter = FilterType::New();
    filter->SetInput(image);
    filter->SetGradientDirections(gradients);
    filter->SetBValue(referenceBValue);
    filter->SetMethod(method);
    filter->SetBThres(170);
    filter->SetS0Thres(0);
    filter->SetDStar(0.06);
    filter->SetFitDStar(true);
    // The filter writes per-voxel state into members, so the callers run it on one work unit.
    filter->SetNumberOfWorkUnits(1);
    return filter;
  }

  static Maps GetMaps(FilterType* filter)
  {
    return {Values(filter->GetOutput(0)), Values(filter->GetOutput(1)), Values(filter->GetOutput(2))};
  }

  VectorImageType::Pointer MakeSyntheticImage()
  {
    std::vector<std::vector<double>> signals;
    for (const auto& p : m_Truth)
      signals.push_back(IvimSignal(m_BValues, 10000, p.f, p.d, p.dStar));
    return MakeImage(signals);
  }

  Maps FitSynthetic(FilterType::IVIM_Method method, double dStar = 0.06)
  {
    auto image = MakeSyntheticImage();
    auto gradients = MakeGradients(m_BValues, 1000);
    auto filter = MakeFilter(image, gradients, 1000, method);
    filter->SetDStar(dStar);
    filter->Update();
    return GetMaps(filter);
  }

public:
  void FitAll_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(FilterType::IVIM_FIT_ALL);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.005);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 2e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void DStarFix_SyntheticVoxels_RecoversFAndD()
  {
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      auto maps = FitSynthetic(FilterType::IVIM_DSTAR_FIX, m_Truth[i].dStar);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.005);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 2e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 1e-7);
    }
  }

  void DThenDStar_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(FilterType::IVIM_D_THEN_DSTAR);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      // D and f come from b > 170 only, where the perfusion signal has not fully decayed for
      // D* = 0.02; D* is searched on a grid of step 0.00149.
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.01);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 3e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void LinearDThenF_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(FilterType::IVIM_LINEAR_D_THEN_F);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.01);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 3e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void Fwe_AllMethods_MatchRecordedValues()
  {
    struct Recorded
    {
      FilterType::IVIM_Method method;
      Parameters mean;
      Parameters center;
    };
    // Recorded from the code before the review #13 fixes, b threshold 170, D* 0.06 for the fixed
    // D* method; voxel (5, 5, 5) is in the center of the 11^3 image.
    const std::vector<Recorded> recorded = {
      {FilterType::IVIM_FIT_ALL, {0.4292342052, 0.0005106752415, 0.01820434694}, {0.8115641475, 0.0005265806685, 0.01223030128}},
      {FilterType::IVIM_DSTAR_FIX, {0.2712541107, 0.0007509498357, 0.05999999866}, {0.7812179327, 0.0006338245585, 0.05999999866}},
      {FilterType::IVIM_D_THEN_DSTAR, {0.271253102, 0.0007509525527, 0.0302335529}, {0.7812152505, 0.000633834512, 0.01441000029}},
      {FilterType::IVIM_LINEAR_D_THEN_F, {0.3726078892, 0.0005149761726, 0.02103719706}, {0.7865120173, 0.000659120793, 0.01441000029}}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const std::size_t center = 5 + 11 * 5 + 121 * 5;
    for (const auto& r : recorded)
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue, r.method);
      filter->Update();
      auto maps = GetMaps(filter);
      MITK_INFO << std::setprecision(10) << "method " << r.method << ": mean f " << Mean(maps.f) << ", mean D "
                << Mean(maps.d) << ", mean D* " << Mean(maps.dStar) << "; center f " << maps.f[center] << ", D "
                << maps.d[center] << ", D* " << maps.dStar[center];
      AssertClose(r.mean.f, Mean(maps.f));
      AssertClose(r.mean.d, Mean(maps.d));
      AssertClose(r.mean.dStar, Mean(maps.dStar));
      AssertClose(r.center.f, maps.f[center]);
      AssertClose(r.center.d, maps.d[center]);
      AssertClose(r.center.dStar, maps.dStar[center]);
    }
  }

  void Defaults_ThresholdsAndDStarNotSet_SameAsViewDefaults()
  {
    auto image = MakeSyntheticImage();
    auto gradients = MakeGradients(m_BValues, 1000);
    const std::vector<FilterType::IVIM_Method> methods = {
      FilterType::IVIM_FIT_ALL, FilterType::IVIM_DSTAR_FIX, FilterType::IVIM_D_THEN_DSTAR, FilterType::IVIM_LINEAR_D_THEN_F};
    for (auto method : methods)
    {
      auto explicitFilter = MakeFilter(image, gradients, 1000, method);
      explicitFilter->Update();

      auto defaultFilter = FilterType::New();
      defaultFilter->SetInput(image);
      defaultFilter->SetGradientDirections(gradients.GetPointer());
      defaultFilter->SetBValue(1000);
      defaultFilter->SetMethod(method);
      defaultFilter->SetNumberOfWorkUnits(1);
      defaultFilter->Update();

      auto expected = GetMaps(explicitFilter);
      auto actual = GetMaps(defaultFilter);
      CPPUNIT_ASSERT(expected.f == actual.f);
      CPPUNIT_ASSERT(expected.d == actual.d);
      CPPUNIT_ASSERT(expected.dStar == actual.dStar);
    }
  }

  void Update_ReferenceBValueNotSet_Throws()
  {
    auto filter = FilterType::New();
    filter->SetInput(MakeSyntheticImage());
    filter->SetGradientDirections(MakeGradients(m_BValues, 1000).GetPointer());
    filter->SetMethod(FilterType::IVIM_FIT_ALL);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionIvimFit)
