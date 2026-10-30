/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTestingMacros.h"
#include <mitkFiberBundle.h>
#include <mitkBaseData.h>
#include <itksys/SystemTools.hxx>
#include <mitkTestingConfig.h>
#include <mitkIOUtil.h>
#include <omp.h>
#include <cstring>
#include <itkFitFibersToImageFilter.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkImageToItk.h>
#include <mitkImageCast.h>
#include <itkImageFileWriter.h>
#include <mitkLocaleSwitch.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

class mitkFiberFitTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkFiberFitTestSuite);
  MITK_TEST(Fit1);
  MITK_TEST(Fit2);
  MITK_TEST(Fit3);
  MITK_TEST(Fit4);
  MITK_TEST(Fit5);
  MITK_TEST(Fit6);
  MITK_TEST(RegularizationGradients);
  MITK_TEST(SameAtThreadCounts);
  CPPUNIT_TEST_SUITE_END();

  typedef itk::Image<float, 3> ItkFloatImgType;

private:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/

  typedef itk::FitFibersToImageFilter FitterType;
  FitterType::Pointer fitter;

public:

  mitk::FiberBundle::Pointer LoadFib(std::string fib_name)
  {
    std::vector<mitk::BaseData::Pointer> fibInfile = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath("FiberFit/" + fib_name));
    mitk::BaseData::Pointer baseData = fibInfile.at(0);
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(baseData.GetPointer());
    return fib;
  }

  void setUp() override
  {
    std::vector<mitk::FiberBundle::Pointer> tracts;
    tracts.push_back(LoadFib("Cluster_0.fib"));
    tracts.push_back(LoadFib("Cluster_1.fib"));
    tracts.push_back(LoadFib("Cluster_2.fib"));
    tracts.push_back(LoadFib("Cluster_3.fib"));
    tracts.push_back(LoadFib("Cluster_4.fib"));

    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Peak Image"}, std::vector<std::string>());
    mitk::PeakImage::Pointer peaks = mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("FiberFit/csd_peak_image.nii.gz"), &functor);


    typedef mitk::ImageToItk< mitk::PeakImage::ItkPeakImageType > CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(peaks);
    caster->Update();
    mitk::PeakImage::ItkPeakImageType::Pointer peak_image = caster->GetOutput();

    fitter = FitterType::New();
    fitter->SetPeakImage(peak_image);
    fitter->SetTractograms(tracts);
  }

  void tearDown() override
  {

  }

  void CompareFibs(mitk::FiberBundle::Pointer test, mitk::FiberBundle::Pointer ref, std::string out_name)
  {
    vtkSmartPointer<vtkFloatArray> weights = test->GetFiberWeights();
    vtkSmartPointer<vtkFloatArray> ref_weights = ref->GetFiberWeights();

    CPPUNIT_ASSERT_MESSAGE("Number of weights should be equal", weights->GetNumberOfValues()==ref_weights->GetNumberOfValues());

    for (int i=0; i<weights->GetNumberOfValues(); ++i)
    {
      if (ref_weights->GetValue(i)>0)
      {
        if (fabs( weights->GetValue(i)/ref_weights->GetValue(i)-1 )>0.01)
        {
          mitk::IOUtil::Save(test, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + out_name);
          CPPUNIT_ASSERT_MESSAGE("Weights should be equal", false);
        }
      }
      else if (weights->GetValue(i)>0)
      {
        mitk::IOUtil::Save(test, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + out_name);
        CPPUNIT_ASSERT_MESSAGE("Weights should be equal", false);
      }
    }
  }

  void CompareImages(mitk::PeakImage::ItkPeakImageType::Pointer testImage, std::string name)
  {
    mitk::LocaleSwitch localeSwitch("C");
    typedef mitk::ImageToItk< mitk::PeakImage::ItkPeakImageType > CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("FiberFit/out/" + name)));
    caster->Update();
    mitk::PeakImage::ItkPeakImageType::Pointer refImage = caster->GetOutput();
    CPPUNIT_ASSERT_MESSAGE("Peak image size should equal the size of " + name,
                           testImage->GetLargestPossibleRegion().GetSize() == refImage->GetLargestPossibleRegion().GetSize());

    itk::ImageRegionConstIterator< mitk::PeakImage::ItkPeakImageType > it1(testImage, testImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator< mitk::PeakImage::ItkPeakImageType > it2(refImage, refImage->GetLargestPossibleRegion());

    // peak components (and the residual, over- and underexplained images) are negative as often as positive
    while(!it1.IsAtEnd())
    {
      if (fabs(it2.Get())>0.0001)
      {
        if (fabs( it1.Get()/it2.Get()-1 )>0.01)
        {
          itk::ImageFileWriter< mitk::PeakImage::ItkPeakImageType >::Pointer writer = itk::ImageFileWriter< mitk::PeakImage::ItkPeakImageType >::New();
          writer->SetInput(testImage);
          writer->SetFileName(mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + name);
          writer->Update();

          MITK_INFO << it1.Get() << " - " << it2.Get();

          CPPUNIT_ASSERT_MESSAGE("Peak images should be equal 1", false);
        }
      }
      else if (fabs(it1.Get())>0.0001)
      {
        itk::ImageFileWriter< mitk::PeakImage::ItkPeakImageType >::Pointer writer = itk::ImageFileWriter< mitk::PeakImage::ItkPeakImageType >::New();
        writer->SetInput(testImage);
        writer->SetFileName(mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + name);
        writer->Update();

        CPPUNIT_ASSERT_MESSAGE("Peak images should be equal 2", false);
      }

      ++it1;
      ++it2;
    }
  }

  /**
   * The L-BFGS-B fit needs the gradient of its cost: compare VnlCostFunction::gradf with central
   * differences of f on a small random problem. VOXEL_VARIANCE is left out: whether its voxel means
   * are differentiated is an open question (the gradient treats them as constant).
   */
  void RegularizationGradients()
  {
    auto randGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    randGen->SetSeed(3);
    const unsigned int rows = 20, cols = 6;
    vnl_sparse_matrix<double> A(rows, cols);
    vnl_vector<double> b(rows);
    for (unsigned int r = 0; r < rows; ++r)
    {
      b[r] = randGen->GetUniformVariate(0, 2);
      for (unsigned int c = 0; c < cols; ++c)
        if (randGen->GetUniformVariate(0, 1) < 0.5)
          A(r, c) = randGen->GetUniformVariate(0, 1);
    }
    vnl_vector<double> x(cols);
    for (unsigned int c = 0; c < cols; ++c)
      x[c] = randGen->GetUniformVariate(0.1, 1.0);

    for (auto regu : { VnlCostFunction::NONE, VnlCostFunction::MSM, VnlCostFunction::VARIANCE, VnlCostFunction::LASSO, VnlCostFunction::GROUP_LASSO, VnlCostFunction::GROUP_VARIANCE })
    {
      VnlCostFunction cost(cols);
      cost.SetProblem(A, b, 0.1, regu);
      cost.SetGroupSizes({ 2, 4 });
      vnl_vector<double> gradient(cols);
      cost.gradf(x, gradient);
      const double h = 1e-6;
      double maxError = 0, maxGradient = 0;
      for (unsigned int c = 0; c < cols; ++c)
      {
        vnl_vector<double> xp(x), xm(x);
        xp[c] += h;
        xm[c] -= h;
        const double numeric = (cost.f(xp) - cost.f(xm)) / (2 * h);
        maxError = std::max(maxError, std::fabs(numeric - gradient[c]));
        maxGradient = std::max(maxGradient, std::fabs(numeric));
      }
      MITK_INFO << "Regularization " << regu << ": gradient error " << maxError << " of " << maxGradient;
      CPPUNIT_ASSERT_MESSAGE("gradient of regularization " + std::to_string(regu), maxError <= 1e-6 * std::max(1.0, maxGradient));
    }
  }

  /** Fiber weights and output images of a fit do not depend on the number of OpenMP threads, bit for bit. */
  void SameAtThreadCounts()
  {
    // restores the OpenMP thread count for the later cases, also when an assertion throws
    struct ThreadCountGuard
    {
      int threads = omp_get_max_threads();
      ~ThreadCountGuard() { omp_set_num_threads(threads); }
    } threadCountGuard;

    auto fit = [this](int threads)
    {
      omp_set_num_threads(threads);
      setUp();   // fresh tractograms: the fit writes its weights into them
      fitter->SetRegularization(VnlCostFunction::VOXEL_VARIANCE);
      fitter->Update();
      return fitter;
    };
    FitterType::Pointer reference = fit(1);
    FitterType::Pointer test = fit(8);

    const auto referenceTracts = reference->GetTractograms();
    const auto testTracts = test->GetTractograms();
    CPPUNIT_ASSERT_EQUAL(referenceTracts.size(), testTracts.size());
    for (std::size_t t = 0; t < testTracts.size(); ++t)
      for (unsigned int i = 0; i < testTracts[t]->GetNumFibers(); ++i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("weight of fiber " + std::to_string(i) + " of bundle " + std::to_string(t),
                                     referenceTracts[t]->GetFiberWeight(i), testTracts[t]->GetFiberWeight(i));

    for (auto getter : {&FitterType::GetFittedImage, &FitterType::GetResidualImage, &FitterType::GetOverexplainedImage, &FitterType::GetUnderexplainedImage})
    {
      auto a = ((*reference).*getter)();
      auto b = ((*test).*getter)();
      CPPUNIT_ASSERT(a->GetLargestPossibleRegion() == b->GetLargestPossibleRegion());
      const auto n = a->GetLargestPossibleRegion().GetNumberOfPixels();
      CPPUNIT_ASSERT_MESSAGE("output image", std::memcmp(a->GetBufferPointer(), b->GetBufferPointer(), n*sizeof(*a->GetBufferPointer())) == 0);
    }
  }

  void Fit1()
  {
    omp_set_num_threads(1);
//    fitter->SetLambda(0.1);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::NONE);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/NONE_fitted_v2.fib");
    CompareFibs(test, ref, "NONE_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "NONE_fitted_image_v3.nrrd");
    CompareImages(fitter->GetResidualImage(), "NONE_residual_image_v3.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "NONE_overexplained_image_v4.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "NONE_underexplained_image_v4.nrrd");
  }

  void Fit2()
  {
    omp_set_num_threads(1);
//    fitter->SetLambda(0.1);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::MSM);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/MSM_fitted_v2.fib");
    CompareFibs(test, ref, "MSM_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "MSM_fitted_image_v3.nrrd");
    CompareImages(fitter->GetResidualImage(), "MSM_residual_image_v3.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "MSM_overexplained_image_v4.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "MSM_underexplained_image_v4.nrrd");
  }

  void Fit3()
  {
    omp_set_num_threads(1);
//    fitter->SetLambda(0.1);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::VARIANCE);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/MSE_fitted_v2.fib");
    CompareFibs(test, ref, "MSE_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "MSE_fitted_image_v4.nrrd");
    CompareImages(fitter->GetResidualImage(), "MSE_residual_image_v4.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "MSE_overexplained_image_v5.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "MSE_underexplained_image_v5.nrrd");
  }

  void Fit4()
  {
    omp_set_num_threads(1);
    fitter->SetLambda(100);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::VOXEL_VARIANCE);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/LocalMSE_fitted_v2.fib");
    CompareFibs(test, ref, "LocalMSE_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "LocalMSE_fitted_image_v3.nrrd");
    CompareImages(fitter->GetResidualImage(), "LocalMSE_residual_image_v3.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "LocalMSE_overexplained_image_v4.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "LocalMSE_underexplained_image_v4.nrrd");
  }

  void Fit5()
  {
    omp_set_num_threads(1);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::GROUP_VARIANCE);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/GroupMSE_fitted_v2.fib");
    CompareFibs(test, ref, "GroupMSE_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "GroupMSE_fitted_image_v3.nrrd");
    CompareImages(fitter->GetResidualImage(), "GroupMSE_residual_image_v3.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "GroupMSE_overexplained_image_v4.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "GroupMSE_underexplained_image_v4.nrrd");
  }

  void Fit6()
  {
    omp_set_num_threads(1);
    fitter->SetLambda(1000);
    fitter->SetFilterOutliers(false);
    fitter->SetRegularization(VnlCostFunction::GROUP_LASSO);
    fitter->Update();

    std::vector< mitk::FiberBundle::Pointer > output_tracts = fitter->GetTractograms();
    mitk::FiberBundle::Pointer test = mitk::FiberBundle::New();
    test = test->AddBundles(output_tracts);
    mitk::FiberBundle::Pointer ref = LoadFib("out/GroupLasso_fitted_v2.fib");
    CompareFibs(test, ref, "GroupLasso_fitted.fib");
    CompareImages(fitter->GetFittedImage(), "GroupLasso_fitted_image_v3.nrrd");
    CompareImages(fitter->GetResidualImage(), "GroupLasso_residual_image_v3.nrrd");
    CompareImages(fitter->GetOverexplainedImage(), "GroupLasso_overexplained_image_v4.nrrd");
    CompareImages(fitter->GetUnderexplainedImage(), "GroupLasso_underexplained_image_v4.nrrd");
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkFiberFit)
