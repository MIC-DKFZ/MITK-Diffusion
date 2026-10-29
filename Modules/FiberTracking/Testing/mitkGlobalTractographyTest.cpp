/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkTensorImage.h>
#include <mitkOdfImage.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkTestingMacros.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkFiberBundle.h>
#include <omp.h>
#include <itksys/SystemTools.hxx>
#include <mitkEqual.h>
#include <itkGibbsTrackingFilter.h>
#include <atomic>
#include <chrono>
#include <thread>

class mitkGlobalTractographyTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkGlobalTractographyTestSuite);
  MITK_TEST(Test_Odf);
  MITK_TEST(Test_Odf_FibersRequestedDuringRun);
  CPPUNIT_TEST_SUITE_END();

  typedef itk::VectorImage< short, 3>   ItkDwiType;

  private:

  public:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/

  typedef itk::Vector<float, ODF_SAMPLING_SIZE> OdfVectorType;
  typedef itk::Image<OdfVectorType,3> OdfVectorImgType;
  typedef itk::Image<float,3> MaskImgType;
  typedef itk::GibbsTrackingFilter<OdfVectorImgType> GibbsTrackingFilterType;

  GibbsTrackingFilterType::Pointer tracker;

  OdfVectorImgType::Pointer itk_odf_image;
  MaskImgType::Pointer itk_mask_image;

  void setUp() override
  {
    omp_set_num_threads(1);

    mitk::Image::Pointer odf_image = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi"));
    mitk::Image::Pointer mask_image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("diffusionImageMask.nrrd"));

    itk_mask_image = MaskImgType::New();
    mitk::CastToItkImage(mask_image, itk_mask_image);


    OdfVectorImgType::Pointer itk_odf = OdfVectorImgType::New();
    mitk::CastToItkImage(odf_image, itk_odf_image);
  }

  mitk::FiberBundle::Pointer LoadReferenceFib(std::string filename)
  {
    mitk::FiberBundle::Pointer fib = nullptr;

    if (itksys::SystemTools::FileExists(mitk::GetDiffusionTestDataFilePath(filename)))
    {
      mitk::BaseData::Pointer baseData = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(filename)).at(0);
      fib = dynamic_cast<mitk::FiberBundle*>(baseData.GetPointer());
    }
    return fib;
  }

  void tearDown() override
  {

  }

  void CheckFibResult(std::string ref_file, mitk::FiberBundle::Pointer test_fib)
  {
    mitk::FiberBundle::Pointer ref = LoadReferenceFib(ref_file);
    if (ref.IsNull())
    {
      mitk::IOUtil::Save(test_fib, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
      CPPUNIT_FAIL("Reference file not found. Saving test file to " + mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
    }
    else
    {
      bool is_equal = ref->Equals(test_fib);
      if (!is_equal)
      {
        mitk::IOUtil::Save(test_fib, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
        CPPUNIT_FAIL("Tractograms are not equal! Saving test file to " + mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
      }
    }
  }


  void Test_Odf()
  {
    tracker = GibbsTrackingFilterType::New();
    tracker->SetOdfImage(itk_odf_image.GetPointer());
    tracker->SetMaskImage(itk_mask_image);
    tracker->SetDuplicateImage(false);
    tracker->SetRandomSeed(1);
    tracker->SetLoadParameterFile(mitk::GetDiffusionTestDataFilePath("gibbsTrackingParameters.gtp"));
    tracker->Update();

    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(tracker->GetFiberBundle());

    // Regenerated when FiberBuilder kept the first particle of every fiber (v2) and when the optimal
    // shift with one neighbour faced that neighbour (v3).
    CheckFibResult("gibbsTractogram_v3.fib", outFib);
  }

  /** The Workbench view asks for intermediate fibers from the GUI thread while Update() runs in a
   *  worker thread. Every request must return, also the ones racing with the end of the run, and
   *  the intermediate fiber builds must not change the result. */
  void Test_Odf_FibersRequestedDuringRun()
  {
    tracker = GibbsTrackingFilterType::New();
    tracker->SetOdfImage(itk_odf_image.GetPointer());
    tracker->SetMaskImage(itk_mask_image);
    tracker->SetDuplicateImage(false);
    tracker->SetRandomSeed(1);
    tracker->SetLoadParameterFile(mitk::GetDiffusionTestDataFilePath("gibbsTrackingParameters.gtp"));

    std::atomic<bool> done(false);
    std::string error;
    std::thread worker([&]() {
      try
      {
        tracker->Update();
      }
      catch (const std::exception &e)
      {
        error = e.what();
      }
      done = true;
    });
    int snapshots = 0;
    while (!done)
    {
      if (tracker->GetFiberBundle() != nullptr)
        ++snapshots;
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    worker.join();
    CPPUNIT_ASSERT_MESSAGE(error, error.empty());
    MITK_INFO << snapshots << " intermediate tractograms";

    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(tracker->GetFiberBundle());
    CheckFibResult("gibbsTractogram_v3.fib", outFib);
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkGlobalTractography)
