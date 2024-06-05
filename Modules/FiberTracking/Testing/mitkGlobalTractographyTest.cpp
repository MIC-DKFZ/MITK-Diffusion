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
#include <mitkTestFixture.h>
#include <mitkFiberBundle.h>
#include <omp.h>
#include <itksys/SystemTools.hxx>
#include <mitkEqual.h>
#include <itkGibbsTrackingFilter.h>

class mitkGlobalTractographyTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkGlobalTractographyTestSuite);
  MITK_TEST(Test_Odf);
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

    mitk::Image::Pointer odf_image = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/qBallImage.qbi"));
    mitk::Image::Pointer mask_image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/diffusionImageMask.nrrd"));

    itk_mask_image = MaskImgType::New();
    mitk::CastToItkImage(mask_image, itk_mask_image);


    OdfVectorImgType::Pointer itk_odf = OdfVectorImgType::New();
    mitk::CastToItkImage(odf_image, itk_odf_image);
  }

  mitk::FiberBundle::Pointer LoadReferenceFib(std::string filename)
  {
    mitk::FiberBundle::Pointer fib = nullptr;

    if (itksys::SystemTools::FileExists(GetTestDataFilePath("DiffusionImaging/" + filename)))
    {
      mitk::BaseData::Pointer baseData = mitk::IOUtil::Load(GetTestDataFilePath("DiffusionImaging/" + filename)).at(0);
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
      mitk::IOUtil::Save(test_fib, mitk::IOUtil::GetTempPath()+ref_file);
      CPPUNIT_FAIL("Reference file not found. Saving test file to " + mitk::IOUtil::GetTempPath() + ref_file);
    }
    else
    {
      bool is_equal = ref->Equals(test_fib);
      if (!is_equal)
      {
        mitk::IOUtil::Save(test_fib, mitk::IOUtil::GetTempPath()+ref_file);
        CPPUNIT_FAIL("Tractograms are not equal! Saving test file to " + mitk::IOUtil::GetTempPath() + ref_file);
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
    tracker->SetLoadParameterFile(GetTestDataFilePath("DiffusionImaging/gibbsTrackingParameters.gtp"));
    tracker->Update();

    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(tracker->GetFiberBundle());

    CheckFibResult("gibbsTractogram.fib", outFib);
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkGlobalTractography)
