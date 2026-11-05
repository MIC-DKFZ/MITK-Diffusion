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
#include <itkStreamlineTrackingFilter.h>
#include <Algorithms/TrackingHandlers/mitkTrackingDataHandler.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.h>
#include <itkDiffusionTensor3D.h>
#include <mitkTestingMacros.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkFiberBundle.h>
#include <mitkOmpThreadGuard.h>
#include <omp.h>
#include <itksys/SystemTools.hxx>
#include <mitkEqual.h>
#include <mitkStreamlineTractographyParameters.h>
#include <mitkStreamlineTractography.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkPeakImage.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionReconstruction.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <algorithm>

class mitkStreamlineTractographyTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkStreamlineTractographyTestSuite);
  MITK_TEST(Test_Peak1);
  MITK_TEST(Test_Peak2);
  MITK_TEST(Test_Tensor1);
  MITK_TEST(Test_Tensor2);
  MITK_TEST(Test_Tensor3);
  MITK_TEST(Test_Odf1);
  MITK_TEST(Test_Odf2);
  MITK_TEST(Test_Odf3);
  MITK_TEST(Test_Odf4);
  MITK_TEST(Test_Odf5);
  MITK_TEST(Test_Odf6);
  MITK_TEST(Test_Tensor1_FourThreads);
  MITK_TEST(Test_Odf5_FourThreads);
  MITK_TEST(TensorInput_FollowsAFitOfTheTestDwi);
  MITK_TEST(SharedSetup_Peaks_SameAsTest_Peak1);
  MITK_TEST(SharedSetup_Tensors_SameAsTest_Tensor1);
  MITK_TEST(SharedSetup_Odf_SameAsTest_Odf1);
  MITK_TEST(SharedSetup_TensorsProbabilistic_TracksTheirOdfs);
  MITK_TEST(SharedSetup_Sh_TracksItsOdfs);
  MITK_TEST(SharedSetup_Prior_FlipsAndNoCutoff);
  MITK_TEST(SharedSetup_UnsupportedInputs_Throw);
  MITK_TEST(Filter_CountsSeedsAndFibers);
  CPPUNIT_TEST_SUITE_END();

  typedef itk::VectorImage< short, 3>   ItkDwiType;

private:

public:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/

  typedef itk::Image<float, 3> ItkFloatImgType;

  mitk::TrackingHandlerOdf::ItkOdfImageType::Pointer itk_odf_image;
  mitk::TrackingHandlerTensor::ItkTensorImageType::ConstPointer itk_tensor_image;
  mitk::TrackingHandlerPeaks::PeakImgType::Pointer itk_peak_image;
  ItkFloatImgType::Pointer itk_seed_image;
  ItkFloatImgType::Pointer itk_mask_image;
  ItkFloatImgType::Pointer itk_gfa_image;

  float gfa_threshold;
  float odf_threshold;
  float peak_threshold;
  std::shared_ptr<mitk::StreamlineTractographyParameters> params;

  itk::StreamlineTrackingFilter::Pointer tracker;
  std::unique_ptr<mitk::OmpThreadGuard> m_Threads;  ///< 1 thread per case, restored in tearDown

  // The images as loaded, for the shared setup of view and app (mitkStreamlineTractography.h).
  mitk::Image::Pointer m_OdfImage;
  mitk::Image::Pointer m_TensorImage;
  mitk::Image::Pointer m_SeedImage;
  mitk::Image::Pointer m_MaskImage;
  mitk::Image::Pointer m_GfaImage;

  void setUp() override
  {
    m_Threads = std::make_unique<mitk::OmpThreadGuard>(1);

    gfa_threshold = 0.2f;
    odf_threshold = 0.1f;
    peak_threshold = 0.1f;

    mitk::Image::Pointer odf_image = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/qball_image.qbi"));
    mitk::Image::Pointer tensor_image = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/tensor_image.dti"));
    mitk::Image::Pointer peak_image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/qball_peak_image.nii.gz"));
    mitk::Image::Pointer seed_image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/seed_image.nii.gz"));
    mitk::Image::Pointer mask_image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/mask_image.nii.gz"));
    mitk::Image::Pointer gfa_image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/gfa_image.nii.gz"));

    m_OdfImage = odf_image;
    m_TensorImage = tensor_image;
    m_SeedImage = seed_image;
    m_MaskImage = mask_image;
    m_GfaImage = gfa_image;

    params = std::make_shared<mitk::StreamlineTractographyParameters>();
    params->m_FixRandomSeed = true;
    params->m_SecondOrder = false;
    params->m_InterpolateRoiImages = false;
    params->SetLoopCheckDeg(-1);  // todo: test loop check

    {
      typedef mitk::ImageToItk< mitk::TrackingHandlerPeaks::PeakImgType > CasterType;
      CasterType::Pointer caster = CasterType::New();
      caster->SetInput(peak_image);
      caster->Update();
      itk_peak_image = caster->GetOutput();
    }

    {
      typedef mitk::ImageToItk< mitk::TrackingHandlerTensor::ItkTensorImageType > CasterType;
      CasterType::Pointer caster = CasterType::New();
      caster->SetInput(tensor_image);
      caster->Update();
      itk_tensor_image = caster->GetOutput();
    }

    {
      typedef mitk::ImageToItk< mitk::TrackingHandlerOdf::ItkOdfImageType > CasterType;
      CasterType::Pointer caster = CasterType::New();
      caster->SetInput(odf_image);
      caster->Update();
      itk_odf_image = caster->GetOutput();
    }

    itk_gfa_image = ItkFloatImgType::New();
    mitk::CastToItkImage(gfa_image, itk_gfa_image);

    itk_seed_image = ItkFloatImgType::New();
    mitk::CastToItkImage(seed_image, itk_seed_image);

    itk_mask_image = ItkFloatImgType::New();
    mitk::CastToItkImage(mask_image, itk_mask_image);
  }

  mitk::FiberBundle::Pointer LoadReferenceFib(std::string filename)
  {
    mitk::FiberBundle::Pointer fib = nullptr;

    if (itksys::SystemTools::FileExists(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/ReferenceFibs/" + filename)))
    {
      mitk::BaseData::Pointer baseData = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/ReferenceFibs/" + filename)).at(0);
      fib = dynamic_cast<mitk::FiberBundle*>(baseData.GetPointer());
    }
    return fib;
  }

  mitk::Image::Pointer LoadReferenceImage(std::string filename)
  {
    mitk::Image::Pointer img = nullptr;

    if (itksys::SystemTools::FileExists(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/ReferenceFibs/" + filename)))
    {
      img = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/ReferenceFibs/" + filename));
    }
    return img;
  }

  void SetupTracker(mitk::TrackingDataHandler* handler)
  {

    tracker = itk::StreamlineTrackingFilter::New();

//    tracker->SetInterpolateMasks(false);
//    tracker->SetNumberOfSamples(0);
//    tracker->SetAngularThreshold(-1);
    tracker->SetMaskImage(itk_mask_image);
    tracker->SetSeedImage(itk_seed_image);
    tracker->SetStoppingRegions(nullptr);
//    tracker->SetSeedsPerVoxel(1);
//    tracker->SetStepSize(0.5);
//    tracker->SetSamplingDistance(0.25);
//    tracker->SetUseStopVotes(true);
//    tracker->SetOnlyForwardSamples(true);
//    tracker->SetMinTractLength(20);
//    tracker->SetMaxNumTracts(-1);
    tracker->SetTrackingHandler(handler);
//    tracker->SetUseOutputProbabilityMap(false);
    tracker->SetParameters(params);
  }

  void tearDown() override
  {
    m_Threads.reset();
  }

  /** The peak image as the Workbench and the app read it: a mitk::PeakImage, not a plain 4D image. */
  static mitk::Image::Pointer LoadPeakImage(const std::string& file)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath(file), &functor).GetPointer();
  }

  /**
   * Releases the ITK views of setUp on the tensor and ODF images: a view holds a read access, and
   * MITK refuses a second access of the same thread to the same image, which the shared setup takes.
   */
  void ReleaseItkViews()
  {
    itk_tensor_image = nullptr;
    itk_odf_image = nullptr;
  }

  /** Tracks with the shared setup: handler from the model images, seed and mask ROIs as mitk images. */
  mitk::FiberBundle::Pointer TrackWithSharedSetup(const std::vector<const mitk::Image*>& modelImages,
                                                  mitk::StreamlineTractographyParameters::MODE mode,
                                                  const mitk::Image* faImage = nullptr)
  {
    params->m_Mode = mode;
    ReleaseItkViews();
    auto handler = mitk::CreateTrackingHandler(modelImages, mode, faImage);

    mitk::StreamlineTractographyRois rois;
    rois.seed = m_SeedImage;
    rois.mask = m_MaskImage;

    tracker = itk::StreamlineTrackingFilter::New();
    mitk::SetTrackingRois(*tracker, rois);
    tracker->SetTrackingHandler(handler.get());
    tracker->SetParameters(params);
    tracker->Update();

    return mitk::FiberBundle::New(tracker->GetFiberPolyData());
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

  void CheckImageResult(std::string ref_file, mitk::Image::Pointer test_img)
  {
    mitk::Image::Pointer ref = LoadReferenceImage(ref_file);
    if (ref.IsNull())
    {
      mitk::IOUtil::Save(test_img, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
      CPPUNIT_FAIL("Reference file not found. Saving test file to " + mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + ref_file);
    }
    else
    {
      MITK_ASSERT_EQUAL(test_img, ref, "Images should be equal");
    }
  }

  void Test_Peak1()
  {
    mitk::TrackingHandlerPeaks* handler = new mitk::TrackingHandlerPeaks();
    handler->SetPeakImage(itk_peak_image);
    params->m_Cutoff = peak_threshold;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Peak1.fib", outFib);

    delete handler;
  }

  void Test_Peak2()
  {
    mitk::TrackingHandlerPeaks* handler = new mitk::TrackingHandlerPeaks();
    handler->SetPeakImage(itk_peak_image);
    params->m_Cutoff = peak_threshold;
    params->m_InterpolateTractographyData = false;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Peak2.fib", outFib);

    delete handler;
  }

  void Test_Tensor1()
  {
    mitk::TrackingHandlerTensor* handler = new mitk::TrackingHandlerTensor();
    handler->SetTensorImage(itk_tensor_image);
    params->m_Cutoff = gfa_threshold;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Tensor1.fib", outFib);

    delete handler;
  }

  void Test_Tensor2()
  {
    mitk::TrackingHandlerTensor* handler = new mitk::TrackingHandlerTensor();
    handler->SetTensorImage(itk_tensor_image);
    params->m_Cutoff = gfa_threshold;
    params->m_InterpolateTractographyData = false;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Tensor2.fib", outFib);

    delete handler;
  }

  void Test_Tensor3()
  {
    mitk::TrackingHandlerTensor* handler = new mitk::TrackingHandlerTensor();
    handler->SetTensorImage(itk_tensor_image);
    params->m_Cutoff = gfa_threshold;
    params->m_InterpolateTractographyData = false;
    params->m_F = 0;
    params->m_G = 1;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Tensor3.fib", outFib);

    delete handler;
  }

  void Test_Odf1()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 1;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Odf1_v2.fib", outFib);  // regenerated when deterministic ODF peaks were refined between the sampling directions

    delete handler;
  }

  void Test_Odf2()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 8;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Odf2_v2.fib", outFib);  // regenerated when deterministic ODF peaks were refined between the sampling directions

    delete handler;
  }


  void Test_Odf3()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 1;
    params->m_InterpolateTractographyData = false;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Odf3_v2.fib", outFib);  // regenerated when deterministic ODF peaks were refined between the sampling directions

    delete handler;
  }

  void Test_Odf4()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 1;
    params->m_SeedsPerVoxel = 3;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    CheckFibResult("Test_Odf4_v2.fib", outFib);  // regenerated when deterministic ODF peaks were refined between the sampling directions

    delete handler;
  }

  void Test_Odf5()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SeedsPerVoxel = 10;
    params->m_SharpenOdfs = 8;
    params->m_Mode = mitk::TrackingDataHandler::MODE::PROBABILISTIC;

    SetupTracker(handler);
    tracker->Update();

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);

    // Regenerated when every seed point got its own random streams (results independent of the thread count).
    CheckFibResult("Test_Odf5_v2.fib", outFib);

    delete handler;
  }

  void Test_Odf6()
  {
    mitk::TrackingHandlerOdf* handler = new mitk::TrackingHandlerOdf();
    handler->SetOdfImage(itk_odf_image);

    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SeedsPerVoxel = 10;
    params->m_SharpenOdfs = 8;
    params->m_Mode = mitk::TrackingDataHandler::MODE::PROBABILISTIC;
    params->m_OutputProbMap = true;

    SetupTracker(handler);
    tracker->Update();

    itk::StreamlineTrackingFilter::ItkDoubleImgType::Pointer outImg = tracker->GetOutputProbabilityMap();
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg.GetPointer());
    img->SetVolume(outImg->GetBufferPointer());

    CheckImageResult("Test_Odf6_v2.nrrd", img);

    delete handler;
  }


  /** Same as Test_Tensor1 at four threads: the fiber order must not depend on the thread scheduling. */
  void Test_Tensor1_FourThreads()
  {
    m_Threads->SetNumberOfThreads(4);
    Test_Tensor1();
  }

  /** Same as Test_Odf5 at four threads: probabilistic tracking must give the single-threaded result. */
  void Test_Odf5_FourThreads()
  {
    m_Threads->SetNumberOfThreads(4);
    Test_Odf5();
  }

  /**
   * tensor_image.dti is the linear least-squares fit of diffusion_image_v2 (the test DWI with FSL
   * bvecs), so its directions follow the phantom's bundles. diffusion_image.bvecs holds the same
   * gradients with x not negated, as earlier MITK versions wrote NIfTI bvecs: read as FSL bvecs,
   * they give directions mirrored in x, 47 degrees (median) off the bundles.
   */
  void TensorInput_FollowsAFitOfTheTestDwi()
  {
    auto dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/diffusion_image_v2.nii.gz"));
    auto fitted = mitk::ReconstructTensors(dwi, mitk::TensorFit::LinearLeastSquares, 0, false);
    auto itkFitted = mitk::ImageToItkImage<mitk::TrackingHandlerTensor::TensorType, 3>(fitted.GetPointer());

    std::vector<double> angles;
    itk::ImageRegionConstIteratorWithIndex<mitk::TrackingHandlerTensor::ItkTensorImageType> it(itk_tensor_image, itk_tensor_image->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      auto input = it.Get();
      if (input.GetFractionalAnisotropy() < 0.4)
        continue;
      mitk::TrackingHandlerTensor::TensorType::EigenValuesArrayType values;
      mitk::TrackingHandlerTensor::TensorType::EigenVectorsMatrixType inputVectors, fittedVectors;
      input.ComputeEigenAnalysis(values, inputVectors);
      itkFitted->GetPixel(it.GetIndex()).ComputeEigenAnalysis(values, fittedVectors);
      double dot = 0;
      for (int c = 0; c < 3; ++c)
        dot += inputVectors(2, c) * fittedVectors(2, c);
      angles.push_back(std::acos(std::min(1.0, std::abs(dot))) * 180.0 / itk::Math::pi);
    }
    CPPUNIT_ASSERT(angles.size() > 500);
    std::nth_element(angles.begin(), angles.begin() + angles.size() / 2, angles.end());
    CPPUNIT_ASSERT_MESSAGE("median angle " + std::to_string(angles[angles.size() / 2]), angles[angles.size() / 2] < 1.0);
  }

  void SharedSetup_Peaks_SameAsTest_Peak1()
  {
    params->m_Cutoff = peak_threshold;
    const auto peaks = LoadPeakImage("StreamlineTractography/qball_peak_image.nii.gz");
    CheckFibResult("Test_Peak1.fib", TrackWithSharedSetup({peaks.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC));
  }

  void SharedSetup_Tensors_SameAsTest_Tensor1()
  {
    params->m_Cutoff = gfa_threshold;
    CheckFibResult("Test_Tensor1.fib", TrackWithSharedSetup({m_TensorImage.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC));
  }

  void SharedSetup_Odf_SameAsTest_Odf1()
  {
    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 1;
    CheckFibResult("Test_Odf1_v2.fib", TrackWithSharedSetup({m_OdfImage.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC));
  }

  void SharedSetup_TensorsProbabilistic_TracksTheirOdfs()
  {
    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SeedsPerVoxel = 2;
    params->m_Mode = mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC;

    // The direct setup of the view and the app before the shared one: the tensors' ODFs with the GFA image.
    ReleaseItkViews();
    mitk::TrackingHandlerOdf direct;
    direct.SetOdfImage(mitk::convert::GetItkOdfFromTensorImage(m_TensorImage));
    direct.SetIsOdfFromTensor(true);
    direct.SetGfaImage(itk_gfa_image);
    SetupTracker(&direct);
    tracker->Update();
    mitk::FiberBundle::Pointer reference = mitk::FiberBundle::New(tracker->GetFiberPolyData());
    CPPUNIT_ASSERT(reference->GetNumFibers() > 0);

    auto handler = mitk::CreateTrackingHandler({m_TensorImage.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC);
    auto* odfHandler = dynamic_cast<mitk::TrackingHandlerOdf*>(handler.get());
    CPPUNIT_ASSERT(odfHandler != nullptr && odfHandler->GetIsOdfFromTensor());

    const auto shared = TrackWithSharedSetup({m_TensorImage.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC, m_GfaImage);
    CPPUNIT_ASSERT(reference->Equals(shared));
  }

  void SharedSetup_Sh_TracksItsOdfs()
  {
    params->m_Cutoff = gfa_threshold;
    params->m_OdfCutoff = 0;
    params->m_SharpenOdfs = 1;

    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    mitk::Image::Pointer sh = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/qball_image_sh.nii.gz"), &functor).GetPointer();
    CPPUNIT_ASSERT(sh.IsNotNull());

    mitk::TrackingHandlerOdf direct;
    direct.SetOdfImage(mitk::convert::GetItkOdfFromShImage(sh));
    direct.SetGfaImage(itk_gfa_image);
    SetupTracker(&direct);
    tracker->Update();
    mitk::FiberBundle::Pointer reference = mitk::FiberBundle::New(tracker->GetFiberPolyData());
    CPPUNIT_ASSERT(reference->GetNumFibers() > 0);

    const auto shared = TrackWithSharedSetup({sh.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC, m_GfaImage);
    CPPUNIT_ASSERT(reference->Equals(shared));
  }

  void SharedSetup_Prior_FlipsAndNoCutoff()
  {
    const auto peaks = LoadPeakImage("StreamlineTractography/qball_peak_image.nii.gz");
    auto prior = mitk::CreateTrackingPriorHandler(peaks);
    CPPUNIT_ASSERT(dynamic_cast<mitk::TrackingHandlerPeaks*>(prior.get()) != nullptr);

    params->m_Cutoff = 0.3f;
    params->m_PriorFlipX = true;
    params->m_FlipY = true;
    params->m_Weight = 0.7f;
    mitk::SetTrackingPriorParameters(*prior, *params);

    const auto priorParameters = prior->GetParameters();
    CPPUNIT_ASSERT(priorParameters->m_FlipX && !priorParameters->m_FlipY && !priorParameters->m_FlipZ);
    CPPUNIT_ASSERT_EQUAL(0.0f, priorParameters->m_Cutoff);
    CPPUNIT_ASSERT_EQUAL(0.7f, priorParameters->m_Weight);
    // The run's own parameters stay as they are.
    CPPUNIT_ASSERT_EQUAL(0.3f, params->m_Cutoff);
    CPPUNIT_ASSERT(params->m_FlipY && !params->m_FlipX);

    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingPriorHandler(m_GfaImage), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingPriorHandler(nullptr), mitk::Exception);
  }

  void SharedSetup_UnsupportedInputs_Throw()
  {
    using Mode = mitk::StreamlineTractographyParameters::MODE;
    ReleaseItkViews();
    const auto peaks = LoadPeakImage("StreamlineTractography/qball_peak_image.nii.gz");
    const mitk::Image* tensors = m_TensorImage;
    // A second tensor image: the same image twice is refused by MITK's image access.
    const mitk::Image::Pointer moreTensors = m_TensorImage->Clone();

    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({}, Mode::DETERMINISTIC), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({nullptr}, Mode::DETERMINISTIC), mitk::Exception);
    // a scalar image is no model
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({m_GfaImage.GetPointer()}, Mode::DETERMINISTIC), mitk::Exception);
    // several images only as deterministic multi-tensor model
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({tensors, moreTensors.GetPointer()}, Mode::PROBABILISTIC), mitk::Exception);
    const mitk::Image::Pointer morePeaks = peaks->Clone();
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({peaks.GetPointer(), morePeaks.GetPointer()}, Mode::DETERMINISTIC), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::CreateTrackingHandler({tensors, m_OdfImage.GetPointer()}, Mode::DETERMINISTIC), mitk::Exception);

    auto multiTensor = mitk::CreateTrackingHandler({tensors, moreTensors.GetPointer()}, Mode::DETERMINISTIC);
    auto* tensorHandler = dynamic_cast<mitk::TrackingHandlerTensor*>(multiTensor.get());
    CPPUNIT_ASSERT(tensorHandler != nullptr);
    CPPUNIT_ASSERT_EQUAL(2, tensorHandler->GetNumTensorImages());
  }

  void Filter_CountsSeedsAndFibers()
  {
    params->m_Cutoff = peak_threshold;
    const auto peaks = LoadPeakImage("StreamlineTractography/qball_peak_image.nii.gz");
    const auto fibers = TrackWithSharedSetup({peaks.GetPointer()}, mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC);

    CPPUNIT_ASSERT(tracker->GetNumSeeds() > 0);
    CPPUNIT_ASSERT_EQUAL(tracker->GetNumSeeds(), tracker->GetProcessedSeeds());
    CPPUNIT_ASSERT_EQUAL(fibers->GetNumFibers(), tracker->GetAcceptedFibers());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkStreamlineTractography)
