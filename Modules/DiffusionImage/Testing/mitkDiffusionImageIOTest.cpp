/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkFileReaderSelector.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include <itksys/SystemTools.hxx>

#include <format>

/** Reading and writing of DWIs as NRRD and NIfTI, on MITK-Data input and on synthetic files. */
class mitkDiffusionImageIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageIOTestSuite);
  MITK_TEST(NrrdReader_MitkFrameKey_ReadsFrame);
  MITK_TEST(NrrdReader_StandardFrameField_ReadsFrame);
  MITK_TEST(NrrdReader_NoFrame_ReadsIdentity);
  MITK_TEST(NrrdReader_BothFrameKinds_ReadsIdentity);
  MITK_TEST(NrrdReader_GradientCountMismatch_Throws);
  MITK_TEST(NrrdReader_MissingBValue_Throws);
  MITK_TEST(NrrdReader_DetachedHeader_Reads);
  MITK_TEST(NrrdWriter_RoundTrip_KeepsDataAndGradients);
  MITK_TEST(NrrdWriter_OnlyModifiedGradients_WritesThem);
  MITK_TEST(NiftiWriter_RoundTrip_KeepsDataAndGradients);
  MITK_TEST(PlainNrrdWriter_Dwi_WritesNoDiffusionKeys);
  MITK_TEST(ReaderSelection_PlainImages_NoDwiReader);
  MITK_TEST(ReaderSelection_DiffusionFiles_DwiReaderIsDefault);
  MITK_TEST(WriterSelection_Dwi_DwiWriterIsDefault);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;

  mitk::Image::Pointer m_Dwi;
  std::string m_TempDir;

  static mitk::PreferenceListReaderOptionsFunctor DwiReaderFunctor()
  {
    return mitk::PreferenceListReaderOptionsFunctor(mitk::PreferenceListReaderOptionsFunctor::ListType({"Diffusion Weighted Images"}),
                                                    mitk::PreferenceListReaderOptionsFunctor::ListType());
  }

  static double MaxDifference(const PropHelper::GradientDirectionsContainerType* a, const PropHelper::GradientDirectionsContainerType* b)
  {
    CPPUNIT_ASSERT_EQUAL(a->Size(), b->Size());
    double result = 0;
    for (unsigned int i = 0; i < a->Size(); ++i)
      result = std::max(result, (a->ElementAt(i) - b->ElementAt(i)).two_norm());
    return result;
  }

  static const std::string& DwiReaderDescription()
  {
    static const std::string description = "Diffusion Weighted Images";
    return description;
  }

  static bool OffersDwiReader(const std::string& path)
  {
    mitk::FileReaderSelector selector(path);
    for (const auto& item : selector.Get())
    {
      if (item.GetDescription() == DwiReaderDescription())
        return true;
    }
    return false;
  }

  static std::string DefaultReader(const std::string& path)
  {
    return mitk::FileReaderSelector(path).GetDefault().GetDescription();
  }

  static PropHelper::MeasurementFrameType Identity()
  {
    PropHelper::MeasurementFrameType identity;
    identity.set_identity();
    return identity;
  }

  /**
   * Writes the pixels of m_Dwi with only the given key/value pairs through plain ITK, the way other
   * tools write DWI NRRDs.
   */
  std::string WriteNrrd(const std::string& name, const itk::MetaDataDictionary& dictionary)
  {
    auto image = PropHelper::ImageType::New();
    image->Graft(PropHelper::GetItkVectorImage(m_Dwi));
    image->SetMetaDataDictionary(dictionary);
    auto writer = itk::ImageFileWriter<PropHelper::ImageType>::New();
    writer->SetInput(image);
    writer->SetImageIO(itk::NrrdImageIO::New());
    writer->UseInputMetaDataDictionaryOn();
    const auto path = m_TempDir + "/" + name;
    writer->SetFileName(path);
    writer->Update();
    return path;
  }

  /** modality, b-value and the first count gradients of m_Dwi. */
  itk::MetaDataDictionary DiffusionKeys(unsigned int count)
  {
    itk::MetaDataDictionary dictionary;
    itk::EncapsulateMetaData<std::string>(dictionary, "modality", "DWMRI");
    itk::EncapsulateMetaData<std::string>(dictionary, "DWMRI_b-value", "3000");
    auto gradients = PropHelper::GetOriginalGradientContainer(m_Dwi);
    for (unsigned int i = 0; i < count; ++i)
    {
      const auto& g = gradients->ElementAt(i);
      itk::EncapsulateMetaData<std::string>(dictionary, std::format("DWMRI_gradient_{:04}", i), std::format("{} {} {}", g[0], g[1], g[2]));
    }
    return dictionary;
  }

  void AssertSameDwi(mitk::Image* expected, mitk::Image* actual, double gradientTolerance)
  {
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(actual));
    CPPUNIT_ASSERT(mitk::Equal(*expected, *actual, mitk::eps, true));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(PropHelper::GetReferenceBValue(expected), PropHelper::GetReferenceBValue(actual), 1e-6);
    CPPUNIT_ASSERT(MaxDifference(PropHelper::GetOriginalGradientContainer(expected), PropHelper::GetOriginalGradientContainer(actual)) <= gradientTolerance);
    CPPUNIT_ASSERT(MaxDifference(PropHelper::GetGradientContainer(expected), PropHelper::GetGradientContainer(actual)) <= gradientTolerance);
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(expected) == PropHelper::GetMeasurementFrame(actual));
    CPPUNIT_ASSERT(PropHelper::GetBValueMap(expected) == PropHelper::GetBValueMap(actual));
  }

public:
  void setUp() override
  {
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("DiffusionImageIOTest_XXXXXX");
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void NrrdReader_MitkFrameKey_ReadsFrame()
  {
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(m_Dwi));
    CPPUNIT_ASSERT_EQUAL(size_t(31), m_Dwi->GetPixelType().GetNumberOfComponents());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3000.0, PropHelper::GetReferenceBValue(m_Dwi), 1e-6);
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(m_Dwi) == Identity());

    auto dictionary = DiffusionKeys(31);
    itk::EncapsulateMetaData<std::string>(dictionary, "measurement frame", " (0,1,0) (-1,0,0) (0,0,1)");
    auto image = mitk::IOUtil::Load<mitk::Image>(WriteNrrd("mitk_frame.dwi", dictionary));
    PropHelper::MeasurementFrameType expected;
    expected.set_row(0, std::vector<double>{0, 1, 0}.data());
    expected.set_row(1, std::vector<double>{-1, 0, 0}.data());
    expected.set_row(2, std::vector<double>{0, 0, 1}.data());
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(image) == expected);
  }

  void NrrdReader_StandardFrameField_ReadsFrame()
  {
    auto dictionary = DiffusionKeys(31);
    const std::vector<std::vector<double>> frame = {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}};
    itk::EncapsulateMetaData<std::vector<std::vector<double>>>(dictionary, "NRRD_measurement frame", frame);
    auto image = mitk::IOUtil::Load<mitk::Image>(WriteNrrd("standard_frame.dwi", dictionary));

    const auto mf = PropHelper::GetMeasurementFrame(image);
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(frame[r][c], mf(r, c), 1e-12);
    CPPUNIT_ASSERT(MaxDifference(PropHelper::GetOriginalGradientContainer(m_Dwi), PropHelper::GetOriginalGradientContainer(image)) < 1e-12);
  }

  void NrrdReader_NoFrame_ReadsIdentity()
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(WriteNrrd("no_frame.dwi", DiffusionKeys(31)));
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(image) == Identity());
    CPPUNIT_ASSERT(MaxDifference(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(image)) < 1e-12);
  }

  void NrrdReader_BothFrameKinds_ReadsIdentity()
  {
    // Has the standard field and MITK's key/value pair.
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(image) == Identity());
  }

  void NrrdReader_GradientCountMismatch_Throws()
  {
    const auto path = WriteNrrd("missing_gradient.dwi", DiffusionKeys(30));
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load<mitk::Image>(path), mitk::Exception);
  }

  void NrrdReader_MissingBValue_Throws()
  {
    auto dictionary = DiffusionKeys(31);
    dictionary.Erase("DWMRI_b-value");
    const auto path = WriteNrrd("missing_bvalue.dwi", dictionary);
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load<mitk::Image>(path), mitk::Exception);
  }

  void NrrdReader_DetachedHeader_Reads()
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/diffusionImage.hdwi"));
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT_EQUAL(size_t(130), image->GetPixelType().GetNumberOfComponents());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(650.0, PropHelper::GetReferenceBValue(image), 1e-6);
  }

  void NrrdWriter_RoundTrip_KeepsDataAndGradients()
  {
    const auto path = m_TempDir + "/roundtrip.dwi";
    mitk::IOUtil::Save(m_Dwi, "DWI_NRRD", path);
    auto image = mitk::IOUtil::Load<mitk::Image>(path);
    AssertSameDwi(m_Dwi, image, 0);
  }

  void NrrdWriter_OnlyModifiedGradients_WritesThem()
  {
    // A filter output with modified gradients but no InitializeImage yet.
    auto source = m_Dwi->Clone();
    source->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(source));

    const auto path = m_TempDir + "/modified_only.dwi";
    mitk::IOUtil::Save(source, "DWI_NRRD", path);
    auto image = mitk::IOUtil::Load<mitk::Image>(path);
    CPPUNIT_ASSERT(MaxDifference(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(image)) == 0);
  }

  void NiftiWriter_RoundTrip_KeepsDataAndGradients()
  {
    auto functor = DwiReaderFunctor();
    auto source = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/diffusion_image.nii.gz"), &functor);
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(source));

    const auto path = m_TempDir + "/roundtrip.nii.gz";
    mitk::IOUtil::Save(source, "DWI_NIFTI", path);
    CPPUNIT_ASSERT(itksys::SystemTools::FileExists(m_TempDir + "/roundtrip.bvals"));
    CPPUNIT_ASSERT(itksys::SystemTools::FileExists(m_TempDir + "/roundtrip.bvecs"));
    auto image = mitk::IOUtil::Load<mitk::Image>(path, &functor);
    // bvals are written as integers, so gradient lengths come back rounded to the b-value.
    AssertSameDwi(source, image, 1e-6);
  }

  void PlainNrrdWriter_Dwi_WritesNoDiffusionKeys()
  {
    const auto path = m_TempDir + "/plain.nrrd";
    mitk::IOUtil::Save(m_Dwi, "application/vnd.mitk.image.nrrd", path);

    auto io = itk::NrrdImageIO::New();
    io->SetFileName(path);
    io->ReadImageInformation();
    for (const auto& key : io->GetMetaDataDictionary().GetKeys())
    {
      CPPUNIT_ASSERT_MESSAGE(key, key != "measurement frame");
      CPPUNIT_ASSERT_MESSAGE(key, key.find("DWMRI") == std::string::npos);
    }
  }

  void ReaderSelection_PlainImages_NoDwiReader()
  {
    CPPUNIT_ASSERT(!OffersDwiReader(GetTestDataFilePath("DiffusionImaging/tensorImage_FA.nrrd")));
    CPPUNIT_ASSERT(!OffersDwiReader(GetTestDataFilePath("DiffusionImaging/diffusionImageMask.nrrd")));
    CPPUNIT_ASSERT(!OffersDwiReader(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/mask_image.nii.gz")));

    // NIfTI file names with more dots and without bvals/bvecs.
    const auto dotted = m_TempDir + "/sub.01.nii.gz";
    itksys::SystemTools::CopyFileAlways(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/mask_image.nii.gz"), dotted);
    CPPUNIT_ASSERT(!OffersDwiReader(dotted));
  }

  void ReaderSelection_DiffusionFiles_DwiReaderIsDefault()
  {
    CPPUNIT_ASSERT_EQUAL(DwiReaderDescription(), DefaultReader(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi")));
    CPPUNIT_ASSERT_EQUAL(DwiReaderDescription(), DefaultReader(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/diffusion_image.nii.gz")));

    const auto nrrd = m_TempDir + "/dwi.nrrd";
    mitk::IOUtil::Save(m_Dwi, "DWI_NRRD", nrrd);
    CPPUNIT_ASSERT_EQUAL(DwiReaderDescription(), DefaultReader(nrrd));

    // HCP layout: data.nii.gz with "bvals" and "bvecs" in the same folder.
    const auto source = GetTestDataFilePath("DiffusionImaging/StreamlineTractography/diffusion_image");
    itksys::SystemTools::CopyFileAlways(source + ".nii.gz", m_TempDir + "/data.nii.gz");
    itksys::SystemTools::CopyFileAlways(source + ".bvals", m_TempDir + "/bvals");
    itksys::SystemTools::CopyFileAlways(source + ".bvecs", m_TempDir + "/bvecs");
    CPPUNIT_ASSERT_EQUAL(DwiReaderDescription(), DefaultReader(m_TempDir + "/data.nii.gz"));
  }

  void WriterSelection_Dwi_DwiWriterIsDefault()
  {
    mitk::IOUtil::Save(m_Dwi, m_TempDir + "/default.nii.gz");
    CPPUNIT_ASSERT(itksys::SystemTools::FileExists(m_TempDir + "/default.bvals"));
    CPPUNIT_ASSERT(itksys::SystemTools::FileExists(m_TempDir + "/default.bvecs"));
    auto functor = DwiReaderFunctor();
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(mitk::IOUtil::Load<mitk::Image>(m_TempDir + "/default.nii.gz", &functor)));

    // Plain images still go to MITK's writers.
    auto plain = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/tensorImage_FA.nrrd"));
    mitk::IOUtil::Save(plain, m_TempDir + "/plain.nrrd");
    CPPUNIT_ASSERT(!OffersDwiReader(m_TempDir + "/plain.nrrd"));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageIO)
