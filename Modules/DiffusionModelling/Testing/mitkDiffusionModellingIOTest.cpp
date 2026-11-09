/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkIOUtil.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkFileReaderSelector.h>
#include <mitkImageReadAccessor.h>
#include <mitkImageWriteAccessor.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkBaseDataSerializer.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionPeaks.h>

#include <mitkITKImageImport.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkMetaDataObject.h>
#include <itkVectorImage.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <itksys/SystemTools.hxx>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_determinant.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class mitkDiffusionModellingIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingIOTestSuite);
  MITK_TEST(ShReader_NoOption_ReadsMrtrix3FileConverted);
  MITK_TEST(TensorReader_StandardFrame_RotatedToWorld);
  MITK_TEST(TensorReader_MitkFrameKey_RotatedToWorld);
  MITK_TEST(TensorReader_PermutationFrame_RotatedToWorld);
  MITK_TEST(TensorReader_NoFrame_Unchanged);
  MITK_TEST(TensorReader_Nrrd_LeavesNoTemporaryFiles);
  MITK_TEST(TensorReader_Nifti_LeavesNoTemporaryFiles);
  MITK_TEST(Writers_UnwritableLocation_Throw);
  MITK_TEST(OdfReader_WrongVectorLength_Throws);
  MITK_TEST(Clone_KeepsTypeAndPixels);
  MITK_TEST(PeakImage_SetCustomColorWithoutPolydata);
  MITK_TEST(ShAndPeakIO_DotsInFileName);
  MITK_TEST(MimeTypes_ShAndPeakReadersOnlyForMatchingFiles);
  MITK_TEST(MimeTypes_ReplacedFile_CheckedAgain);
  MITK_TEST(ShWriter_NewNiftiFile_ReadsBackAsSh);
  MITK_TEST(ShIO_ConventionSurvivesSaving);
  MITK_TEST(PeakReader_MrtrixPeaks_ConvertedToLps);
  MITK_TEST(ShIO_OriginKeepsDoublePrecision);
  MITK_TEST(TensorReader_NiftiWrongComponentCount_Throws);
  MITK_TEST(Readers_RecordedContent);
  MITK_TEST(TensorAndOdfIO_RoundTrip);
  MITK_TEST(PeakImage_PolydataHasOneLinePerPeak);
  MITK_TEST(Serializers_RoundTrip);
  MITK_TEST(TensorReader_NiftiLayouts_ElementMapping);
  MITK_TEST(RgbCache_FollowsPixelChanges);
  MITK_TEST(RgbCache_FollowsShConvention);
  MITK_TEST(PeakImage_PolydataFollowsPixelChanges);
  MITK_TEST(ShReader_Mrtrix3_ConvertedToMitkBasis);
  MITK_TEST(ShReader_Dipy_ConvertedFromBvecFrame);
  MITK_TEST(PeakWriter_Nifti_RasWithMarker);
  MITK_TEST(ShReader_VectorLayout_OpensBothWays);
  MITK_TEST(ShReader_VectorLayout_MitkFslBasis);
  MITK_TEST(MimeTypes_VectorLayoutDwi_NoShReader);
  MITK_TEST(MimeTypes_Sh4dNifti_ImageReaderIsDefault);
  MITK_TEST(PeakReader_MrtrixKnownDirections_Lps);
  CPPUNIT_TEST_SUITE_END();

  using Matrix3 = vnl_matrix_fixed<double, 3, 3>;

  std::string m_TempDir;

  /** Writes a one-voxel ASCII NRRD tensor file (elements xx, xy, xz, yy, yz, zz) with the given header lines. */
  std::string WriteTensorNrrd(const std::string &name, const double (&elements)[6], const std::string &frameLines)
  {
    const std::string path = m_TempDir + "/" + name;
    std::ofstream file(path);
    file << "NRRD0004\n"
         << "type: float\n"
         << "dimension: 4\n"
         << "space: left-posterior-superior\n"
         << "sizes: 6 1 1 1\n"
         << "space directions: none (1,0,0) (0,1,0) (0,0,1)\n"
         << "kinds: 3D-symmetric-matrix domain domain domain\n"
         << "endian: little\n"
         << "encoding: ascii\n"
         << "space origin: (0,0,0)\n"
         << frameLines << (frameLines.empty() ? "" : "\n") << "\n";
    for (double e : elements)
      file << e << " ";
    file << "\n";
    return path;
  }

  Matrix3 ReadTensor(const std::string &path)
  {
    auto image = mitk::IOUtil::Load<mitk::TensorImage>(path);
    mitk::ImagePixelReadAccessor<mitk::TensorImage::PixelType, 3> access(image.GetPointer());
    itk::Index<3> idx = {{0, 0, 0}};
    auto tensor = access.GetPixelByIndex(idx);
    Matrix3 m;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        m(r, c) = tensor(r, c);
    return m;
  }

  static Matrix3 ToMatrix(const double (&e)[6])
  {
    Matrix3 m;
    m(0, 0) = e[0]; m(0, 1) = e[1]; m(0, 2) = e[2];
    m(1, 0) = e[1]; m(1, 1) = e[3]; m(1, 2) = e[4];
    m(2, 0) = e[2]; m(2, 1) = e[4]; m(2, 2) = e[5];
    return m;
  }

  static void AssertEqual(const Matrix3 &expected, const Matrix3 &actual)
  {
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected(r, c), actual(r, c), 1e-5);
  }

  /** The frame vectors are the rows of F; the world tensor is F^T T F, as for DWI gradients. */
  void CheckFrame(const std::string &name, const Matrix3 &frameRows, bool mitkKey)
  {
    const double elements[6] = {3, 0.2, 0, 1, 0, 0.5};
    std::ostringstream frame;
    frame.precision(17);
    if (mitkKey)
      frame << "measurement frame:=";
    else
      frame << "measurement frame: ";
    for (unsigned int r = 0; r < 3; ++r)
      frame << (mitkKey ? " ( " : "(") << frameRows(r, 0) << (mitkKey ? " , " : ",") << frameRows(r, 1) << (mitkKey ? " , " : ",")
            << frameRows(r, 2) << (mitkKey ? " )" : ") ");

    const Matrix3 expected = frameRows.transpose() * ToMatrix(elements) * frameRows;
    AssertEqual(expected, ReadTensor(WriteTensorNrrd(name, elements, frame.str())));
  }

  static Matrix3 RotationZ30()
  {
    const double c = std::cos(itk::Math::pi / 6);
    const double s = std::sin(itk::Math::pi / 6);
    Matrix3 f;
    f(0, 0) = c;  f(0, 1) = s; f(0, 2) = 0;
    f(1, 0) = -s; f(1, 1) = c; f(1, 2) = 0;
    f(2, 0) = 0;  f(2, 1) = 0; f(2, 2) = 1;
    return f;
  }

  /** The tensor reader tries non-NRRD files as NIfTI through a copy named temp_dti*. */
  static unsigned int CountTemporaryTensorCopies()
  {
    unsigned int count = 0;
    for (const auto &entry : std::filesystem::directory_iterator(mitk::IOUtil::GetTempPath()))
    {
      if (entry.path().filename().string().starts_with("temp_dti"))
        ++count;
    }
    return count;
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("mitkDiffusionModellingIOTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void TensorReader_StandardFrame_RotatedToWorld()
  {
    CheckFrame("standard.dti", RotationZ30(), false);
  }

  void TensorReader_MitkFrameKey_RotatedToWorld()
  {
    CheckFrame("mitk.dti", RotationZ30(), true);
  }

  void TensorReader_PermutationFrame_RotatedToWorld()
  {
    Matrix3 f;
    f.fill(0);
    f(0, 1) = 1;
    f(1, 2) = 1;
    f(2, 0) = 1;
    CheckFrame("permutation.dti", f, false);
  }

  void TensorReader_NoFrame_Unchanged()
  {
    const double elements[6] = {3, 0.2, 0, 1, 0, 0.5};
    AssertEqual(ToMatrix(elements), ReadTensor(WriteTensorNrrd("noframe.dti", elements, "")));
  }

  void TensorReader_Nrrd_LeavesNoTemporaryFiles()
  {
    const unsigned int before = CountTemporaryTensorCopies();
    mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti"));
    CPPUNIT_ASSERT_EQUAL(before, CountTemporaryTensorCopies());
  }

  void TensorReader_Nifti_LeavesNoTemporaryFiles()
  {
    using VectorImageType = itk::VectorImage<float, 3>;
    auto image = VectorImageType::New();
    VectorImageType::RegionType region;
    region.SetSize({{2, 1, 1}});
    image->SetRegions(region);
    image->SetVectorLength(6);
    image->Allocate();
    itk::VariableLengthVector<float> pixel(6);
    for (unsigned int i = 0; i < 6; ++i)
      pixel[i] = i + 1;
    image->FillBuffer(pixel);

    const std::string niftiPath = m_TempDir + "/tensor.nii.gz";
    mitk::IOUtil::Save(mitk::GrabItkImageMemory(image.GetPointer()), niftiPath);
    const std::string dtiPath = m_TempDir + "/tensor.dti";
    itksys::SystemTools::CopyAFile(niftiPath, dtiPath);

    const unsigned int before = CountTemporaryTensorCopies();
    const double elements[6] = {1, 2, 3, 4, 5, 6};
    AssertEqual(ToMatrix(elements), ReadTensor(dtiPath));
    CPPUNIT_ASSERT_EQUAL(before, CountTemporaryTensorCopies());
  }

  void Writers_UnwritableLocation_Throw()
  {
    const std::string missingDir = m_TempDir + "/missing/";

    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(tensorImage, missingDir + "tensor.dti"), std::exception);

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_QA0.qbi"));
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(odfImage, missingDir + "odf.odf"), std::exception);

    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &functor);
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(shImage, missingDir + "sh.shi"), std::exception);

    CPPUNIT_ASSERT(!itksys::SystemTools::FileExists(missingDir));
  }

  // The reader copied 252 values per voxel whatever the file held.
  void OdfReader_WrongVectorLength_Throws()
  {
    const std::string path = m_TempDir + "/short.odf";
    std::ofstream file(path);
    file << "NRRD0004\ntype: float\ndimension: 4\nspace: left-posterior-superior\nsizes: 251 1 1 1\n"
         << "space directions: none (1,0,0) (0,1,0) (0,0,1)\nkinds: vector domain domain domain\n"
         << "endian: little\nencoding: ascii\nspace origin: (0,0,0)\n\n";
    for (unsigned int i = 0; i < 251; ++i)
      file << "0.1 ";
    file << "\n";
    file.close();

    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load(path), std::exception);
  }

  /** Clones through the mitk::Image interface, as views do, and compares class and pixel buffer. */
  static void CheckClone(mitk::Image *image)
  {
    mitk::Image::Pointer clone = image->Clone();
    CPPUNIT_ASSERT(clone.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(std::string(image->GetNameOfClass()), std::string(clone->GetNameOfClass()));
    CPPUNIT_ASSERT(clone.GetPointer() != image);

    const auto size = image->GetPixelType().GetSize() * image->GetDimension(0) * image->GetDimension(1) * image->GetDimension(2);
    mitk::ImageReadAccessor original(image);
    mitk::ImageReadAccessor copy(clone);
    CPPUNIT_ASSERT(std::memcmp(original.GetData(), copy.GetData(), size) == 0);
  }

  // Clone() returned nullptr for these types, and cloning through mitk::Image lost the type.
  void Clone_KeepsTypeAndPixels()
  {
    CheckClone(mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti")));
    CheckClone(mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_QA0.qbi")));

    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    auto peakImage = mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &peakFunctor);
    CheckClone(peakImage);
    CPPUNIT_ASSERT(peakImage->Clone().IsNotNull());

    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    shImage->SetShConvention(mitk::ShImage::FSL);
    CheckClone(shImage);
    mitk::ShImage::Pointer shClone = shImage->Clone();
    CPPUNIT_ASSERT(shClone.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(mitk::ShImage::FSL, shClone->GetShConvention());
  }

  // A peak image created by a view has no polydata until it is rendered; colouring it crashed.
  void PeakImage_SetCustomColorWithoutPolydata()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    auto loaded = mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &functor);
    auto itkPeaks = loaded->GetItkPeakImage();

    auto peaks = mitk::PeakImage::New();
    peaks->InitializeByItk(itkPeaks.GetPointer());
    peaks->SetVolume(itkPeaks->GetBufferPointer());
    peaks->SetCustomColor(255, 0, 0);

    auto colors = dynamic_cast<vtkUnsignedCharArray *>(peaks->GetPolyData()->GetPointData()->GetArray("FIBER_COLORS"));
    CPPUNIT_ASSERT(colors != nullptr);
    CPPUNIT_ASSERT(colors->GetNumberOfTuples() > 0);
    CPPUNIT_ASSERT_EQUAL(peaks->GetPolyData()->GetNumberOfPoints(), colors->GetNumberOfTuples());
    CPPUNIT_ASSERT_EQUAL(255, static_cast<int>(colors->GetTypedComponent(0, 0)));
    CPPUNIT_ASSERT_EQUAL(0, static_cast<int>(colors->GetTypedComponent(0, 1)));
  }

  static void AssertEqualPixels(mitk::Image *expected, mitk::Image *actual)
  {
    CPPUNIT_ASSERT_EQUAL(std::string(expected->GetNameOfClass()), std::string(actual->GetNameOfClass()));
    const auto size = expected->GetPixelType().GetSize() * expected->GetDimension(0) * expected->GetDimension(1) * expected->GetDimension(2) * (expected->GetDimension() > 3 ? expected->GetDimension(3) : 1);
    CPPUNIT_ASSERT_EQUAL(size, actual->GetPixelType().GetSize() * actual->GetDimension(0) * actual->GetDimension(1) * actual->GetDimension(2) * (actual->GetDimension() > 3 ? actual->GetDimension(3) : 1));
    mitk::ImageReadAccessor expectedAccess(expected);
    mitk::ImageReadAccessor actualAccess(actual);
    CPPUNIT_ASSERT(std::memcmp(expectedAccess.GetData(), actualAccess.GetData(), size) == 0);
  }

  /**
   * SH files hold MITK's coefficients in MRtrix 3's basis (m != 0 divided by sqrt(2)), so a saved
   * and reloaded image agrees to float rounding, not bit for bit.
   */
  static void AssertEqualShCoefficients(mitk::Image *expected, mitk::Image *actual)
  {
    CPPUNIT_ASSERT_EQUAL(std::string(expected->GetNameOfClass()), std::string(actual->GetNameOfClass()));
    const auto count = expected->GetPixelType().GetNumberOfComponents() * expected->GetDimension(0) * expected->GetDimension(1) * expected->GetDimension(2);
    CPPUNIT_ASSERT_EQUAL(count, actual->GetPixelType().GetNumberOfComponents() * actual->GetDimension(0) * actual->GetDimension(1) * actual->GetDimension(2));
    mitk::ImageReadAccessor expectedAccess(expected);
    mitk::ImageReadAccessor actualAccess(actual);
    const auto *expectedValues = static_cast<const float *>(expectedAccess.GetData());
    const auto *actualValues = static_cast<const float *>(actualAccess.GetData());
    for (std::size_t i = 0; i < count; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedValues[i], actualValues[i], 3e-7 * std::fabs(expectedValues[i]));
  }

  // The readers and the SH writer chose their ITK IO by everything after the first dot.
  void ShAndPeakIO_DotsInFileName()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    const std::string shPath = m_TempDir + "/sh.v2.shi";
    mitk::IOUtil::Save(shImage, shPath);
    AssertEqualShCoefficients(shImage, mitk::IOUtil::Load<mitk::ShImage>(shPath, &shFunctor));

    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    auto peakImage = mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &peakFunctor);
    const std::string peakPath = m_TempDir + "/peaks.v2.nrrd";
    mitk::IOUtil::Save(peakImage, peakPath);
    AssertEqualPixels(peakImage, mitk::IOUtil::Load<mitk::PeakImage>(peakPath, &peakFunctor));
  }

  static std::set<std::string> ReaderDescriptions(const std::string &path)
  {
    std::set<std::string> descriptions;
    for (const auto &item : mitk::FileReaderSelector(path).Get())
      descriptions.insert(item.GetDescription());
    return descriptions;
  }

  // The SH and peak mime types were registered as plain copies without their content checks, so
  // both readers were offered for every NIfTI and NRRD file.
  void MimeTypes_ShAndPeakReadersOnlyForMatchingFiles()
  {
    for (const std::string file : {"mean_FA.nii.gz", "tensorImage_FA.nrrd"})
    {
      const auto descriptions = ReaderDescriptions(mitk::GetDiffusionTestDataFilePath(file));
      CPPUNIT_ASSERT_MESSAGE(file, descriptions.count("SH Image") == 0);
      CPPUNIT_ASSERT_MESSAGE(file, descriptions.count("Peak Image") == 0);
    }

    const auto shDescriptions = ReaderDescriptions(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"));
    CPPUNIT_ASSERT(shDescriptions.count("SH Image") == 1);
    CPPUNIT_ASSERT(shDescriptions.count("Peak Image") == 1);
  }

  // The SH and peak mime types share the header of the last file checked; a file replaced under
  // the same name is read again.
  void MimeTypes_ReplacedFile_CheckedAgain()
  {
    const std::string path = m_TempDir + "/replaced.nii.gz";
    std::filesystem::copy_file(mitk::GetDiffusionTestDataFilePath("mean_FA.nii.gz"), path);
    CPPUNIT_ASSERT(ReaderDescriptions(path).count("SH Image") == 0);
    std::filesystem::copy_file(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), path, std::filesystem::copy_options::overwrite_existing);
    const auto descriptions = ReaderDescriptions(path);
    CPPUNIT_ASSERT(descriptions.count("SH Image") == 1);
    CPPUNIT_ASSERT(descriptions.count("Peak Image") == 1);
  }

  // Writers are chosen through the mime type of the output path, which does not exist yet.
  void ShWriter_NewNiftiFile_ReadsBackAsSh()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    const std::string path = m_TempDir + "/sh.nii.gz";
    mitk::IOUtil::Save(shImage, path);
    AssertEqualShCoefficients(shImage, mitk::IOUtil::Load<mitk::ShImage>(path, &shFunctor));
  }

  // An image read as FSL/Dipy SH was saved without its convention and read back as MRtrix.
  void ShIO_ConventionSurvivesSaving()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    for (auto convention : {mitk::ShImage::FSL, mitk::ShImage::MRTRIX})
    {
      shImage->SetShConvention(convention);
      for (const std::string extension : {".nii.gz", ".shi", ".nrrd"})
      {
        const std::string path = m_TempDir + "/convention" + std::to_string(convention) + extension;
        mitk::IOUtil::Save(shImage, path);
        auto reloaded = mitk::IOUtil::Load<mitk::ShImage>(path, &shFunctor);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(path, convention, reloaded->GetShConvention());
        AssertEqualShCoefficients(shImage, reloaded);
      }
    }
  }

  static mitk::PeakImage::Pointer LoadPeaks(const std::string &path, std::optional<bool> convertMrtrix)
  {
    mitk::IFileReader::Options options;
    if (convertMrtrix)
      options["Convert MRtrix peaks (RAS) to MITK (LPS)"] = *convertMrtrix;
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, options);
    return mitk::IOUtil::Load<mitk::PeakImage>(path, &functor);
  }

  static std::vector<float> Values(mitk::PeakImage *peaks)
  {
    auto itkPeaks = peaks->GetItkPeakImage();
    return std::vector<float>(itkPeaks->GetBufferPointer(), itkPeaks->GetBufferPointer() + itkPeaks->GetLargestPossibleRegion().GetNumberOfPixels());
  }

  // MRtrix writes RAS directions; MITK's peaks are LPS.
  void PeakReader_MrtrixPeaks_ConvertedToLps()
  {
    const std::string mrtrixFile = mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_peak_image.nii.gz");
    auto raw = LoadPeaks(mrtrixFile, false);
    auto converted = LoadPeaks(mrtrixFile, true);

    const auto rawValues = Values(raw);
    const auto convertedValues = Values(converted);
    CPPUNIT_ASSERT_EQUAL(rawValues.size(), convertedValues.size());
    const unsigned int numVoxels = raw->GetDimension(0) * raw->GetDimension(1) * raw->GetDimension(2);
    unsigned int nonZero = 0;
    for (std::size_t i = 0; i < rawValues.size(); ++i)
    {
      // The file holds NaN for missing peaks; niftilib reads them as 0.
      CPPUNIT_ASSERT(std::isfinite(rawValues[i]));
      const bool isZ = (i / numVoxels) % 3 == 2;
      CPPUNIT_ASSERT_EQUAL(isZ ? rawValues[i] : -rawValues[i], convertedValues[i]);
      if (rawValues[i] != 0)
        ++nonZero;
    }
    CPPUNIT_ASSERT(nonZero > 0);
    CPPUNIT_ASSERT(Values(LoadPeaks(mrtrixFile, std::nullopt)) == convertedValues);

    vtkPoints *points = raw->GetPolyData()->GetPoints();
    CPPUNIT_ASSERT(points->GetNumberOfPoints() > 0);
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
      CPPUNIT_ASSERT(std::isfinite(points->GetPoint(i)[0]) && std::isfinite(points->GetPoint(i)[1]) && std::isfinite(points->GetPoint(i)[2]));

    // Peak files without the MRtrix description, e.g. NIfTI files of earlier MITK versions, read unchanged.
    const std::string mitkFile = mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz");
    CPPUNIT_ASSERT(Values(LoadPeaks(mitkFile, false)) == Values(LoadPeaks(mitkFile, true)));
  }

  // The SH importer and exporter converted origins through float.
  void ShIO_OriginKeepsDoublePrecision()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    mitk::Point3D origin;
    origin[0] = 0.1;
    origin[1] = -20.3;
    origin[2] = 1000.7;
    shImage->GetGeometry()->SetOrigin(origin);

    const std::string path = m_TempDir + "/origin.shi";
    mitk::IOUtil::Save(shImage, path);
    auto reloaded = mitk::IOUtil::Load<mitk::ShImage>(path, &shFunctor);
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(origin[i], reloaded->GetGeometry()->GetOrigin()[i], 1e-12);
  }

  // A NIfTI file with neither 6 nor 9 tensor components gave an uninitialized tensor image.
  void TensorReader_NiftiWrongComponentCount_Throws()
  {
    using VectorImageType = itk::VectorImage<float, 3>;
    auto image = VectorImageType::New();
    VectorImageType::RegionType region;
    region.SetSize({{2, 1, 1}});
    image->SetRegions(region);
    image->SetVectorLength(4);
    image->Allocate();
    itk::VariableLengthVector<float> pixel(4);
    pixel.Fill(1);
    image->FillBuffer(pixel);

    const std::string niftiPath = m_TempDir + "/four.nii.gz";
    mitk::IOUtil::Save(mitk::GrabItkImageMemory(image.GetPointer()), niftiPath);
    const std::string dtiPath = m_TempDir + "/four.dti";
    itksys::SystemTools::CopyAFile(niftiPath, dtiPath);

    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load(dtiPath), std::exception);
  }

  /** FNV-1a hash of the whole pixel buffer of an image. */
  static std::uint64_t Hash(mitk::Image *image)
  {
    std::size_t size = image->GetPixelType().GetSize();
    for (unsigned int d = 0; d < image->GetDimension(); ++d)
      size *= image->GetDimension(d);
    mitk::ImageReadAccessor access(image);
    auto bytes = static_cast<const unsigned char *>(access.GetData());
    std::uint64_t hash = 14695981039346656037ull;
    for (std::size_t i = 0; i < size; ++i)
      hash = (hash ^ bytes[i]) * 1099511628211ull;
    return hash;
  }

  static void AssertEqualGeometry(mitk::Image *expected, mitk::Image *actual)
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(expected->GetDimension(i), actual->GetDimension(i));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetGeometry()->GetOrigin()[i], actual->GetGeometry()->GetOrigin()[i], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetGeometry()->GetSpacing()[i], actual->GetGeometry()->GetSpacing()[i], 1e-6);
      for (unsigned int j = 0; j < 3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetGeometry()->GetIndexToWorldTransform()->GetMatrix()[i][j],
                                     actual->GetGeometry()->GetIndexToWorldTransform()->GetMatrix()[i][j], 1e-6);
    }
  }

  // Pixel buffers and geometry of mitk-diffusion-test-data files as read before the readers were restructured.
  void Readers_RecordedContent()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    // The hash pins the file's raw coefficients; ShReader_Mrtrix3_ConvertedToMitkBasis checks
    // that this MRtrix 3 file is converted by default.
    mitk::IFileReader::Options legacyOptions;
    legacyOptions["SH convention"] = std::string("MITK legacy (MRtrix 0.2)");
    mitk::PreferenceListReaderOptionsFunctor legacyShFunctor({"SH Image"}, legacyOptions);
    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    struct Recorded
    {
      mitk::Image::Pointer image;
      std::uint64_t hash;
      double origin[3];
      double spacing[3];
    };
    const Recorded recorded[] = {
      {mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti")).GetPointer(),
       910238346892541248ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("FWE/ReferenceTensorImage.dti")).GetPointer(),
       14634070166575176342ull, {-14.1639, 13.4743, -12.6538}, {2, 2, 2}},
      {mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_QA0.qbi")).GetPointer(),
       17363730383617297091ull, {1, 1, 1}, {2, 2, 2}},
      {mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_image.qbi")).GetPointer(),
       13068465665225140656ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor).GetPointer(),
       13341037711982211669ull, {0.615, 1.0005, 0.7045}, {1.23, 2.001, 1.409}},
      {mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_image_sh.nii.gz"), &legacyShFunctor).GetPointer(),
       8054354499704950156ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &peakFunctor).GetPointer(),
       1038741785062553984ull, {0.615, 1.0005, 0.7045}, {1.23, 2.001, 1.409}},
      {mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("FiberFit/csd_peak_image.nii.gz"), &peakFunctor).GetPointer(),
       14391242496329318739ull, {0, 0, 0}, {3, 3, 3}}};
    for (const auto &[image, hash, origin, spacing] : recorded)
    {
      CPPUNIT_ASSERT_EQUAL(hash, Hash(image));
      for (unsigned int i = 0; i < 3; ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(origin[i], image->GetGeometry()->GetOrigin()[i], 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(spacing[i], image->GetGeometry()->GetSpacing()[i], 1e-6);
      }
    }
  }

  // Tensor and ODF images survive saving in each of their formats bit for bit.
  void TensorAndOdfIO_RoundTrip()
  {
    for (const std::string file : {"tensorImage.dti", "FWE/ReferenceTensorImage.dti"})
    {
      auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath(file));
      const std::string path = m_TempDir + "/tensor.dti";
      mitk::IOUtil::Save(tensorImage, path);
      auto reloaded = mitk::IOUtil::Load<mitk::TensorImage>(path);
      AssertEqualPixels(tensorImage, reloaded);
      AssertEqualGeometry(tensorImage, reloaded);
    }

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_QA0.qbi"));
    for (const std::string extension : {".odf", ".qbi"})
    {
      const std::string path = m_TempDir + "/odf" + extension;
      mitk::IOUtil::Save(odfImage, path);
      auto reloaded = mitk::IOUtil::Load<mitk::OdfImage>(path);
      AssertEqualPixels(odfImage, reloaded);
      AssertEqualGeometry(odfImage, reloaded);
    }
  }

  // The peak polydata holds one two-point line per peak longer than 0.001, with a colour per point.
  void PeakImage_PolydataHasOneLinePerPeak()
  {
    for (const std::string file : {"peak_image_test.nii.gz", "FiberFit/csd_peak_image.nii.gz"})
    {
      auto peakImage = LoadPeaks(mitk::GetDiffusionTestDataFilePath(file), std::nullopt);
      auto itkPeaks = peakImage->GetItkPeakImage();
      const auto size = itkPeaks->GetLargestPossibleRegion().GetSize();
      vtkIdType expectedLines = 0;
      for (itk::IndexValueType z = 0; z < static_cast<itk::IndexValueType>(size[2]); ++z)
        for (itk::IndexValueType y = 0; y < static_cast<itk::IndexValueType>(size[1]); ++y)
          for (itk::IndexValueType x = 0; x < static_cast<itk::IndexValueType>(size[0]); ++x)
            for (itk::IndexValueType p = 0; p < static_cast<itk::IndexValueType>(size[3]); p += 3)
            {
              double squaredLength = 0;
              for (itk::IndexValueType c = 0; c < 3; ++c)
              {
                const double v = itkPeaks->GetPixel({{x, y, z, p + c}});
                squaredLength += v * v;
              }
              if (std::sqrt(squaredLength) >= 0.001)
                ++expectedLines;
            }

      auto polyData = peakImage->GetPolyData();
      CPPUNIT_ASSERT_MESSAGE(file, expectedLines > 0);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(file, expectedLines, polyData->GetNumberOfLines());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(file, 2 * expectedLines, polyData->GetNumberOfPoints());
      auto colors = polyData->GetPointData()->GetArray("FIBER_COLORS");
      CPPUNIT_ASSERT(colors != nullptr);
      CPPUNIT_ASSERT_EQUAL(polyData->GetNumberOfPoints(), colors->GetNumberOfTuples());
    }
  }

  /** Runs the scene serializer registered for the data's class into the temp directory. */
  std::string Serialize(mitk::BaseData *data, const std::string &hint)
  {
    const std::string serializerName = std::string(data->GetNameOfClass()) + "Serializer";
    auto candidates = itk::ObjectFactoryBase::CreateAllInstance(serializerName.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE(serializerName, std::size_t(1), candidates.size());
    auto serializer = dynamic_cast<mitk::BaseDataSerializer *>(candidates.front().GetPointer());
    CPPUNIT_ASSERT(serializer != nullptr);
    serializer->SetData(data);
    serializer->SetFilenameHint(hint);
    serializer->SetWorkingDirectory(m_TempDir);
    const std::string filename = serializer->Serialize();
    CPPUNIT_ASSERT_MESSAGE(serializerName, !filename.empty());
    return m_TempDir + "/" + filename;
  }

  // Scene files keep tensor, ODF, SH and peak images with their type and pixels.
  void Serializers_RoundTrip()
  {
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    auto tensorReloaded = mitk::IOUtil::Load<mitk::TensorImage>(Serialize(tensorImage, "tensor"));
    AssertEqualPixels(tensorImage, tensorReloaded);
    AssertEqualGeometry(tensorImage, tensorReloaded);

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_QA0.qbi"));
    auto odfReloaded = mitk::IOUtil::Load<mitk::OdfImage>(Serialize(odfImage, "odf"));
    AssertEqualPixels(odfImage, odfReloaded);
    AssertEqualGeometry(odfImage, odfReloaded);

    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &shFunctor);
    auto shReloaded = mitk::IOUtil::Load<mitk::ShImage>(Serialize(shImage, "sh"), &shFunctor);
    AssertEqualShCoefficients(shImage, shReloaded);
    AssertEqualGeometry(shImage, shReloaded);

    auto peakImage = LoadPeaks(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), std::nullopt);
    auto peakReloaded = LoadPeaks(Serialize(peakImage, "peaks"), std::nullopt);
    AssertEqualPixels(peakImage, peakReloaded);
    AssertEqualGeometry(peakImage, peakReloaded);
  }

  /** Writes the image with ITK as NIfTI and renames the file to .dti, as users do with FSL tensors. */
  template <class TImage>
  std::string WriteNiftiAsDti(TImage *image, const std::string &name)
  {
    const std::string niftiPath = m_TempDir + "/" + name + ".nii.gz";
    auto writer = itk::ImageFileWriter<TImage>::New();
    writer->SetInput(image);
    writer->SetFileName(niftiPath);
    writer->Update();
    const std::string dtiPath = m_TempDir + "/" + name + ".dti";
    itksys::SystemTools::CopyAFile(niftiPath, dtiPath);
    return dtiPath;
  }

  // NIfTI tensors: 6 components xx, xy, xz, yy, yz, zz; 9 components a full row-major matrix;
  // or the same as 6 or 9 volumes of a 4D image.
  void TensorReader_NiftiLayouts_ElementMapping()
  {
    const double fromSix[6] = {1, 2, 3, 4, 5, 6};
    const double fromNine[6] = {1, 2, 3, 5, 6, 9};

    using VectorImageType = itk::VectorImage<float, 3>;
    auto vectorImage = VectorImageType::New();
    VectorImageType::RegionType region;
    region.SetSize({{2, 1, 1}});
    vectorImage->SetRegions(region);
    vectorImage->SetVectorLength(9);
    vectorImage->Allocate();
    itk::VariableLengthVector<float> pixel(9);
    for (unsigned int i = 0; i < 9; ++i)
      pixel[i] = i + 1;
    vectorImage->FillBuffer(pixel);
    AssertEqual(ToMatrix(fromNine), ReadTensor(WriteNiftiAsDti(vectorImage.GetPointer(), "nine")));

    using VolumeImageType = itk::Image<float, 4>;
    for (unsigned int numVolumes : {6u, 9u})
    {
      auto volumes = VolumeImageType::New();
      VolumeImageType::RegionType volumeRegion;
      volumeRegion.SetSize({{2, 1, 1, numVolumes}});
      volumes->SetRegions(volumeRegion);
      volumes->Allocate();
      for (itk::IndexValueType v = 0; v < numVolumes; ++v)
        for (itk::IndexValueType x = 0; x < 2; ++x)
          volumes->SetPixel({{x, 0, 0, v}}, v + 1);
      AssertEqual(ToMatrix(numVolumes == 6 ? fromSix : fromNine),
                  ReadTensor(WriteNiftiAsDti(volumes.GetPointer(), "volumes" + std::to_string(numVolumes))));
    }
  }

  /** RGBA values of a colour map as shown by MITK's image mappers. */
  static std::vector<unsigned char> Colours(mitk::Image *image)
  {
    vtkImageData *rgb = image->GetVtkImageData();
    CPPUNIT_ASSERT(rgb != nullptr);
    CPPUNIT_ASSERT_EQUAL(4, rgb->GetNumberOfScalarComponents());
    const auto *begin = static_cast<const unsigned char *>(rgb->GetScalarPointer());
    return std::vector<unsigned char>(begin, begin + 4 * rgb->GetNumberOfPoints());
  }

  // The colour map was built once and kept when the pixels changed later.
  void RgbCache_FollowsPixelChanges()
  {
    auto image = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    const auto before = Colours(image);
    CPPUNIT_ASSERT(image->GetVtkImageData() == image->GetVtkImageData());

    {
      // Every tensor diag(1, 0, 0): FA 1, principal direction x, so opaque red.
      mitk::ImageWriteAccessor accessor(image);
      auto *tensors = static_cast<float *>(accessor.GetData());
      const std::size_t numberOfVoxels = image->GetDimension(0) * image->GetDimension(1) * image->GetDimension(2);
      for (std::size_t v = 0; v < numberOfVoxels; ++v)
      {
        std::fill(tensors + 6 * v, tensors + 6 * v + 6, 0.0f);
        tensors[6 * v] = 1.0f;
      }
    }
    image->Modified();

    const auto after = Colours(image);
    CPPUNIT_ASSERT(before != after);
    for (std::size_t i = 0; i < after.size(); i += 4)
    {
      CPPUNIT_ASSERT(after[i] >= 254 && after[i + 3] >= 254);
      CPPUNIT_ASSERT_EQUAL(0, int(after[i + 1]));
      CPPUNIT_ASSERT_EQUAL(0, int(after[i + 2]));
    }
  }

  // The colours of an SH image depend on its convention.
  void RgbCache_FollowsShConvention()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto image = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &functor);
    CPPUNIT_ASSERT_EQUAL(mitk::ShImage::MRTRIX, image->GetShConvention());
    const auto mrtrix = Colours(image);

    image->SetShConvention(mitk::ShImage::FSL);
    const auto fsl = Colours(image);
    CPPUNIT_ASSERT(mrtrix != fsl);

    auto reference = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &functor);
    reference->SetShConvention(mitk::ShImage::FSL);
    CPPUNIT_ASSERT(fsl == Colours(reference));
  }

  // The peak lines were built once and kept when the peaks changed later; a rebuilt polydata
  // keeps the custom colour.
  void PeakImage_PolydataFollowsPixelChanges()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    auto peaks = mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &functor);
    vtkSmartPointer<vtkPolyData> before = peaks->GetPolyData();
    CPPUNIT_ASSERT(before->GetNumberOfLines() > 1);
    peaks->SetCustomColor(0, 255, 0);
    CPPUNIT_ASSERT(before == peaks->GetPolyData());

    {
      // Only the first peak of the first voxel stays.
      mitk::ImageWriteAccessor accessor(peaks);
      auto *values = static_cast<float *>(accessor.GetData());
      const std::size_t numberOfValues = peaks->GetDimension(0) * peaks->GetDimension(1) * peaks->GetDimension(2) * peaks->GetDimension(3);
      std::fill(values, values + numberOfValues, 0.0f);
      values[0] = 1.0f;
    }
    peaks->Modified();

    vtkSmartPointer<vtkPolyData> after = peaks->GetPolyData();
    CPPUNIT_ASSERT(before != after);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(1), after->GetNumberOfLines());
    auto colors = dynamic_cast<vtkUnsignedCharArray *>(after->GetPointData()->GetArray("FIBER_COLORS"));
    CPPUNIT_ASSERT(colors != nullptr);
    CPPUNIT_ASSERT_EQUAL(0, static_cast<int>(colors->GetTypedComponent(0, 0)));
    CPPUNIT_ASSERT_EQUAL(255, static_cast<int>(colors->GetTypedComponent(0, 1)));
  }

  static mitk::ShImage::Pointer LoadSh(const std::string &path, const std::string &convention = "auto (file marker, else MITK legacy; vector layout: MITK FSL)")
  {
    mitk::IFileReader::Options options;
    options["SH convention"] = convention;
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, options);
    return mitk::IOUtil::Load<mitk::ShImage>(path, &functor);
  }

  /** Coefficients of voxel (x, y, z) of an SH image. */
  static vnl_vector<float> ShCoefficients(mitk::ShImage *image, unsigned int x, unsigned int y, unsigned int z)
  {
    const unsigned int count = image->GetPixelType().GetNumberOfComponents();
    mitk::ImageReadAccessor access(image);
    const auto *values = static_cast<const float *>(access.GetData());
    const std::size_t voxel = x + image->GetDimension(0) * (y + static_cast<std::size_t>(image->GetDimension(1)) * z);
    return vnl_vector<float>(values + voxel * count, count);
  }

  /** Writes a one-voxel SH file in MRtrix's layout with the given coefficients, NIfTI descrip and direction. */
  std::string WriteShNifti(const std::string &name, const std::vector<float> &coefficients, const std::string &descrip, const Matrix3 &direction)
  {
    using ImageType = mitk::ShImage::ShOnDiskType;
    auto image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize({{1, 1, 1, coefficients.size()}});
    image->SetRegions(region);
    ImageType::DirectionType imageDirection;
    imageDirection.SetIdentity();
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        imageDirection[r][c] = direction(r, c);
    image->SetDirection(imageDirection);
    image->Allocate();
    for (unsigned int j = 0; j < coefficients.size(); ++j)
      image->SetPixel({{0, 0, 0, j}}, coefficients[j]);
    if (!descrip.empty())
      itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "ITK_FileNotes", descrip);

    const std::string path = m_TempDir + "/" + name;
    auto writer = itk::ImageFileWriter<ImageType>::New();
    writer->SetInput(image);
    writer->SetFileName(path);
    writer->Update();
    return path;
  }

  // The reader is constructed while the library is loaded; its option labels must be valid then,
  // or the default "auto" is empty and MRtrix 3 files are not converted.
  void ShReader_NoOption_ReadsMrtrix3FileConverted()
  {
    const std::string path = mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_image_sh.nii.gz");
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto csd = mitk::IOUtil::Load<mitk::ShImage>(path, &functor);
    // MRtrix 3.0.8 sh2amp at voxel (30, 30, 1), RAS (-1, 0, 0).
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.20321192, mitk::sh::GetValue(ShCoefficients(csd, 30, 30, 1), 8, vnl_vector_fixed<double, 3>(1, 0, 0), true), 1e-6);
  }

  // MITK's basis is MRtrix 0.2's; MRtrix 3 (every FOD since 3.0) uses an orthonormal basis whose
  // functions with m != 0 are sqrt(2) times MITK's. Reading MRtrix 3 files as MITK's puts 23 % of
  // the ODF maximum into wrong directions, and MRtrix reads MITK's legacy files wrongly.
  void ShReader_Mrtrix3_ConvertedToMitkBasis()
  {
    const vnl_vector_fixed<double, 3> x(1, 0, 0);
    Matrix3 identity;
    identity.set_identity();

    // Only (l = 2, m = 2): MRtrix 3 gives sqrt(15 / (16 pi)) at world (1, 0, 0), i.e. RAS (-1, 0, 0).
    const double mrtrix3Value = std::sqrt(15.0 / (16.0 * itk::Math::pi));
    std::vector<float> coefficients(6, 0);
    coefficients[5] = 1;
    const auto marked = WriteShNifti("mrtrix3.nii.gz", coefficients, "MRtrix version: 3.0.8", identity);
    const auto unmarked = WriteShNifti("unmarked.nii.gz", coefficients, "", identity);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrix3Value, mitk::sh::GetValue(ShCoefficients(LoadSh(marked), 0, 0, 0), 2, x, true), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrix3Value, mitk::sh::GetValue(ShCoefficients(LoadSh(unmarked, "MRtrix 3"), 0, 0, 0), 2, x, true), 1e-6);
    // Files without a marker, which includes every SH file of earlier MITK versions, are read in the legacy basis.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrix3Value / std::sqrt(2.0), mitk::sh::GetValue(ShCoefficients(LoadSh(unmarked), 0, 0, 0), 2, x, true), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrix3Value / std::sqrt(2.0), mitk::sh::GetValue(ShCoefficients(LoadSh(marked, "MITK legacy (MRtrix 0.2)"), 0, 0, 0), 2, x, true), 1e-6);

    // MITK writes MRtrix 3's basis with its own marker, which decides over the option.
    for (const std::string extension : {".nii.gz", ".shi"})
    {
      const std::string path = m_TempDir + "/written" + extension;
      mitk::IOUtil::Save(LoadSh(marked), path);
      // ITK recognizes a NRRD file by its extension.
      const std::string rawPath = extension == ".shi" ? path + ".nrrd" : path;
      if (extension == ".shi")
        std::filesystem::copy_file(path, rawPath);
      auto reader = itk::ImageFileReader<mitk::ShImage::ShOnDiskType>::New();
      reader->SetFileName(rawPath);
      reader->Update();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, reader->GetOutput()->GetPixel({{0, 0, 0, 5}}), 1e-7);
      std::string marker;
      if (extension == ".shi")
        itk::ExposeMetaData<std::string>(reader->GetOutput()->GetMetaDataDictionary(), "MITK_SH_convention", marker);
      else
        itk::ExposeMetaData<std::string>(reader->GetOutput()->GetMetaDataDictionary(), "ITK_FileNotes", marker);
      CPPUNIT_ASSERT_EQUAL(std::string(extension == ".shi" ? "MRtrix3" : "MITK_SH_convention=MRtrix3"), marker);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrix3Value, mitk::sh::GetValue(ShCoefficients(LoadSh(path, "MITK legacy (MRtrix 0.2)"), 0, 0, 0), 2, x, true), 1e-6);
    }

    // The MRtrix 3 FOD in mitk-diffusion-test-data, voxel (30, 30, 1): amplitudes of MRtrix 3.0.8 sh2amp at
    // RAS (-1, 0, 0), (0, -1, 0), (0, 0, 1), (-0.48, 0.6, 0.64), given here as LPS directions.
    auto csd = LoadSh(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_image_sh.nii.gz"));
    const vnl_vector<float> csdCoefficients = ShCoefficients(csd, 30, 30, 1);
    const double directions[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.48, -0.6, 0.64}};
    const double amplitudes[4] = {0.20321192, -0.00769596, -0.007680845, 0.023118196};
    for (unsigned int i = 0; i < 4; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(amplitudes[i], mitk::sh::GetValue(csdCoefficients, 8, vnl_vector_fixed<double, 3>(directions[i]), true), 1e-6);
  }

  // "FSL/Dipy" read the coefficients in MITK's own FSL definition at world directions, which is
  // neither dipy's default basis nor the frame of dipy's coefficients (the bvecs it was given):
  // the first peaks of a dipy CSA fit were 47 degrees off (median).
  void ShReader_Dipy_ConvertedFromBvecFrame()
  {
    // dipy 1.12.1 sf_to_sh(..., sh_order_max=4, basis_type='descoteaux07', legacy=True) of
    // f(v) = (v.d)^4 on repulsion724, d = (1, 2, 3) / |(1, 2, 3)| in the bvec frame.
    const std::vector<float> dipyCoefficients = {0.70898154f, -0.168114716f, -0.336229432f, 0.420597999f, -0.672458865f,
      0.224152955f, -0.00713331864f, 0.0951157676f, -0.113237923f, -0.0686325455f, -0.0614932969f, -0.137265091f,
      0.150983897f, 0.0172937759f, -0.0244570925f};
    vnl_vector_fixed<double, 3> d(1, 2, 3);
    d.normalize();

    // LAS (dcm2niix's layout, negative determinant: bvec frame = voxel axes), identity and an
    // oblique matrix (positive determinant: x of the voxel frame negated).
    Matrix3 las;
    las.set_identity();
    las(1, 1) = -1;
    Matrix3 identity;
    identity.set_identity();
    const double angle = 0.4;
    Matrix3 oblique;
    oblique.set_identity();
    oblique(0, 0) = std::cos(angle);
    oblique(0, 1) = -std::sin(angle);
    oblique(1, 0) = std::sin(angle);
    oblique(1, 1) = std::cos(angle);

    const vnl_matrix<double> directions = SphereDirections(200);
    unsigned int index = 0;
    for (const Matrix3 &direction : {las, identity, oblique})
    {
      Matrix3 flip;
      flip.set_identity();
      if (vnl_det(direction) > 0)
        flip(0, 0) = -1;
      const Matrix3 frame = (direction * flip).transpose(); // v = frame * u

      const auto path = WriteShNifti("dipy" + std::to_string(index++) + ".nii.gz", dipyCoefficients, "", direction);
      auto image = LoadSh(path, "FSL/Dipy");
      CPPUNIT_ASSERT_EQUAL(mitk::ShImage::MRTRIX, image->GetShConvention());
      const vnl_vector<float> coefficients = ShCoefficients(image, 0, 0, 0);
      for (unsigned int i = 0; i < directions.cols(); ++i)
      {
        const vnl_vector_fixed<double, 3> u = directions.get_column(i);
        const double expected = std::pow(dot_product(frame * u, d), 4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, mitk::sh::GetValue(coefficients, 4, u, true), 1e-5);
      }
    }

    // Files that MITK marked with its own FSL definition keep it.
    auto image = LoadSh(WriteShNifti("marked.nii.gz", dipyCoefficients, "MITK_SH_convention=FSL", identity), "FSL/Dipy");
    CPPUNIT_ASSERT_EQUAL(mitk::ShImage::FSL, image->GetShConvention());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dipyCoefficients[1], ShCoefficients(image, 0, 0, 0)[1], 1e-9);
  }

  /** Unit vectors spread evenly over the sphere (Fibonacci lattice), as columns. */
  static vnl_matrix<double> SphereDirections(unsigned int count)
  {
    vnl_matrix<double> directions(3, count);
    const double golden = itk::Math::pi * (3.0 - std::sqrt(5.0));
    for (unsigned int i = 0; i < count; ++i)
    {
      const double z = 1.0 - (2.0 * i + 1.0) / count;
      const double r = std::sqrt(1.0 - z * z);
      directions(0, i) = r * std::cos(golden * i);
      directions(1, i) = r * std::sin(golden * i);
      directions(2, i) = z;
    }
    return directions;
  }

  // MITK wrote NIfTI peaks with its LPS vectors and no marker, so MRtrix and the tools that use its
  // peak files read them mirrored in x and y.
  void PeakWriter_Nifti_RasWithMarker()
  {
    auto peaks = LoadPeaks(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), std::nullopt);
    const auto lps = Values(peaks);
    const unsigned int numVoxels = peaks->GetDimension(0) * peaks->GetDimension(1) * peaks->GetDimension(2);
    unsigned int nonZero = 0;
    for (const auto value : lps)
      nonZero += value != 0;
    CPPUNIT_ASSERT(nonZero > 0);

    for (const std::string extension : {".nii.gz", ".nrrd"})
    {
      const std::string path = m_TempDir + "/written_peaks" + extension;
      mitk::IOUtil::Save(peaks, path);
      CPPUNIT_ASSERT(Values(peaks) == lps);

      auto reader = itk::ImageFileReader<mitk::PeakImage::ItkPeakImageType>::New();
      reader->SetFileName(path);
      reader->Update();
      const auto *raw = reader->GetOutput()->GetBufferPointer();
      const bool nifti = extension == ".nii.gz";
      for (std::size_t i = 0; i < lps.size(); ++i)
      {
        const bool isZ = (i / numVoxels) % 3 == 2;
        CPPUNIT_ASSERT_EQUAL(nifti && !isZ ? -lps[i] : lps[i], raw[i]);
      }
      std::string marker;
      itk::ExposeMetaData<std::string>(reader->GetOutput()->GetMetaDataDictionary(), "ITK_FileNotes", marker);
      if (nifti)
        CPPUNIT_ASSERT_EQUAL(std::string("MITK_peak_convention=RAS"), marker);

      // MITK's marker converts back whatever the MRtrix option says.
      CPPUNIT_ASSERT(Values(LoadPeaks(path, false)) == lps);
      CPPUNIT_ASSERT(Values(LoadPeaks(path, true)) == lps);
    }
  }

  /** Earlier MITK versions saved SH coefficients as one vector per voxel of a 3D image. */
  static std::vector<std::pair<std::string, unsigned int>> VectorLayoutShFiles()
  {
    return {{"qBallImage_SHCoeffs.nrrd", 4}, {"StreamlineTractography/qball_image_sh.nii.gz", 6}};
  }

  // MITK's image reader stays the default for SH coefficients in vector layout, as the component
  // count alone does not tell SH coefficients from other vectors; the SH reader is offered and
  // reads them unchanged. Saved by the SH writer and read back, they give the same coefficients,
  // so the same glyphs, and the same peaks.
  void ShReader_VectorLayout_OpensBothWays()
  {
    for (const auto &[file, order] : VectorLayoutShFiles())
    {
      const std::string path = mitk::GetDiffusionTestDataFilePath(file);
      mitk::Image::Pointer plain = dynamic_cast<mitk::Image *>(mitk::IOUtil::Load(path).at(0).GetPointer());
      CPPUNIT_ASSERT_MESSAGE(file, plain.IsNotNull() && dynamic_cast<mitk::ShImage *>(plain.GetPointer()) == nullptr);
      const std::size_t count = (order + 1) * (order + 2) / 2;
      CPPUNIT_ASSERT_EQUAL(count, plain->GetPixelType().GetNumberOfComponents());
      CPPUNIT_ASSERT_MESSAGE(file, ReaderDescriptions(path).count("SH Image") == 1);

      auto shImage = LoadSh(path);
      CPPUNIT_ASSERT_EQUAL(order, shImage->ShOrder());
      CPPUNIT_ASSERT_EQUAL(mitk::ShImage::FSL, shImage->GetShConvention());
      AssertEqualGeometry(plain, shImage);
      const std::size_t values = count * plain->GetDimension(0) * plain->GetDimension(1) * plain->GetDimension(2);
      mitk::ImageReadAccessor plainAccess(plain.GetPointer());
      mitk::ImageReadAccessor shAccess(shImage.GetPointer());
      CPPUNIT_ASSERT(std::memcmp(plainAccess.GetData(), shAccess.GetData(), values * sizeof(float)) == 0);

      const std::string savedPath = m_TempDir + "/vector_layout_sh" + std::to_string(order) + ".nii.gz";
      mitk::IOUtil::Save(shImage, savedPath);
      auto saved = LoadSh(savedPath);
      CPPUNIT_ASSERT_EQUAL(mitk::ShImage::FSL, saved->GetShConvention());
      AssertEqualPixels(shImage, saved);

      // Peaks in a block of 10 x 10 voxels of the middle slice, to keep the Debug test short.
      auto mask = mitk::Image::New();
      mask->Initialize(mitk::MakeScalarPixelType<unsigned char>(), *shImage->GetGeometry());
      {
        mitk::ImageWriteAccessor maskAccess(mask);
        auto *maskValues = static_cast<unsigned char *>(maskAccess.GetData());
        const unsigned int nx = mask->GetDimension(0), ny = mask->GetDimension(1);
        std::fill(maskValues, maskValues + static_cast<std::size_t>(nx) * ny * mask->GetDimension(2), 0);
        for (unsigned int y = 27; y < 37; ++y)
          for (unsigned int x = 27; x < 37; ++x)
            maskValues[x + nx * (y + static_cast<std::size_t>(ny))] = 1;
      }
      const auto peaks = mitk::ExtractPeaks(shImage, mitk::PeakParameters(), mask);
      AssertEqualPixels(peaks.peaks, mitk::ExtractPeaks(saved, mitk::PeakParameters(), mask).peaks);
      mitk::ImageReadAccessor numberAccess(peaks.numberOfPeaks);
      const auto *numbers = static_cast<const unsigned char *>(numberAccess.GetData());
      CPPUNIT_ASSERT_MESSAGE(file, std::accumulate(numbers, numbers + values / count, 0u) >= 100u);
    }
  }

  /** Mean Pearson correlation of the ODFs of an SH image with those of an ODF image, over voxels with an ODF. */
  static double MeanOdfCorrelation(mitk::ShImage *shImage, const std::string &odfFile)
  {
    auto shOdfs = mitk::convert::GetItkOdfFromShImage(shImage);
    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath(odfFile));
    auto odfs = mitk::convert::GetItkOdfFromOdfImage(odfImage);
    CPPUNIT_ASSERT(shOdfs->GetLargestPossibleRegion() == odfs->GetLargestPossibleRegion());
    double sum = 0;
    unsigned int voxels = 0;
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> itSh(shOdfs, shOdfs->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> it(odfs, odfs->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it, ++itSh)
    {
      const unsigned int n = ODF_SAMPLING_SIZE;
      double meanA = 0, meanB = 0;
      for (unsigned int i = 0; i < n; ++i)
      {
        meanA += itSh.Get()[i] / n;
        meanB += it.Get()[i] / n;
      }
      double ab = 0, aa = 0, bb = 0;
      for (unsigned int i = 0; i < n; ++i)
      {
        ab += (itSh.Get()[i] - meanA) * (it.Get()[i] - meanB);
        aa += (itSh.Get()[i] - meanA) * (itSh.Get()[i] - meanA);
        bb += (it.Get()[i] - meanB) * (it.Get()[i] - meanB);
      }
      if (aa > 0 && bb > 0)
      {
        sum += ab / std::sqrt(aa * bb);
        ++voxels;
      }
    }
    CPPUNIT_ASSERT(voxels > 0);
    return sum / voxels;
  }

  // The vector-layout files hold coefficients in MITK's FSL definition (the basis of MITK's
  // earlier Q-ball filters): read so, their ODFs follow the Q-ball ODFs saved from the same data;
  // read in MITK's own basis, they do not.
  void ShReader_VectorLayout_MitkFslBasis()
  {
    const std::string odfFiles[] = {"qBallImage.qbi", "StreamlineTractography/qball_image.qbi"};
    const auto files = VectorLayoutShFiles();
    for (unsigned int f = 0; f < files.size(); ++f)
    {
      auto shImage = LoadSh(mitk::GetDiffusionTestDataFilePath(files[f].first));
      const double fsl = MeanOdfCorrelation(shImage, odfFiles[f]);
      shImage->SetShConvention(mitk::ShImage::MRTRIX);
      const double mitkBasis = MeanOdfCorrelation(shImage, odfFiles[f]);
      CPPUNIT_ASSERT_MESSAGE(files[f].first, fsl > 0.95);
      CPPUNIT_ASSERT_MESSAGE(files[f].first, mitkBasis < 0.5);
    }
  }

  // A DWI in NRRD has one vector of volumes per voxel too; with 15 volumes it is no SH image.
  void MimeTypes_VectorLayoutDwi_NoShReader()
  {
    using VectorImageType = itk::VectorImage<float, 3>;
    for (const bool dwi : {false, true})
    {
      auto image = VectorImageType::New();
      VectorImageType::RegionType region;
      region.SetSize({{2, 2, 2}});
      image->SetRegions(region);
      image->SetVectorLength(15);
      image->Allocate();
      VectorImageType::PixelType pixel(15);
      pixel.Fill(1);
      image->FillBuffer(pixel);
      if (dwi)
        itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "modality", "DWMRI");
      const std::string path = m_TempDir + (dwi ? "/dwi.nrrd" : "/vectors.nrrd");
      auto writer = itk::ImageFileWriter<VectorImageType>::New();
      writer->SetInput(image);
      writer->SetFileName(path);
      writer->Update();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(path, dwi ? std::size_t(0) : std::size_t(1), ReaderDescriptions(path).count("SH Image"));
    }
  }

  // A 4D NIfTI file whose volume count fits SH and peak images is opened by MITK's image reader
  // unless the SH or peak reader is chosen; this rests on the mime-type rank and MITK's reader order.
  void MimeTypes_Sh4dNifti_ImageReaderIsDefault()
  {
    const std::string path = mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz");
    const auto descriptions = ReaderDescriptions(path);
    CPPUNIT_ASSERT(descriptions.count("SH Image") == 1 && descriptions.count("Peak Image") == 1);
    const std::string defaultReader = mitk::FileReaderSelector(path).GetDefault().GetDescription();
    CPPUNIT_ASSERT_MESSAGE(defaultReader, defaultReader != "SH Image" && defaultReader != "Peak Image");

    auto data = mitk::IOUtil::Load(path).at(0);
    CPPUNIT_ASSERT(dynamic_cast<mitk::Image *>(data.GetPointer()) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::ShImage *>(data.GetPointer()) == nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::PeakImage *>(data.GetPointer()) == nullptr);
  }

  // MRtrix (sh2peaks, 3.0.8 checked) writes peaks as scanner RAS vectors, whatever the image
  // orientation, and names itself in the NIfTI descrip; MITK reads them as world LPS vectors.
  void PeakReader_MrtrixKnownDirections_Lps()
  {
    using ImageType = mitk::PeakImage::ItkPeakImageType;
    auto image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize({{2, 1, 1, 6}});
    image->SetRegions(region);
    // An oblique image: the peaks are world vectors, the image rotation does not apply.
    const Matrix3 rotation = RotationZ30();
    ImageType::DirectionType direction;
    direction.SetIdentity();
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        direction[r][c] = rotation(r, c);
    image->SetDirection(direction);
    image->Allocate();
    // RAS peaks per voxel: (1, 0, 0) and (0, 0.6, 0.8); (0.6, 0.8, 0) and none.
    const float ras[2][6] = {{1, 0, 0, 0, 0.6f, 0.8f}, {0.6f, 0.8f, 0, 0, 0, 0}};
    for (unsigned int x = 0; x < 2; ++x)
      for (unsigned int j = 0; j < 6; ++j)
        image->SetPixel({{x, 0, 0, j}}, ras[x][j]);
    itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "ITK_FileNotes", "MRtrix version: 3.0.8");
    const std::string path = m_TempDir + "/mrtrix_peaks.nii.gz";
    auto writer = itk::ImageFileWriter<ImageType>::New();
    writer->SetInput(image);
    writer->SetFileName(path);
    writer->Update();

    auto peaks = LoadPeaks(path, std::nullopt)->GetItkPeakImage();
    for (unsigned int x = 0; x < 2; ++x)
      for (unsigned int j = 0; j < 6; ++j)
        CPPUNIT_ASSERT_EQUAL(j % 3 == 2 ? ras[x][j] : -ras[x][j], peaks->GetPixel({{x, 0, 0, j}}));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingIO)
