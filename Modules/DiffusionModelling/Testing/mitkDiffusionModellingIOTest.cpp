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
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkBaseDataSerializer.h>

#include <mitkITKImageImport.h>

#include <itkImageFileWriter.h>
#include <itkVectorImage.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <itksys/SystemTools.hxx>
#include <vnl/vnl_matrix_fixed.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>

class mitkDiffusionModellingIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingIOTestSuite);
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
    mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/tensorImage.dti"));
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

    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(tensorImage, missingDir + "tensor.dti"), std::exception);

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi"));
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(odfImage, missingDir + "odf.odf"), std::exception);

    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
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
    CheckClone(mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti")));
    CheckClone(mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi")));

    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    auto peakImage = mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &peakFunctor);
    CheckClone(peakImage);
    CPPUNIT_ASSERT(peakImage->Clone().IsNotNull());

    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
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
    auto loaded = mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &functor);
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

  // The readers and the SH writer chose their ITK IO by everything after the first dot.
  void ShAndPeakIO_DotsInFileName()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
    const std::string shPath = m_TempDir + "/sh.v2.shi";
    mitk::IOUtil::Save(shImage, shPath);
    AssertEqualPixels(shImage, mitk::IOUtil::Load<mitk::ShImage>(shPath, &shFunctor));

    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    auto peakImage = mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &peakFunctor);
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
    for (const std::string file : {"DiffusionImaging/mean_FA.nii.gz", "DiffusionImaging/tensorImage_FA.nrrd"})
    {
      const auto descriptions = ReaderDescriptions(GetTestDataFilePath(file));
      CPPUNIT_ASSERT_MESSAGE(file, descriptions.count("SH Image") == 0);
      CPPUNIT_ASSERT_MESSAGE(file, descriptions.count("Peak Image") == 0);
    }

    const auto shDescriptions = ReaderDescriptions(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"));
    CPPUNIT_ASSERT(shDescriptions.count("SH Image") == 1);
    CPPUNIT_ASSERT(shDescriptions.count("Peak Image") == 1);
  }

  // The SH and peak mime types share the header of the last file checked; a file replaced under
  // the same name is read again.
  void MimeTypes_ReplacedFile_CheckedAgain()
  {
    const std::string path = m_TempDir + "/replaced.nii.gz";
    std::filesystem::copy_file(GetTestDataFilePath("DiffusionImaging/mean_FA.nii.gz"), path);
    CPPUNIT_ASSERT(ReaderDescriptions(path).count("SH Image") == 0);
    std::filesystem::copy_file(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), path, std::filesystem::copy_options::overwrite_existing);
    const auto descriptions = ReaderDescriptions(path);
    CPPUNIT_ASSERT(descriptions.count("SH Image") == 1);
    CPPUNIT_ASSERT(descriptions.count("Peak Image") == 1);
  }

  // Writers are chosen through the mime type of the output path, which does not exist yet.
  void ShWriter_NewNiftiFile_ReadsBackAsSh()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
    const std::string path = m_TempDir + "/sh.nii.gz";
    mitk::IOUtil::Save(shImage, path);
    AssertEqualPixels(shImage, mitk::IOUtil::Load<mitk::ShImage>(path, &shFunctor));
  }

  // An image read as FSL/Dipy SH was saved without its convention and read back as MRtrix.
  void ShIO_ConventionSurvivesSaving()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
    for (auto convention : {mitk::ShImage::FSL, mitk::ShImage::MRTRIX})
    {
      shImage->SetShConvention(convention);
      for (const std::string extension : {".nii.gz", ".shi"})
      {
        const std::string path = m_TempDir + "/convention" + std::to_string(convention) + extension;
        mitk::IOUtil::Save(shImage, path);
        auto reloaded = mitk::IOUtil::Load<mitk::ShImage>(path, &shFunctor);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(path, convention, reloaded->GetShConvention());
        AssertEqualPixels(shImage, reloaded);
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
    const std::string mrtrixFile = GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_peak_image.nii.gz");
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

    // Peak files without the MRtrix description, e.g. all MITK-written ones, read unchanged.
    const std::string mitkFile = GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz");
    CPPUNIT_ASSERT(Values(LoadPeaks(mitkFile, false)) == Values(LoadPeaks(mitkFile, true)));
  }

  // The SH importer and exporter converted origins through float.
  void ShIO_OriginKeepsDoublePrecision()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
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

  // Pixel buffers and geometry of MITK-Data files as read before the readers were restructured.
  void Readers_RecordedContent()
  {
    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    mitk::PreferenceListReaderOptionsFunctor peakFunctor({"Peak Image"}, std::vector<std::string>());
    struct Recorded
    {
      mitk::Image::Pointer image;
      std::uint64_t hash;
      double origin[3];
      double spacing[3];
    };
    const Recorded recorded[] = {
      {mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/tensorImage.dti")).GetPointer(),
       910238346892541248ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/FWE/ReferenceTensorImage.dti")).GetPointer(),
       14634070166575176342ull, {-14.1639, 13.4743, -12.6538}, {2, 2, 2}},
      {mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi")).GetPointer(),
       17363730383617297091ull, {1, 1, 1}, {2, 2, 2}},
      {mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image.qbi")).GetPointer(),
       13068465665225140656ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor).GetPointer(),
       13341037711982211669ull, {0.615, 1.0005, 0.7045}, {1.23, 2.001, 1.409}},
      {mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image_sh.nii.gz"), &shFunctor).GetPointer(),
       8054354499704950156ull, {0, 0, 0}, {3, 3, 3}},
      {mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &peakFunctor).GetPointer(),
       1038741785062553984ull, {0.615, 1.0005, 0.7045}, {1.23, 2.001, 1.409}},
      {mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/FiberFit/csd_peak_image.nii.gz"), &peakFunctor).GetPointer(),
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
    for (const std::string file : {"DiffusionImaging/tensorImage.dti", "DiffusionImaging/FWE/ReferenceTensorImage.dti"})
    {
      auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath(file));
      const std::string path = m_TempDir + "/tensor.dti";
      mitk::IOUtil::Save(tensorImage, path);
      auto reloaded = mitk::IOUtil::Load<mitk::TensorImage>(path);
      AssertEqualPixels(tensorImage, reloaded);
      AssertEqualGeometry(tensorImage, reloaded);
    }

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi"));
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
    for (const std::string file : {"DiffusionImaging/peak_image_test.nii.gz", "DiffusionImaging/FiberFit/csd_peak_image.nii.gz"})
    {
      auto peakImage = LoadPeaks(GetTestDataFilePath(file), std::nullopt);
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
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    auto tensorReloaded = mitk::IOUtil::Load<mitk::TensorImage>(Serialize(tensorImage, "tensor"));
    AssertEqualPixels(tensorImage, tensorReloaded);
    AssertEqualGeometry(tensorImage, tensorReloaded);

    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi"));
    auto odfReloaded = mitk::IOUtil::Load<mitk::OdfImage>(Serialize(odfImage, "odf"));
    AssertEqualPixels(odfImage, odfReloaded);
    AssertEqualGeometry(odfImage, odfReloaded);

    mitk::PreferenceListReaderOptionsFunctor shFunctor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &shFunctor);
    auto shReloaded = mitk::IOUtil::Load<mitk::ShImage>(Serialize(shImage, "sh"), &shFunctor);
    AssertEqualPixels(shImage, shReloaded);
    AssertEqualGeometry(shImage, shReloaded);

    auto peakImage = LoadPeaks(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), std::nullopt);
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
    auto image = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
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
    auto image = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    CPPUNIT_ASSERT_EQUAL(mitk::ShImage::MRTRIX, image->GetShConvention());
    const auto mrtrix = Colours(image);

    image->SetShConvention(mitk::ShImage::FSL);
    const auto fsl = Colours(image);
    CPPUNIT_ASSERT(mrtrix != fsl);

    auto reference = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    reference->SetShConvention(mitk::ShImage::FSL);
    CPPUNIT_ASSERT(fsl == Colours(reference));
  }

  // The peak lines were built once and kept when the peaks changed later; a rebuilt polydata
  // keeps the custom colour.
  void PeakImage_PolydataFollowsPixelChanges()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    auto peaks = mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &functor);
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingIO)
