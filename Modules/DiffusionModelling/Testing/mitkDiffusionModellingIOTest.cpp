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
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkITKImageImport.h>

#include <itkVectorImage.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <itksys/SystemTools.hxx>
#include <vnl/vnl_matrix_fixed.h>

#include <cmath>
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
  MITK_TEST(ShWriter_NewNiftiFile_ReadsBackAsSh);
  MITK_TEST(ShIO_ConventionSurvivesSaving);
  MITK_TEST(PeakReader_MrtrixPeaks_ConvertedToLps);
  MITK_TEST(ShIO_OriginKeepsDoublePrecision);
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingIO)
