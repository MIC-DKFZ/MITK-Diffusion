/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkIOMimeTypes.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itksys/SystemTools.hxx>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <vector>

/** Regression tests for the fiber bundle readers and writers. Files are written to a temporary
 *  directory, nothing is added to MITK-Data. */
class mitkFiberBundleIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleIOTestSuite);
  MITK_TEST(VtkWriter_AllOptionCombinations_LeaveBundleUnchanged);
  MITK_TEST(DicomWriter_ReferencedImageUids_TractogramGetsOwnUids);
  MITK_TEST(TckReader_AllFloatTypes_ReadSamePoints);
  MITK_TEST(TckReader_ReadOnlyOrUnterminatedFile_Loads);
  MITK_TEST(TckReader_HeaderWithoutEnd_Throws);
  CPPUNIT_TEST_SUITE_END();

  std::string m_TempDir;

  mitk::FiberBundle::Pointer LoadFibers(const std::string &relativePath)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/" + relativePath));
  }

  /** original.fib with a different weight per fiber */
  mitk::FiberBundle::Pointer LoadWeightedFibers()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      fib->SetFiberWeight(i, 0.5f + static_cast<float>(i) / fib->GetNumFibers());
    return fib;
  }

  static std::set<std::string> ArrayNames(vtkPolyData *polyData)
  {
    std::set<std::string> names;
    for (int i = 0; i < polyData->GetPointData()->GetNumberOfArrays(); ++i)
      names.insert(std::string("point ") + polyData->GetPointData()->GetArrayName(i));
    for (int i = 0; i < polyData->GetCellData()->GetNumberOfArrays(); ++i)
      names.insert(std::string("cell ") + polyData->GetCellData()->GetArrayName(i));
    return names;
  }

  using Fiber = std::vector<std::array<double, 3>>;

  /** Two fibers in RAS millimetres, as MRtrix stores them */
  static std::vector<Fiber> TckFibers()
  {
    return {{{1.5, 2.0, 3.0}, {4.0, 5.25, 6.0}, {7.0, 8.0, 9.5}}, {{-1.0, -2.0, -3.0}, {-4.0, -5.0, 10.0}}};
  }

  template <typename T>
  static void WriteValue(std::ostream &out, double value, bool bigEndian)
  {
    std::array<char, sizeof(T)> bytes;
    const T typed = static_cast<T>(value);
    std::memcpy(bytes.data(), &typed, sizeof(T));
    if (bigEndian != (std::endian::native == std::endian::big))
      std::reverse(bytes.begin(), bytes.end());
    out.write(bytes.data(), bytes.size());
  }

  /** Writes an MRtrix .tck file: header, data at offset 1024, NaN after each track, Inf at the end. */
  template <typename T>
  static void WriteTck(const std::string &path,
                       const std::string &datatype,
                       const std::string &extraHeader = "",
                       bool terminator = true)
  {
    const bool bigEndian = datatype.substr(datatype.size() - 2) == "BE";
    const auto fibers = TckFibers();
    std::string header = "mrtrix tracks\n" + extraHeader + "datatype: " + datatype +
                         "\ncount: " + std::to_string(fibers.size()) + "\nfile: . 1024\nEND\n";
    header.resize(1024, '\0');
    std::ofstream out(path, std::ios::binary);
    out.write(header.data(), header.size());
    for (const auto &fiber : fibers)
    {
      for (const auto &p : fiber)
        for (double c : p)
          WriteValue<T>(out, c, bigEndian);
      for (int c = 0; c < 3; ++c)
        WriteValue<T>(out, std::numeric_limits<double>::quiet_NaN(), bigEndian);
    }
    if (terminator)
      for (int c = 0; c < 3; ++c)
        WriteValue<T>(out, std::numeric_limits<double>::infinity(), bigEndian);
  }

  /** The bundle holds TckFibers() in LPS (x and y negated). */
  static void AssertTckFibers(mitk::FiberBundle *fib, const std::string &message)
  {
    const auto fibers = TckFibers();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<unsigned int>(fibers.size()), fib->GetNumFibers());
    for (unsigned int i = 0; i < fibers.size(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<vtkIdType>(fibers[i].size()), cell->GetNumberOfPoints());
      for (unsigned int j = 0; j < fibers[i].size(); ++j)
      {
        const double *p = cell->GetPoints()->GetPoint(j);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -fibers[i][j][0], p[0], 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -fibers[i][j][1], p[1], 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, fibers[i][j][2], p[2], 1e-6);
      }
    }
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("FiberBundleIOTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void VtkWriter_AllOptionCombinations_LeaveBundleUnchanged()
  {
    auto fib = LoadWeightedFibers();
    const auto arraysBefore = ArrayNames(fib->GetFiberPolyData());
    const std::string path = m_TempDir + "/fibers.fib";

    for (bool saveWeights : {true, false})
    {
      for (bool saveColors : {true, false})
      {
        mitk::IFileWriter::Options options;
        options["Save fiber weights"] = saveWeights;
        options["Save color information"] = saveColors;
        mitk::IOUtil::Save(fib, path, options);

        CPPUNIT_ASSERT_MESSAGE("Saving does not change the arrays of the bundle",
                               arraysBefore == ArrayNames(fib->GetFiberPolyData()));

        auto reloaded = mitk::IOUtil::Load<mitk::FiberBundle>(path);
        CPPUNIT_ASSERT(reloaded->Equals(fib));
        for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
        {
          const float expected = saveWeights ? fib->GetFiberWeight(i) : 1.0f;
          CPPUNIT_ASSERT_EQUAL(expected, reloaded->GetFiberWeight(i));
        }
      }
    }
  }

  static std::string StringProperty(const mitk::FiberBundle *fib, const std::string &name)
  {
    std::string value;
    fib->GetPropertyList()->GetStringProperty(name.c_str(), value);
    return value;
  }

  void DicomWriter_ReferencedImageUids_TractogramGetsOwnUids()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string mrImageStorage = "1.2.840.10008.5.1.4.1.1.4";
    fib->GetPropertyList()->SetStringProperty("DICOM.patient_id", "TEST0001");
    fib->GetPropertyList()->SetStringProperty("DICOM.patient_name", "Test^Tractogram");
    fib->GetPropertyList()->SetStringProperty("DICOM.study_instance_uid", "1.2.826.0.1.3680043.2.1125.1");
    fib->GetPropertyList()->SetStringProperty("DICOM.series_instance_uid", "1.2.826.0.1.3680043.2.1125.2");
    fib->GetPropertyList()->SetStringProperty("DICOM.sop_instance_uid", "1.2.826.0.1.3680043.2.1125.3");
    fib->GetPropertyList()->SetStringProperty("DICOM.sop_class_uid", mrImageStorage.c_str());
    fib->GetPropertyList()->SetStringProperty("DICOM.frame_of_reference_uid", "1.2.826.0.1.3680043.2.1125.4");

    // The mime type helper class of the module is not exported; this is its DICOM name.
    const std::string mimeType = mitk::IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.dcm";
    const std::string path1 = m_TempDir + "/tractogram1.dcm";
    mitk::IOUtil::Save(fib, mimeType, path1);
    auto first = mitk::IOUtil::Load<mitk::FiberBundle>(path1);
    CPPUNIT_ASSERT(first->Equals(fib, 0.001));

    // The properties name the referenced image, as before saving
    for (const std::string name : {"DICOM.patient_id", "DICOM.study_instance_uid", "DICOM.series_instance_uid",
                                   "DICOM.sop_instance_uid", "DICOM.sop_class_uid", "DICOM.frame_of_reference_uid"})
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, StringProperty(fib, name), StringProperty(first, name));

    // The tractogram's own instance is neither the image nor its SOP class
    const auto ownInstance = StringProperty(first, "DICOM.tractogram.sop_instance_uid");
    const auto ownSeries = StringProperty(first, "DICOM.tractogram.series_instance_uid");
    CPPUNIT_ASSERT(!ownInstance.empty() && ownInstance != "0");
    CPPUNIT_ASSERT(ownInstance != StringProperty(fib, "DICOM.sop_instance_uid"));
    CPPUNIT_ASSERT(ownInstance != "1.2.840.10008.5.1.4.1.1.66.6");
    CPPUNIT_ASSERT(ownSeries != StringProperty(fib, "DICOM.series_instance_uid"));

    // Saving the loaded bundle again refers to the same image under a new instance UID
    const std::string path2 = m_TempDir + "/tractogram2.dcm";
    mitk::IOUtil::Save(first, mimeType, path2);
    auto second = mitk::IOUtil::Load<mitk::FiberBundle>(path2);
    CPPUNIT_ASSERT_EQUAL(StringProperty(fib, "DICOM.sop_instance_uid"), StringProperty(second, "DICOM.sop_instance_uid"));
    CPPUNIT_ASSERT(StringProperty(second, "DICOM.tractogram.sop_instance_uid") != ownInstance);
  }

  void TckReader_AllFloatTypes_ReadSamePoints()
  {
    const std::string path = m_TempDir + "/fibers.tck";
    WriteTck<float>(path, "Float32LE");
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "Float32LE");
    WriteTck<float>(path, "Float32BE");
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "Float32BE");
    WriteTck<double>(path, "Float64LE");
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "Float64LE");
    WriteTck<double>(path, "Float64BE");
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "Float64BE");

    // "END" inside a header value does not end the header
    WriteTck<float>(path, "Float32LE", "command_history: tckgen /data/LEGEND/BACKEND.mif\n");
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "END in a header value");
  }

  void TckReader_ReadOnlyOrUnterminatedFile_Loads()
  {
    const std::string path = m_TempDir + "/readonly.tck";
    WriteTck<float>(path, "Float32LE");
    std::filesystem::permissions(path, std::filesystem::perms::owner_read | std::filesystem::perms::group_read |
                                         std::filesystem::perms::others_read);
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(path), "read-only file");

    const std::string unterminated = m_TempDir + "/unterminated.tck";
    WriteTck<float>(unterminated, "Float32LE", "", false);
    AssertTckFibers(mitk::IOUtil::Load<mitk::FiberBundle>(unterminated), "file without end marker");
  }

  void TckReader_HeaderWithoutEnd_Throws()
  {
    const std::string path = m_TempDir + "/truncated.tck";
    {
      std::ofstream out(path, std::ios::binary);
      out << "mrtrix tracks\ndatatype: Float32LE\nfile: . 1024\n";
    }
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load(path), std::exception);

    {
      std::ofstream out(path, std::ios::binary);
      out << "mrtrix tracks\ndatatype: Float32LE\nEND\n";
    }
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load(path), std::exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleIO)
