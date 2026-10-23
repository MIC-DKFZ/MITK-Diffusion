/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkFiberBundleIOHelper.h>
#include <mitkIOMimeTypes.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itksys/SystemTools.hxx>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
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
  MITK_TEST(TrkReader_ScalarsAndProperties_AreSkipped);
  MITK_TEST(TrkReader_InvalidFiles_Throw);
  MITK_TEST(TrkReader_NonRasVoxelOrders_FollowTrackVisDefinition);
  MITK_TEST(TrkReader_LegacyMitkLpsFile_ReadAsBefore);
  MITK_TEST(TrkReader_NoMatrixRecorded_UsesVoxelOrder);
  MITK_TEST(TrkWriter_RoundTrip_KeepsPoints);
  MITK_TEST(TrkWriter_Header_IsReproducibleAndMatchesData);
  MITK_TEST(Writers_WriteErrors_Throw);
  MITK_TEST(VtkReader_NoFiberPolyData_ThrowsMitkException);
  MITK_TEST(VtkReader_ColoursOfUnusedAndReorderedPoints_StayWithTheirPoints);
  MITK_TEST(DicomMimeType_NewDcmFile_FindsWriter);
  MITK_TEST(IOHelper_BrokenFileAndDottedNames_CoutRestored);
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

  using TrkHeader = mitk::FiberBundle::TrackVis_header;
  static_assert(sizeof(TrkHeader) == 1000, "TrackVis header layout");

  /** Two tracks in TrackVis voxmm coordinates (corner based) */
  static std::vector<Fiber> TrkFibers()
  {
    return {{{3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}, {9.0, 10.0, 11.5}}, {{1.0, 1.0, 1.0}, {2.0, 3.0, 4.0}}};
  }

  /** Version 2 header with the given voxel order and vox_to_ras, 2 mm voxels */
  static TrkHeader MakeTrkHeader(const std::string &voxelOrder, const std::array<std::array<float, 4>, 4> &voxToRas)
  {
    TrkHeader header;
    std::memset(&header, 0, sizeof(header));
    std::memcpy(header.id_string, "TRACK", 6);
    for (int i = 0; i < 3; ++i)
    {
      header.dim[i] = 20;
      header.voxel_size[i] = 2.0f;
    }
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        header.vox_to_ras[i][j] = voxToRas[i][j];
    std::memcpy(header.voxel_order, voxelOrder.c_str(), 3);
    header.version = 2;
    header.hdr_size = 1000;
    return header;
  }

  static std::array<std::array<float, 4>, 4> RasAffine()
  {
    return {{{2, 0, 0, -10}, {0, 2, 0, -20}, {0, 0, 2, 5}, {0, 0, 0, 1}}};
  }

  static void WriteTrk(const std::string &path, TrkHeader header, int numScalars = 0, int numProperties = 0)
  {
    header.n_scalars = static_cast<short>(numScalars);
    header.n_properties = static_cast<short>(numProperties);
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    for (const auto &fiber : TrkFibers())
    {
      const std::int32_t numPoints = static_cast<std::int32_t>(fiber.size());
      out.write(reinterpret_cast<const char *>(&numPoints), 4);
      for (const auto &p : fiber)
      {
        for (double c : p)
        {
          const float value = static_cast<float>(c);
          out.write(reinterpret_cast<const char *>(&value), 4);
        }
        const float scalar = 42.0f;
        for (int s = 0; s < numScalars; ++s)
          out.write(reinterpret_cast<const char *>(&scalar), 4);
      }
      const float property = -7.0f;
      for (int s = 0; s < numProperties; ++s)
        out.write(reinterpret_cast<const char *>(&property), 4);
    }
  }

  /** The bundle holds TrkFibers() mapped as the TrackVis format defines it,
   *  RAS = vox_to_ras * (voxmm / voxel_size - 0.5), converted to LPS. */
  static void AssertTrkFibersSpecMap(mitk::FiberBundle *fib, const TrkHeader &header, const std::string &message)
  {
    const auto fibers = TrkFibers();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<unsigned int>(fibers.size()), fib->GetNumFibers());
    for (unsigned int i = 0; i < fibers.size(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<vtkIdType>(fibers[i].size()), cell->GetNumberOfPoints());
      for (unsigned int j = 0; j < fibers[i].size(); ++j)
      {
        double ras[3];
        for (int r = 0; r < 3; ++r)
        {
          ras[r] = header.vox_to_ras[r][3];
          for (int c = 0; c < 3; ++c)
            ras[r] += header.vox_to_ras[r][c] * (fibers[i][j][c] / header.voxel_size[c] - 0.5);
        }
        const double *p = cell->GetPoints()->GetPoint(j);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -ras[0], p[0], 1e-4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -ras[1], p[1], 1e-4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, ras[2], p[2], 1e-4);
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
        std::string firstLine;
        std::getline(std::ifstream(path), firstLine);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Legacy format readable by VTK 8", std::string("# vtk DataFile Version 4.2"), firstLine);
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

  void TrkReader_ScalarsAndProperties_AreSkipped()
  {
    const auto header = MakeTrkHeader("RAS", RasAffine());
    const std::string path = m_TempDir + "/fibers.trk";
    WriteTrk(path, header);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "coordinates only");
    WriteTrk(path, header, 2, 1);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "2 scalars, 1 property");
  }

  void TrkReader_InvalidFiles_Throw()
  {
    const std::string valid = m_TempDir + "/valid.trk";
    WriteTrk(valid, MakeTrkHeader("RAS", RasAffine()));
    const auto validSize = std::filesystem::file_size(valid);

    const std::string path = m_TempDir + "/invalid.trk";
    std::filesystem::copy_file(valid, path);
    std::filesystem::resize_file(path, 500);
    CPPUNIT_ASSERT_THROW_MESSAGE("shorter than the header", mitk::IOUtil::Load(path), std::exception);

    std::filesystem::remove(path);
    std::filesystem::copy_file(valid, path);
    std::filesystem::resize_file(path, validSize - 6);
    CPPUNIT_ASSERT_THROW_MESSAGE("truncated in a point", mitk::IOUtil::Load(path), std::exception);

    auto header = MakeTrkHeader("RAS", RasAffine());
    std::memcpy(header.id_string, "TRICK", 6);
    WriteTrk(path, header);
    CPPUNIT_ASSERT_THROW_MESSAGE("wrong id_string", mitk::IOUtil::Load(path), std::exception);

    header = MakeTrkHeader("RAS", RasAffine());
    header.hdr_size = static_cast<int>(0xE8030000u); // 1000 with swapped bytes
    WriteTrk(path, header);
    CPPUNIT_ASSERT_THROW_MESSAGE("big-endian file", mitk::IOUtil::Load(path), std::exception);

    {
      std::ofstream out(path, std::ios::binary);
      const auto valid_header = MakeTrkHeader("RAS", RasAffine());
      out.write(reinterpret_cast<const char *>(&valid_header), sizeof(valid_header));
      const std::int32_t numPoints = 0;
      out.write(reinterpret_cast<const char *>(&numPoints), 4);
    }
    CPPUNIT_ASSERT_THROW_MESSAGE("track without points", mitk::IOUtil::Load(path), std::exception);
  }

  void TrkReader_NonRasVoxelOrders_FollowTrackVisDefinition()
  {
    const std::string path = m_TempDir + "/fibers.trk";

    // FSL/MNI orientation, e.g. dipy saving against an LAS image
    auto header = MakeTrkHeader("LAS", {{{-2, 0, 0, 30}, {0, 2, 0, -20}, {0, 0, 2, 5}, {0, 0, 0, 1}}});
    header.n_count = 2;
    WriteTrk(path, header);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "LAS");

    header = MakeTrkHeader("LPS", {{{-2, 0, 0, 30}, {0, -2, 0, 40}, {0, 0, 2, 5}, {0, 0, 0, 1}}});
    header.n_count = 2;
    WriteTrk(path, header);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "LPS");
  }

  void TrkReader_LegacyMitkLpsFile_ReadAsBefore()
  {
    // As MITK wrote by default from 2019-12 to 2023-01: voxel order LPS, MITK's LPS
    // index-to-world matrix as vox_to_ras, n_count 0, default image_orientation_patient
    auto header = MakeTrkHeader("LPS", RasAffine());
    header.image_orientation_patient[0] = 1;
    header.image_orientation_patient[4] = 1;
    const std::string path = m_TempDir + "/legacy.trk";
    WriteTrk(path, header);
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(path);

    const auto fibers = TrkFibers();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(fibers.size()), fib->GetNumFibers());
    for (unsigned int i = 0; i < fibers.size(); ++i)
      for (unsigned int j = 0; j < fibers[i].size(); ++j)
      {
        const double *p = fib->GetFiberPolyData()->GetCell(i)->GetPoints()->GetPoint(j);
        for (int r = 0; r < 3; ++r)
        {
          double lps = header.vox_to_ras[r][3];
          for (int c = 0; c < 3; ++c)
            lps += header.vox_to_ras[r][c] * (fibers[i][j][c] / header.voxel_size[c] - 0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(lps, p[r], 1e-4);
        }
      }
  }

  void TrkReader_NoMatrixRecorded_UsesVoxelOrder()
  {
    // vox_to_ras[3][3] == 0 means "not recorded"; points are voxmm, flipped to LPS by voxel order
    auto header = MakeTrkHeader("RAS", {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}});
    const std::string path = m_TempDir + "/nomatrix.trk";
    WriteTrk(path, header);
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(path);

    const auto fibers = TrkFibers();
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(fibers.size()), fib->GetNumFibers());
    for (unsigned int i = 0; i < fibers.size(); ++i)
      for (unsigned int j = 0; j < fibers[i].size(); ++j)
      {
        const double *p = fib->GetFiberPolyData()->GetCell(i)->GetPoints()->GetPoint(j);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-(fibers[i][j][0] - 1.0), p[0], 1e-4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-(fibers[i][j][1] - 1.0), p[1], 1e-4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(fibers[i][j][2] - 1.0, p[2], 1e-4);
      }
  }

  void TrkWriter_RoundTrip_KeepsPoints()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string path = m_TempDir + "/roundtrip.trk";
    for (bool useLps : {false, true})
    {
      mitk::IFileWriter::Options options;
      options["Use LPS coordinate system (default is RAS)"] = useLps;
      options["Print header"] = false;
      mitk::IOUtil::Save(fib, path, options);
      CPPUNIT_ASSERT_MESSAGE(useLps ? "LPS" : "RAS", fib->Equals(mitk::IOUtil::Load<mitk::FiberBundle>(path), 0.001));
    }

    // A bundle carrying the header of a third-party LAS file is written in that convention
    const std::string lasPath = m_TempDir + "/las.trk";
    auto header = MakeTrkHeader("LAS", {{{-2, 0, 0, 30}, {0, 2, 0, -20}, {0, 0, 2, 5}, {0, 0, 0, 1}}});
    header.n_count = 2;
    WriteTrk(lasPath, header);
    auto las = mitk::IOUtil::Load<mitk::FiberBundle>(lasPath);
    mitk::IOUtil::Save(las, path);
    CPPUNIT_ASSERT_MESSAGE("LAS header", las->Equals(mitk::IOUtil::Load<mitk::FiberBundle>(path), 0.001));
  }

  static std::string FileContent(const std::string &path)
  {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }

  void TrkWriter_Header_IsReproducibleAndMatchesData()
  {
    mitk::IFileWriter::Options options;
    options["Print header"] = false;

    // Two bundles loaded separately give byte-identical files
    auto first = LoadFibers("FiberProcessing/original.fib");
    auto second = LoadFibers("FiberProcessing/original.fib");
    const std::string path1 = m_TempDir + "/first.trk";
    const std::string path2 = m_TempDir + "/second.trk";
    mitk::IOUtil::Save(first, path1, options);
    mitk::IOUtil::Save(second, path2, options);
    const auto content = FileContent(path1);
    CPPUNIT_ASSERT_MESSAGE("Same bundle, same file", content == FileContent(path2));

    TrkHeader header;
    std::memcpy(&header, content.data(), sizeof(header));
    CPPUNIT_ASSERT_EQUAL(static_cast<int>(first->GetNumFibers()), header.n_count);
    CPPUNIT_ASSERT_EQUAL(static_cast<short>(0), header.n_scalars);
    CPPUNIT_ASSERT_EQUAL(static_cast<short>(0), header.n_properties);
    CPPUNIT_ASSERT_MESSAGE("Reserved bytes are zero",
                           std::all_of(header.reserved, header.reserved + sizeof(header.reserved), [](char c) { return c == 0; }));

    // A header taken over from a file with scalars and properties describes the written data
    const std::string withScalars = m_TempDir + "/scalars.trk";
    auto scalarHeader = MakeTrkHeader("RAS", RasAffine());
    scalarHeader.n_count = 2;
    WriteTrk(withScalars, scalarHeader, 2, 1);
    auto loaded = mitk::IOUtil::Load<mitk::FiberBundle>(withScalars);
    mitk::IOUtil::Save(loaded, path1, options);
    CPPUNIT_ASSERT(loaded->Equals(mitk::IOUtil::Load<mitk::FiberBundle>(path1), 0.001));
  }

  void Writers_WriteErrors_Throw()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string base = mitk::IOMimeTypes::DEFAULT_BASE_NAME();
    mitk::IFileWriter::Options options;
    options["Print header"] = false;

    const std::string missingDir = m_TempDir + "/missing/";
    CPPUNIT_ASSERT_THROW_MESSAGE(".fib, missing directory", mitk::IOUtil::Save(fib, missingDir + "fibers.fib"), std::exception);
    CPPUNIT_ASSERT_THROW_MESSAGE(".trk, missing directory", mitk::IOUtil::Save(fib, missingDir + "fibers.trk", options), std::exception);
    CPPUNIT_ASSERT_THROW_MESSAGE(".dcm, missing directory",
                                 mitk::IOUtil::Save(fib, base + ".FiberBundle.dcm", missingDir + "fibers.dcm"), std::exception);

    // A device on which every write fails (disk full); opening it succeeds.
    if (std::filesystem::exists("/dev/full"))
    {
      CPPUNIT_ASSERT_THROW_MESSAGE(".fib, disk full",
                                   mitk::IOUtil::Save(fib, base + ".FiberBundle.vtk", "/dev/full", options, false), std::exception);
      CPPUNIT_ASSERT_THROW_MESSAGE(".trk, disk full",
                                   mitk::IOUtil::Save(fib, base + ".FiberBundle.trk", "/dev/full", options, false), std::exception);
    }
  }

  void VtkReader_NoFiberPolyData_ThrowsMitkException()
  {
    const std::string grid = m_TempDir + "/grid.fib";
    {
      std::ofstream out(grid);
      out << "# vtk DataFile Version 3.0\ngrid\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 1 float\n0 0 0\n"
             "CELLS 1 2\n1 0\nCELL_TYPES 1\n1\n";
    }
    CPPUNIT_ASSERT_THROW_MESSAGE("unstructured grid", mitk::IOUtil::Load(grid), mitk::Exception);

    const std::string surface = m_TempDir + "/surface.fib";
    {
      std::ofstream out(surface);
      out << "# vtk DataFile Version 3.0\nsurface\nASCII\nDATASET POLYDATA\nPOINTS 3 float\n0 0 0 1 0 0 0 1 0\n"
             "POLYGONS 1 4\n3 0 1 2\n";
    }
    CPPUNIT_ASSERT_THROW_MESSAGE("polygons", mitk::IOUtil::Load(surface), mitk::Exception);
  }

  void VtkReader_ColoursOfUnusedAndReorderedPoints_StayWithTheirPoints()
  {
    // Point 0 is unused; fiber 0 uses points 3 and 4, fiber 1 uses points 1 and 2. Every point
    // has its own RGB colour (3 components, as other tools write them) and each fiber a weight.
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < 5; ++i)
      points->InsertNextPoint(i, 2.0 * i, 0.0);
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    const vtkIdType fiber0[2] = {3, 4};
    const vtkIdType fiber1[2] = {1, 2};
    lines->InsertNextCell(2, fiber0);
    lines->InsertNextCell(2, fiber1);
    auto colours = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colours->SetName("FIBER_COLORS");
    colours->SetNumberOfComponents(3);
    for (int i = 0; i < 5; ++i)
    {
      const unsigned char rgb[3] = {static_cast<unsigned char>(10 * i), static_cast<unsigned char>(10 * i + 1),
                                    static_cast<unsigned char>(10 * i + 2)};
      colours->InsertNextTypedTuple(rgb);
    }
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    weights->SetName("FIBER_WEIGHTS");
    weights->InsertNextValue(0.25f);
    weights->InsertNextValue(0.75f);
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    polyData->GetPointData()->AddArray(colours);
    polyData->GetCellData()->AddArray(weights);

    const std::string path = m_TempDir + "/reordered.fib";
    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(polyData);
    writer->SetFileName(path.c_str());
    writer->Write();

    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(path);
    CPPUNIT_ASSERT_EQUAL(2u, fib->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(0.25f, fib->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(0.75f, fib->GetFiberWeight(1));
    auto fibColours = fib->GetFiberColors();
    for (unsigned int f = 0; f < 2; ++f)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(f);
      for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
      {
        const vtkIdType id = cell->GetPointId(j);
        // the point's original index follows from its x coordinate
        const int original = static_cast<int>(std::lround(fib->GetFiberPolyData()->GetPoint(id)[0]));
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_EQUAL(10.0 * original + c, fibColours->GetComponent(id, c));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("RGB colours get an opaque alpha", 255.0, fibColours->GetComponent(id, 3));
      }
    }
  }

  void DicomMimeType_NewDcmFile_FindsWriter()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string path = m_TempDir + "/new.dcm";
    mitk::IOUtil::Save(fib, path);
    CPPUNIT_ASSERT(fib->Equals(mitk::IOUtil::Load<mitk::FiberBundle>(path), 0.001));

    // A file too short for the DICOM preamble is no tractogram
    const std::string shortFile = m_TempDir + "/short.dcm";
    std::ofstream(shortFile) << "DICM";
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load<mitk::FiberBundle>(shortFile), std::exception);
  }

  void IOHelper_BrokenFileAndDottedNames_CoutRestored()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    mitk::IOUtil::Save(fib, m_TempDir + "/sub-01.CST.fib");
    mitk::IOUtil::Save(fib, m_TempDir + "/plain.fib");

    const auto files = mitk::FiberBundleIOHelper::get_file_list(m_TempDir);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Names with several dots are listed", static_cast<std::size_t>(2), files.size());

    const std::string broken = m_TempDir + "/broken.trk";
    std::ofstream(broken) << "no tractogram";
    std::streambuf *coutBuffer = std::cout.rdbuf();
    CPPUNIT_ASSERT_THROW(mitk::FiberBundleIOHelper::load_fibs({m_TempDir + "/plain.fib", broken}), std::exception);
    CPPUNIT_ASSERT_MESSAGE("std::cout is restored after an error", coutBuffer == std::cout.rdbuf());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleIO)
