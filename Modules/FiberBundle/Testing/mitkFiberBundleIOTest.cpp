/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkFiberBundleIOHelper.h>
#include <mitkFiberBundleDicomReferenceTags.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkImage.h>
#include <mitkTemporoSpatialStringProperty.h>
#include <mitkIOMimeTypes.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkCustomMimeType.h>
#include <mitkBaseDataSerializer.h>
#include <usModule.h>
#include <usModuleContext.h>
#include <usModuleRegistry.h>

#include <itkObjectFactoryBase.h>
#include <itksys/SystemTools.hxx>
#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
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
#include <clocale>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
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
  MITK_TEST(TrkReader_VoxelOrderDisagreesWithMatrix_ReadAsNibabel);
  MITK_TEST(TrkReader_LegacyMitkRasFile_ReadAsBefore);
  MITK_TEST(TrkReader_NoMatrixRecorded_UsesVoxelOrder);
  MITK_TEST(TrkWriter_RoundTrip_KeepsPoints);
  MITK_TEST(TrkWriter_Header_IsReproducibleAndMatchesData);
  MITK_TEST(Writers_WriteErrors_Throw);
  MITK_TEST(VtkReader_NoFiberPolyData_ThrowsMitkException);
  MITK_TEST(VtkReader_TruncatedFile_Throws);
  MITK_TEST(VtkReader_ColoursOfUnusedAndReorderedPoints_StayWithTheirPoints);
  MITK_TEST(DicomMimeType_NewDcmFile_FindsWriter);
  MITK_TEST(IOHelper_BrokenFileAndDottedNames_CoutRestored);
  MITK_TEST(CopyDicomReferenceTags_BothNamingsAndMissingTags);
  MITK_TEST(TrkWriter_Header_DescribesReferenceGrid);
  MITK_TEST(TrkWriter_LpsOptionAndOldMitkFiles_WriteConsistentHeaders);
  MITK_TEST(DicomReader_NoAlgorithmIdentification_LoadsAndKeepsLocale);
  MITK_TEST(PlanarFigureComposite_RoundTripAndBrokenFile_KeepLocale);
  MITK_TEST(PlanarFigureComposite_UnknownOperation_Throws);
  MITK_TEST(MimeTypes_RegisteredOncePerName);
  MITK_TEST(VtkWriter_BinaryAsciiAndXml_KeepPointsWeightsAndColours);
  MITK_TEST(VtkReaderWriter_RasOption_RoundTrip);
  MITK_TEST(PlanarFigureComposite_AllOperations_RoundTrip);
  MITK_TEST(Serializers_BundleAndComposite_RoundTrip);
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
    // padded with spaces as MRtrix writes its header lines
    std::string header = "mrtrix tracks    \n" + extraHeader + "datatype: " + datatype +
                         "\ncount: " + std::to_string(fibers.size()) + "    \nfile: . 1024\nEND\n";
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

  /** RAS of a TrackVis point as nibabel and dipy read it (get_affine_trackvis_to_rasmm), for
   *  axis-aligned matrices: the voxel index voxmm / voxel_size - 0.5 is taken from the axes the voxel
   *  order names to those of vox_to_ras (component i of the matrix index is the voxel-order index
   *  of the same world axis, negated within dim when the directions differ), then mapped by vox_to_ras. */
  static std::array<double, 3> SpecRas(const TrkHeader &header, const float voxmm[3])
  {
    auto matrixAxis = [&header](int c)
    {
      int r = 0;
      for (int k = 1; k < 3; ++k)
        if (std::fabs(header.vox_to_ras[k][c]) > std::fabs(header.vox_to_ras[r][c]))
          r = k;
      return std::make_pair(r, header.vox_to_ras[r][c] < 0 ? -1 : 1);
    };
    auto orderAxis = [&header](int i)
    {
      const auto pos = std::string("RLAPSI").find(header.voxel_order[i]);
      return std::make_pair(static_cast<int>(pos / 2), pos % 2 == 0 ? 1 : -1);
    };
    double index[3];
    for (int i = 0; i < 3; ++i)
      index[i] = voxmm[i] / header.voxel_size[i] - 0.5;
    double matrixIndex[3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (matrixAxis(j).first == orderAxis(i).first)
        {
          const bool flip = matrixAxis(j).second != orderAxis(i).second;
          matrixIndex[i] = flip ? header.dim[i] - 1 - index[j] : index[j];
        }
    std::array<double, 3> ras;
    for (int r = 0; r < 3; ++r)
    {
      ras[r] = header.vox_to_ras[r][3];
      for (int c = 0; c < 3; ++c)
        ras[r] += header.vox_to_ras[r][c] * matrixIndex[c];
    }
    return ras;
  }

  /** The bundle holds TrkFibers() mapped as nibabel and dipy read the header (SpecRas), converted to LPS. */
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
        const float voxmm[3] = {static_cast<float>(fibers[i][j][0]), static_cast<float>(fibers[i][j][1]), static_cast<float>(fibers[i][j][2])};
        const auto ras = SpecRas(header, voxmm);
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

  void DicomReader_NoAlgorithmIdentification_LoadsAndKeepsLocale()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string mimeType = mitk::IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.dcm";
    const std::string path = m_TempDir + "/tractogram.dcm";
    mitk::IOUtil::Save(fib, mimeType, path);

    // The Tracking Algorithm Identification Sequence is type 1, but DCMTK reads a file without it.
    DcmFileFormat file;
    CPPUNIT_ASSERT(file.loadFile(path.c_str()).good());
    DcmItem *trackSet = nullptr;
    CPPUNIT_ASSERT(file.getDataset()->findAndGetSequenceItem(DCM_TrackSetSequence, trackSet, 0).good());
    CPPUNIT_ASSERT(trackSet->findAndDeleteElement(DCM_TrackingAlgorithmIdentificationSequence).good());
    const std::string stripped = m_TempDir + "/no_algorithm.dcm";
    CPPUNIT_ASSERT(file.saveFile(stripped.c_str(), EXS_LittleEndianExplicit).good());

    // A locale other than "C" must survive the reader, whether it succeeds or not.
    const std::string before = std::setlocale(LC_ALL, nullptr);
    const char *testLocale = std::setlocale(LC_ALL, "C.UTF-8");
    const std::string expected = testLocale != nullptr ? testLocale : before;
    mitk::FiberBundle::Pointer loaded;
    try
    {
      loaded = mitk::IOUtil::Load<mitk::FiberBundle>(stripped);
    }
    catch (...)
    {
    }
    const std::string after = std::setlocale(LC_ALL, nullptr);
    std::setlocale(LC_ALL, before.c_str());

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Locale restored", expected, after);
    CPPUNIT_ASSERT_MESSAGE("File without algorithm identification loads", loaded.IsNotNull());
    CPPUNIT_ASSERT(loaded->Equals(fib, 0.001));
    CPPUNIT_ASSERT_EQUAL(std::string("0"), StringProperty(loaded, "DICOM.algo_family_code.value"));
  }

  void PlanarFigureComposite_RoundTripAndBrokenFile_KeepLocale()
  {
    auto composite = mitk::PlanarFigureComposite::New();
    composite->setOperationType(mitk::PlanarFigureComposite::NOT);
    const std::string path = m_TempDir + "/composite.pfc";
    const std::string broken = m_TempDir + "/broken.pfc";
    std::ofstream(broken) << "no xml";

    const std::string before = std::setlocale(LC_ALL, nullptr);
    const char *testLocale = std::setlocale(LC_ALL, "C.UTF-8");
    const std::string expected = testLocale != nullptr ? testLocale : before;
    mitk::IOUtil::Save(composite, path);
    const std::string afterWriting = std::setlocale(LC_ALL, nullptr);
    auto loaded = mitk::IOUtil::Load<mitk::PlanarFigureComposite>(path);
    bool brokenThrew = false;
    try
    {
      mitk::IOUtil::Load(broken);
    }
    catch (...)
    {
      brokenThrew = true;
    }
    const std::string afterBroken = std::setlocale(LC_ALL, nullptr);
    std::setlocale(LC_ALL, before.c_str());

    CPPUNIT_ASSERT_EQUAL(mitk::PlanarFigureComposite::NOT, loaded->getOperationType());
    CPPUNIT_ASSERT(brokenThrew);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Locale after writing", expected, afterWriting);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Locale after a failed read", expected, afterBroken);
  }

  void PlanarFigureComposite_UnknownOperation_Throws()
  {
    const std::string path = m_TempDir + "/unknown.pfc";
    std::ofstream(path) << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<comptype>7</comptype>\n";
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Load(path), mitk::Exception);

    std::ofstream(path) << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<comptype>1</comptype>\n";
    CPPUNIT_ASSERT_EQUAL(mitk::PlanarFigureComposite::OR,
                         mitk::IOUtil::Load<mitk::PlanarFigureComposite>(path)->getOperationType());
  }

  void MimeTypes_RegisteredOncePerName()
  {
    // Loading a file loads the module and its activator. The test driver has no module context
    // of its own; every context sees all services.
    LoadFibers("FiberProcessing/original.fib");
    auto module = us::ModuleRegistry::GetModule("MitkFiberBundle");
    CPPUNIT_ASSERT(module != nullptr);
    auto context = module->GetModuleContext();
    std::map<std::string, int> registrations;
    for (const auto &reference : context->GetServiceReferences<mitk::CustomMimeType>())
    {
      // Only this module's registrations are dereferenced: MITK's MitkDICOMRTIO leaves a
      // registered mime type whose object is no longer valid (seen while writing this test).
      if (reference.GetModule() == nullptr || reference.GetModule()->GetName() != "MitkFiberBundle")
        continue;
      const auto *mimeType = context->GetService(reference);
      if (mimeType != nullptr)
        ++registrations[mimeType->GetName()];
      context->UngetService(reference);
    }
    const std::string base = mitk::IOMimeTypes::DEFAULT_BASE_NAME();
    for (const std::string suffix : {".pfc", ".FiberBundle.vtk", ".FiberBundle.trk", ".FiberBundle.tck", ".FiberBundle.dcm"})
      CPPUNIT_ASSERT_EQUAL_MESSAGE(base + suffix, 1, registrations[base + suffix]);
  }

  /** Same fibers, weights and colours; eps for coordinates and weights (ASCII files round them) */
  static void AssertSameBundle(mitk::FiberBundle *expected, mitk::FiberBundle *actual, double eps, const std::string &message)
  {
    CPPUNIT_ASSERT_MESSAGE(message, expected->Equals(actual, eps));
    for (unsigned int i = 0; i < expected->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, expected->GetFiberWeight(i), actual->GetFiberWeight(i), eps);
    auto expectedColours = expected->GetFiberColors();
    auto actualColours = actual->GetFiberColors();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, expectedColours->GetNumberOfValues(), actualColours->GetNumberOfValues());
    // colours belong to points; compare them along the fibers
    for (unsigned int i = 0; i < expected->GetNumFibers(); ++i)
    {
      vtkCell *e = expected->GetFiberPolyData()->GetCell(i);
      vtkCell *a = actual->GetFiberPolyData()->GetCell(i);
      for (vtkIdType j = 0; j < e->GetNumberOfPoints(); ++j)
        for (int c = 0; c < 4; ++c)
          CPPUNIT_ASSERT_EQUAL_MESSAGE(message, expectedColours->GetComponent(e->GetPointId(j), c),
                                       actualColours->GetComponent(a->GetPointId(j), c));
    }
  }

  void VtkWriter_BinaryAsciiAndXml_KeepPointsWeightsAndColours()
  {
    auto fib = LoadWeightedFibers();
    fib->ColorFibersByLength(false, false, mitk::LookupTable::JET);
    for (bool xml : {false, true})
    {
      for (bool binary : {true, false})
      {
        const std::string message = std::string(xml ? "XML" : "legacy") + (binary ? " binary" : " ASCII");
        const std::string path = m_TempDir + "/fibers_" + (xml ? "xml" : "legacy") + (binary ? "_binary" : "_ascii") + ".fib";
        mitk::IFileWriter::Options options;
        options["Save as xml file (vtp style)"] = xml;
        options["Save as binary file"] = binary;
        mitk::IOUtil::Save(fib, path, options);
        AssertSameBundle(fib, mitk::IOUtil::Load<mitk::FiberBundle>(path), binary ? 0.0 : 1e-4, message);
      }
    }
  }

  void VtkReaderWriter_RasOption_RoundTrip()
  {
    const std::string input = GetTestDataFilePath("DiffusionImaging/FiberProcessing/original.fib");
    auto lps = mitk::IOUtil::Load<mitk::FiberBundle>(input);
    mitk::IFileReader::Options rasOption;
    rasOption["Use RAS space"] = true;
    mitk::FiberBundle::Pointer ras = dynamic_cast<mitk::FiberBundle *>(mitk::IOUtil::Load(input, rasOption).at(0).GetPointer());
    CPPUNIT_ASSERT(ras.IsNotNull());
    CPPUNIT_ASSERT(ras->IsRAS());
    CPPUNIT_ASSERT(!lps->IsRAS());
    CPPUNIT_ASSERT_EQUAL(lps->GetNumFibers(), ras->GetNumFibers());
    for (vtkIdType i = 0; i < lps->GetFiberPolyData()->GetNumberOfPoints(); i += 97)
    {
      const double *p = lps->GetFiberPolyData()->GetPoint(i);
      const double *q = ras->GetFiberPolyData()->GetPoint(i);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-p[0], q[0], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-p[1], q[1], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(p[2], q[2], 1e-6);
    }

    // a RAS bundle is written back in RAS: the file holds the coordinates of the input file
    const std::string path = m_TempDir + "/ras.fib";
    mitk::IOUtil::Save(ras, path);
    CPPUNIT_ASSERT(lps->Equals(mitk::IOUtil::Load<mitk::FiberBundle>(path), 1e-5));
    mitk::FiberBundle::Pointer rasAgain = dynamic_cast<mitk::FiberBundle *>(mitk::IOUtil::Load(path, rasOption).at(0).GetPointer());
    CPPUNIT_ASSERT(ras->Equals(rasAgain, 1e-5));
  }

  void PlanarFigureComposite_AllOperations_RoundTrip()
  {
    for (auto operation : {mitk::PlanarFigureComposite::AND, mitk::PlanarFigureComposite::OR, mitk::PlanarFigureComposite::NOT})
    {
      auto composite = mitk::PlanarFigureComposite::New();
      composite->setOperationType(operation);
      const std::string path = m_TempDir + "/composite" + std::to_string(operation) + ".pfc";
      mitk::IOUtil::Save(composite, path);
      CPPUNIT_ASSERT_EQUAL(operation, mitk::IOUtil::Load<mitk::PlanarFigureComposite>(path)->getOperationType());
    }
  }

  /** The scene serializer registered under the given name, as mitk::SceneIO creates it */
  static mitk::BaseDataSerializer::Pointer CreateSerializer(const std::string &name)
  {
    for (const auto &object : itk::ObjectFactoryBase::CreateAllInstance(name.c_str()))
      if (auto *serializer = dynamic_cast<mitk::BaseDataSerializer *>(object.GetPointer()))
        return serializer;
    return nullptr;
  }

  void Serializers_BundleAndComposite_RoundTrip()
  {
    auto fib = LoadWeightedFibers();
    fib->ColorFibersByLength(false, false, mitk::LookupTable::JET);
    auto serializer = CreateSerializer("FiberBundleSerializer");
    CPPUNIT_ASSERT(serializer.IsNotNull());
    serializer->SetData(fib);
    serializer->SetFilenameHint("bundle");
    serializer->SetWorkingDirectory(m_TempDir);
    const std::string fibFile = serializer->Serialize();
    CPPUNIT_ASSERT(!fibFile.empty());
    AssertSameBundle(fib, mitk::IOUtil::Load<mitk::FiberBundle>(m_TempDir + "/" + fibFile), 0.0, "FiberBundleSerializer");

    auto composite = mitk::PlanarFigureComposite::New();
    composite->setOperationType(mitk::PlanarFigureComposite::OR);
    serializer = CreateSerializer("PlanarFigureCompositeSerializer");
    CPPUNIT_ASSERT(serializer.IsNotNull());
    serializer->SetData(composite);
    serializer->SetFilenameHint("composite");
    serializer->SetWorkingDirectory(m_TempDir);
    const std::string pfcFile = serializer->Serialize();
    CPPUNIT_ASSERT(!pfcFile.empty());
    CPPUNIT_ASSERT_EQUAL(mitk::PlanarFigureComposite::OR,
                         mitk::IOUtil::Load<mitk::PlanarFigureComposite>(m_TempDir + "/" + pfcFile)->getOperationType());

    // a serializer given other data writes nothing
    serializer = CreateSerializer("FiberBundleSerializer");
    serializer->SetData(composite);
    serializer->SetWorkingDirectory(m_TempDir);
    CPPUNIT_ASSERT(serializer->Serialize().empty());
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

    header = MakeTrkHeader("RAS", RasAffine());
    header.voxel_size[1] = 0;
    WriteTrk(path, header);
    CPPUNIT_ASSERT_THROW_MESSAGE("voxel size 0", mitk::IOUtil::Load(path), std::exception);

    header = MakeTrkHeader("RAS", {{{2, 0, 0, -10}, {0, 2, 0, -20}, {0, 0, 0, 5}, {0, 0, 0, 1}}});
    WriteTrk(path, header);
    CPPUNIT_ASSERT_THROW_MESSAGE("singular vox_to_ras", mitk::IOUtil::Load(path), std::exception);

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

  void TrkReader_VoxelOrderDisagreesWithMatrix_ReadAsNibabel()
  {
    const std::string path = m_TempDir + "/fibers.trk";

    // nibabel writes voxel order "RAS" when a script passes only the affine of an LAS image
    auto header = MakeTrkHeader("RAS", {{{-2, 0, 0, 30}, {0, 2, 0, -20}, {0, 0, 2, 5}, {0, 0, 0, 1}}});
    header.n_count = 2;
    WriteTrk(path, header);
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(path);
    AssertTrkFibersSpecMap(fib, header, "RAS order, LAS matrix");
    // voxmm (3, 4, 5) is index (1, 1.5, 2) of the RAS order, index (18, 1.5, 2) of the LAS matrix
    const double *p = fib->GetFiberPolyData()->GetPoint(0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, p[0], 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(17.0, p[1], 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0, p[2], 1e-4);

    header = MakeTrkHeader("LPS", RasAffine());
    header.n_count = 2;
    WriteTrk(path, header);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "LPS order, RAS matrix");

    header = MakeTrkHeader("ASR", RasAffine());
    header.n_count = 2;
    WriteTrk(path, header);
    AssertTrkFibersSpecMap(mitk::IOUtil::Load<mitk::FiberBundle>(path), header, "permuted order");
  }

  void TrkReader_LegacyMitkRasFile_ReadAsBefore()
  {
    // As MITK wrote by default from 2023 to 2026: voxel order RAS with MITK's LPS index-to-world
    // matrix (here of an LAS image, whose axes nibabel reads as RPS), n_count 0, default
    // image_orientation_patient; read with the matrix and x, y negated, without reorientation.
    auto header = MakeTrkHeader("RAS", {{{2, 0, 0, -30}, {0, -2, 0, 20}, {0, 0, 2, 5}, {0, 0, 0, 1}}});
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
          double ras = header.vox_to_ras[r][3];
          for (int c = 0; c < 3; ++c)
            ras += header.vox_to_ras[r][c] * (fibers[i][j][c] / header.voxel_size[c] - 0.5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(r < 2 ? -ras : ras, p[r], 1e-4);
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

  void VtkReader_TruncatedFile_Throws()
  {
    auto fib = LoadWeightedFibers();
    for (bool binary : {true, false})
    {
      const std::string full = m_TempDir + (binary ? "/full_binary.fib" : "/full_ascii.fib");
      mitk::IFileWriter::Options options;
      options["Save as binary file"] = binary;
      mitk::IOUtil::Save(fib, full, options);
      CPPUNIT_ASSERT_EQUAL(fib->GetNumFibers(), mitk::IOUtil::Load<mitk::FiberBundle>(full)->GetNumFibers());

      // ends in the points, the lines, the weights and the colours
      const auto size = std::filesystem::file_size(full);
      for (double fraction : {0.05, 0.6, 0.8, 0.95})
      {
        const std::string path = m_TempDir + "/truncated.fib";
        std::filesystem::copy_file(full, path, std::filesystem::copy_options::overwrite_existing);
        std::filesystem::resize_file(path, static_cast<std::uintmax_t>(static_cast<double>(size) * fraction));
        CPPUNIT_ASSERT_THROW_MESSAGE((binary ? "binary, " : "ASCII, ") + std::to_string(fraction), mitk::IOUtil::Load(path), std::exception);
      }
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

  void CopyDicomReferenceTags_BothNamingsAndMissingTags()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");

    // names of the diffusion DICOM import
    auto dwi = mitk::Image::New();
    dwi->GetPropertyList()->SetStringProperty("DICOM.patient_id", "TEST0001");
    dwi->GetPropertyList()->SetStringProperty("DICOM.series_instance_uid", "1.2.826.0.1.3680043.2.1125.2");
    auto missing = mitk::CopyDicomReferenceTags(dwi, fib);
    CPPUNIT_ASSERT_EQUAL(std::string("TEST0001"), StringProperty(fib, "DICOM.patient_id"));
    CPPUNIT_ASSERT_EQUAL(std::string("1.2.826.0.1.3680043.2.1125.2"), StringProperty(fib, "DICOM.series_instance_uid"));
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(5), missing.size());

    // MITK's DICOM tag properties
    auto anatomical = mitk::Image::New();
    anatomical->SetProperty("DICOM.0010.0020", mitk::TemporoSpatialStringProperty::New("TEST0002"));
    anatomical->SetProperty("DICOM.0008.0016", mitk::TemporoSpatialStringProperty::New("1.2.840.10008.5.1.4.1.1.4"));
    missing = mitk::CopyDicomReferenceTags(anatomical, fib);
    CPPUNIT_ASSERT_EQUAL(std::string("TEST0002"), StringProperty(fib, "DICOM.patient_id"));
    CPPUNIT_ASSERT_EQUAL(std::string("1.2.840.10008.5.1.4.1.1.4"), StringProperty(fib, "DICOM.sop_class_uid"));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("A tag the image lacks stays as it was", std::string("1.2.826.0.1.3680043.2.1125.2"),
                                 StringProperty(fib, "DICOM.series_instance_uid"));

    // nothing to copy
    missing = mitk::CopyDicomReferenceTags(mitk::Image::New(), fib);
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(7), missing.size());
    CPPUNIT_ASSERT_EQUAL(std::string("TEST0002"), StringProperty(fib, "DICOM.patient_id"));
  }

  /** Reads a .trk file as nibabel and dipy do (SpecRas); inside tells whether all points lie in the
   *  declared volume (dipy's bounding box check). */
  static std::vector<std::array<double, 3>> ReadTrkAsSpec(const std::string &path, TrkHeader &header, bool &inside)
  {
    std::ifstream in(path, std::ios::binary);
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    std::vector<std::array<double, 3>> ras;
    inside = true;
    std::int32_t numPoints = 0;
    while (in.read(reinterpret_cast<char *>(&numPoints), 4))
    {
      for (int j = 0; j < numPoints; ++j)
      {
        float voxmm[3];
        in.read(reinterpret_cast<char *>(voxmm), 12);
        const auto p = SpecRas(header, voxmm);
        for (int c = 0; c < 3; ++c)
          inside = inside && voxmm[c] >= 0 && voxmm[c] <= header.dim[c] * header.voxel_size[c];
        ras.push_back(p);
      }
    }
    return ras;
  }

  static void AssertSpecPointsMatchBundle(mitk::FiberBundle *fib, const std::vector<std::array<double, 3>> &ras, const std::string &message)
  {
    std::size_t k = 0;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j, ++k)
      {
        CPPUNIT_ASSERT_MESSAGE(message, k < ras.size());
        const double *p = cell->GetPoints()->GetPoint(j);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -p[0], ras[k][0], 1e-3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, -p[1], ras[k][1], 1e-3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message, p[2], ras[k][2], 1e-3);
      }
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, ras.size(), k);
  }

  /** The voxmm-to-LPS map with which MITK read and wrote .trk files from 2019 to 2026: vox_to_ras
   *  divided by the voxel size, rows negated by the voxel order letters R, A and I. */
  static vnl_matrix_fixed<double, 4, 4> Mitk2023Map(const TrkHeader &header)
  {
    vnl_matrix_fixed<double, 4, 4> map;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        map[i][j] = j < 3 ? header.vox_to_ras[i][j] / header.voxel_size[j] : header.vox_to_ras[i][j];
    const bool negate[3] = {header.voxel_order[0] == 'R', header.voxel_order[1] == 'A', header.voxel_order[2] == 'I'};
    for (int i = 0; i < 3; ++i)
      if (negate[i])
        for (int j = 0; j < 4; ++j)
          map[i][j] = -map[i][j];
    return map;
  }

  /** Writes the bundle as MITK did from 2019 to 2026 against the reference geometry: vox_to_ras its
   *  LPS index-to-world matrix, n_count 0, default image_orientation_patient. */
  static void WriteTrkAsMitk2023(const std::string &path, mitk::FiberBundle *fib, const mitk::BaseGeometry *geometry, const std::string &voxelOrder)
  {
    TrkHeader header;
    std::memset(&header, 0, sizeof(header));
    std::memcpy(header.id_string, "TRACK", 6);
    const vtkMatrix4x4 *matrix = geometry->GetVtkMatrix();
    for (int i = 0; i < 3; ++i)
    {
      header.dim[i] = static_cast<short>(geometry->GetExtent(i));
      header.voxel_size[i] = static_cast<float>(geometry->GetSpacing()[i]);
      header.origin[i] = static_cast<float>(geometry->GetOrigin()[i]);
    }
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        header.vox_to_ras[i][j] = static_cast<float>(matrix->GetElement(i, j));
    std::memcpy(header.voxel_order, voxelOrder.c_str(), 3);
    header.image_orientation_patient[0] = 1;
    header.image_orientation_patient[4] = 1;
    header.version = 2;
    header.hdr_size = 1000;

    const auto map = Mitk2023Map(header);
    const auto inverse = vnl_inverse(map.extract(3, 3));
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      const std::int32_t numPoints = static_cast<std::int32_t>(cell->GetNumberOfPoints());
      out.write(reinterpret_cast<const char *>(&numPoints), 4);
      for (vtkIdType j = 0; j < numPoints; ++j)
      {
        const double *p = cell->GetPoints()->GetPoint(j);
        const vnl_vector_fixed<double, 3> lps(p[0] - map[0][3], p[1] - map[1][3], p[2] - map[2][3]);
        const auto voxmm = inverse * lps;
        for (int c = 0; c < 3; ++c)
        {
          const float value = static_cast<float>(voxmm[c] + header.voxel_size[c] / 2);
          out.write(reinterpret_cast<const char *>(&value), 4);
        }
      }
    }
  }

  /** Reads a .trk file as MITK did from 2019 to 2026 (Mitk2023Map), as RAS points. */
  static std::vector<std::array<double, 3>> ReadTrkAsMitk2023(const std::string &path)
  {
    std::ifstream in(path, std::ios::binary);
    TrkHeader header;
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    const auto map = Mitk2023Map(header);
    std::vector<std::array<double, 3>> ras;
    std::int32_t numPoints = 0;
    while (in.read(reinterpret_cast<char *>(&numPoints), 4))
      for (int j = 0; j < numPoints; ++j)
      {
        float voxmm[3];
        in.read(reinterpret_cast<char *>(voxmm), 12);
        std::array<double, 3> p;
        for (int r = 0; r < 3; ++r)
        {
          p[r] = map[r][3];
          for (int c = 0; c < 3; ++c)
            p[r] += map[r][c] * (voxmm[c] - header.voxel_size[c] / 2);
        }
        ras.push_back({-p[0], -p[1], p[2]});
      }
    return ras;
  }

  void TrkWriter_LpsOptionAndOldMitkFiles_WriteConsistentHeaders()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/mask_image.nii.gz"));
    mitk::IFileWriter::Options ras;
    ras["Print header"] = false;
    mitk::IFileWriter::Options lps = ras;
    lps["Use LPS coordinate system (default is RAS)"] = true;
    const std::string path = m_TempDir + "/written.trk";
    TrkHeader header;
    bool inside = false;

    // "Use LPS" writes a grid whose voxel axes point to L, P and S
    mitk::IOUtil::Save(fib, path, lps);
    AssertSpecPointsMatchBundle(fib, ReadTrkAsSpec(path, header, inside), "Use LPS");
    CPPUNIT_ASSERT_MESSAGE("Use LPS: inside the declared volume", inside);
    CPPUNIT_ASSERT_EQUAL(std::string("LPS"), std::string(header.voxel_order, 3));
    CPPUNIT_ASSERT(header.vox_to_ras[0][0] < 0 && header.vox_to_ras[1][1] < 0 && header.vox_to_ras[2][2] > 0);

    // Files that MITK wrote from 2019 to 2026 against the image, by default ("RAS") and with "Use LPS"
    for (const std::string order : {"RAS", "LPS"})
    {
      const std::string old = m_TempDir + "/mitk2023_" + order + ".trk";
      WriteTrkAsMitk2023(old, fib, image->GetGeometry(), order);
      auto loaded = mitk::IOUtil::Load<mitk::FiberBundle>(old);
      CPPUNIT_ASSERT_MESSAGE(order, fib->Equals(loaded, 1e-3));

      // Saved again: nibabel and MITK 2023 to 2026 read the points, which lie in the image grid
      mitk::IOUtil::Save(loaded, path, ras);
      AssertSpecPointsMatchBundle(fib, ReadTrkAsSpec(path, header, inside), order + " saved again");
      CPPUNIT_ASSERT_MESSAGE(order + " saved again: inside the image grid", inside);
      CPPUNIT_ASSERT_EQUAL(std::string("RAS"), std::string(header.voxel_order, 3));
      AssertSpecPointsMatchBundle(fib, ReadTrkAsMitk2023(path), order + " saved again, read as MITK 2023 to 2026");

      mitk::IOUtil::Save(loaded, path, lps);
      AssertSpecPointsMatchBundle(fib, ReadTrkAsSpec(path, header, inside), order + " saved again with Use LPS");
      CPPUNIT_ASSERT_MESSAGE(order + " saved again with Use LPS: inside the image grid", inside);
    }
  }

  void TrkWriter_Header_DescribesReferenceGrid()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    mitk::IFileWriter::Options options;
    options["Print header"] = false;
    const std::string path = m_TempDir + "/spec.trk";
    TrkHeader header;
    bool inside = false;

    // Without reference image: a grid around the bundle
    mitk::IOUtil::Save(fib, path, options);
    AssertSpecPointsMatchBundle(fib, ReadTrkAsSpec(path, header, inside), "no reference");
    CPPUNIT_ASSERT_MESSAGE("All points inside the declared volume", inside);
    CPPUNIT_ASSERT_EQUAL(std::string("RAS"), std::string(header.voxel_order, 3));

    // Reference images stored in LPS and in LAS orientation
    for (const std::string name : {"StreamlineTractography/mask_image.nii.gz", "mean_FA.nii.gz"})
    {
      auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/" + name));
      auto withReference = fib->GetDeepCopy();
      withReference->SetTrackVisHeader(image->GetGeometry());
      mitk::IOUtil::Save(withReference, path, options);
      AssertSpecPointsMatchBundle(fib, ReadTrkAsSpec(path, header, inside), name);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, std::string("RAS"), std::string(header.voxel_order, 3));

      // The header grid is the image grid: every corner voxel of the image is a voxel of the header grid.
      vnl_matrix_fixed<double, 3, 3> a;
      vnl_vector_fixed<double, 3> t;
      for (int r = 0; r < 3; ++r)
      {
        t[r] = header.vox_to_ras[r][3];
        for (int c = 0; c < 3; ++c)
        {
          a[r][c] = header.vox_to_ras[r][c];
          if (r == c)
            CPPUNIT_ASSERT_MESSAGE(name + ": voxel axes point to R, A, S", a[r][c] > 0);
        }
      }
      const auto inverse = vnl_inverse(a);
      auto geometry = image->GetGeometry();
      for (int corner = 0; corner < 8; ++corner)
      {
        mitk::Point3D index;
        for (int c = 0; c < 3; ++c)
          index[c] = (corner >> c & 1) ? geometry->GetExtent(c) - 1 : 0;
        mitk::Point3D world;
        geometry->IndexToWorld(index, world);
        vnl_vector_fixed<double, 3> ras(-world[0], -world[1], world[2]);
        const auto headerIndex = inverse * (ras - t);
        for (int c = 0; c < 3; ++c)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(name, std::round(headerIndex[c]), headerIndex[c], 1e-3);
          CPPUNIT_ASSERT_MESSAGE(name, headerIndex[c] > -0.5 && headerIndex[c] < header.dim[c] - 0.5);
        }
      }
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleIO)
