/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionTestData.h>
#include <mitkFiberBundle.h>
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkTrxFile.h>

#include <itksys/SystemTools.hxx>
#include <nlohmann/json.hpp>
#include <Poco/Zip/ZipArchive.h>
#include <Poco/Zip/ZipLocalFileHeader.h>
#include <Poco/Zip/ZipStream.h>
#include <vtkCell.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  /** The values as little-endian bytes, as TRX stores them. */
  template <typename T>
  std::string Raw(const std::vector<T>& values)
  {
    std::string bytes;
    for (const T& value : values)
    {
      std::array<char, sizeof(T)> raw;
      std::memcpy(raw.data(), &value, sizeof(T));
      if constexpr (std::endian::native == std::endian::big)
        std::reverse(raw.begin(), raw.end());
      bytes.append(raw.data(), raw.size());
    }
    return bytes;
  }

  template <typename T>
  std::vector<T> ReadRaw(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<T> values(bytes.size() / sizeof(T));
    for (std::size_t i = 0; i < values.size(); ++i)
    {
      std::array<char, sizeof(T)> raw;
      std::memcpy(raw.data(), bytes.data() + i * sizeof(T), sizeof(T));
      if constexpr (std::endian::native == std::endian::big)
        std::reverse(raw.begin(), raw.end());
      std::memcpy(&values[i], raw.data(), sizeof(T));
    }
    return values;
  }

  using Point = std::array<double, 3>;
  using Fiber = std::vector<Point>;

  /** The fibers of ThirdPartyTractograms/README_v2.md in RAS millimetres; all coordinates are
   *  exact in float16. */
  const std::vector<Fiber> RasFibers = {{{1.5, 2.0, 3.0}, {4.0, 5.25, 5.0}, {7.0, 8.0, 6.5}},
                                        {{-1.0, -2.0, -3.0}, {-4.0, -5.0, 0.5}},
                                        {{0.0, 0.0, 0.0}, {0.5, -1.0, 1.0}, {1.0, -2.0, 2.0}, {1.5, -3.0, 3.0}}};

  const std::string ReferenceHeader =
    R"({"DIMENSIONS": [10, 12, 8], "VOXEL_TO_RASMM": [[-2.0, 0.0, 0.0, 9.0], [0.0, 2.0, 0.0, -11.0], [0.0, 0.0, 2.0, -7.0], [0.0, 0.0, 0.0, 1.0]], )";

  std::string Header(std::uint64_t streamlines, std::uint64_t vertices)
  {
    return ReferenceHeader + "\"NB_STREAMLINES\": " + std::to_string(streamlines) + ", \"NB_VERTICES\": " + std::to_string(vertices) + "}";
  }

  /** The fibers' points in file order and their offsets. */
  std::vector<double> FlatPoints(const std::vector<Fiber>& fibers)
  {
    std::vector<double> values;
    for (const auto& fiber : fibers)
      for (const auto& point : fiber)
        values.insert(values.end(), point.begin(), point.end());
    return values;
  }

  template <typename T>
  std::vector<T> Offsets(const std::vector<Fiber>& fibers)
  {
    std::vector<T> offsets = {0};
    for (const auto& fiber : fibers)
      offsets.push_back(static_cast<T>(offsets.back() + fiber.size()));
    return offsets;
  }

  template <typename T>
  std::vector<T> Cast(const std::vector<double>& values)
  {
    std::vector<T> result;
    for (double v : values)
      result.push_back(static_cast<T>(v));
    return result;
  }

  std::vector<std::uint16_t> Halves(const std::vector<double>& values)
  {
    std::vector<std::uint16_t> result;
    for (double v : values)
      result.push_back(mitk::TrxFile::FloatToHalf(static_cast<float>(v)));
    return result;
  }

  /** The points of fiber i of the bundle. */
  Fiber FiberPoints(mitk::FiberBundle* fib, unsigned int i)
  {
    Fiber points;
    vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
    for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
    {
      const double* p = cell->GetPoints()->GetPoint(j);
      points.push_back({p[0], p[1], p[2]});
    }
    return points;
  }

  /** The RGBA colour of point j of fiber i. */
  std::array<int, 4> FiberColor(mitk::FiberBundle* fib, unsigned int i, vtkIdType j)
  {
    const vtkIdType id = fib->GetFiberPolyData()->GetCell(i)->GetPointId(j);
    std::array<int, 4> rgba;
    for (int c = 0; c < 4; ++c)
      rgba[static_cast<std::size_t>(c)] = static_cast<int>(fib->GetFiberColors()->GetComponent(id, c));
    return rgba;
  }

  void AssertSamePoints(mitk::FiberBundle* expected, mitk::FiberBundle* actual, const std::string& message)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, expected->GetNumFibers(), actual->GetNumFibers());
    for (unsigned int i = 0; i < expected->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_MESSAGE(message + " fiber " + std::to_string(i), FiberPoints(expected, i) == FiberPoints(actual, i));
  }
}

/** The TRX reader and writer (mitkTrxFile.h): trx-python's test files and the generated ones in
 *  ThirdPartyTractograms/TRX of mitk-diffusion-test-data, and files written here. */
class mitkFiberBundleTrxTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleTrxTestSuite);
  MITK_TEST(HalfConversion_AllHalves_RoundTrip);
  MITK_TEST(HalfConversion_Rounding_TiesToEvenAndOverflow);
  MITK_TEST(Reader_GoldStandard_PointsAndColours);
  MITK_TEST(Reader_GoldStandard_HeaderDescribesReferenceImage);
  MITK_TEST(Reader_StoredDeflatedAndFolder_GiveTheSameBundle);
  MITK_TEST(Reader_Groups_OptionOffAndOn);
  MITK_TEST(Reader_PositionAndOffsetTypes_GiveTheSamePoints);
  MITK_TEST(Reader_GeneratedWeightsAndColours);
  MITK_TEST(Reader_StreamlinesWithoutPoints_DroppedWithTheirData);
  MITK_TEST(Reader_MalformedFiles_ThrowNamingTheFile);
  MITK_TEST(Writer_Float32RoundTrip_KeepsFibersWeightsColoursAndGrid);
  MITK_TEST(Writer_TrxFromTrx_KeepsVoxelToRasmmAndDimensions);
  MITK_TEST(Writer_Float16_WithinHalfAStep);
  MITK_TEST(Writer_Options_EntriesAndCompression);
  MITK_TEST(Writer_ForcedZip64_ReadBack);
  MITK_TEST(Writer_EmptyBundle_HeaderOnly);
  MITK_TEST(Writer_WriteErrors_Throw);
  CPPUNIT_TEST_SUITE_END();

  std::string m_TempDir;
  const std::string GROUPS_OPTION = "Load groups as separate bundles";

  static std::string Data(const std::string& name)
  {
    return mitk::GetDiffusionTestDataFilePath("ThirdPartyTractograms/TRX/" + name);
  }

  /** A TRX folder name.trx in the temporary directory with the given entries. */
  std::string WriteFolder(const std::string& name, const std::map<std::string, std::string>& entries)
  {
    const std::filesystem::path folder = std::filesystem::path(m_TempDir) / (name + ".trx");
    std::filesystem::remove_all(folder);
    std::filesystem::create_directories(folder);
    for (const auto& [entry, bytes] : entries)
    {
      const std::filesystem::path file = folder / entry;
      std::filesystem::create_directories(file.parent_path());
      std::ofstream(file, std::ios::binary) << bytes;
    }
    return folder.string();
  }

  /** The entries of a valid TRX folder with the fibers of RasFibers. */
  static std::map<std::string, std::string> ValidEntries()
  {
    return {{"header.json", Header(3, 9)},
            {"positions.3.float32", Raw(Cast<float>(FlatPoints(RasFibers)))},
            {"offsets.uint32", Raw(Offsets<std::uint32_t>(RasFibers))}};
  }

  std::vector<mitk::BaseData::Pointer> LoadAll(const std::string& path, bool groups)
  {
    mitk::IFileReader::Options options;
    options[GROUPS_OPTION] = groups;
    return mitk::IOUtil::Load(path, options);
  }

  void AssertRasFibers(mitk::FiberBundle* fib, const std::string& message)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<unsigned int>(RasFibers.size()), fib->GetNumFibers());
    for (unsigned int i = 0; i < RasFibers.size(); ++i)
    {
      Fiber expected;
      for (const auto& p : RasFibers[i])
        expected.push_back({-p[0], -p[1], p[2]});
      CPPUNIT_ASSERT_MESSAGE(message + " fiber " + std::to_string(i), expected == FiberPoints(fib, i));
    }
  }

  struct ZipEntry
  {
    int method;
    std::string content;
  };

  /** The entries of the ZIP archive at path with their compression method and content. */
  static std::map<std::string, ZipEntry> ZipEntries(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    Poco::Zip::ZipArchive archive(in);
    std::map<std::string, ZipEntry> entries;
    for (auto it = archive.headerBegin(); it != archive.headerEnd(); ++it)
    {
      std::ifstream file(path, std::ios::binary);
      Poco::Zip::ZipInputStream zipIn(file, it->second, true);
      const std::string content((std::istreambuf_iterator<char>(zipIn)), std::istreambuf_iterator<char>());
      entries[it->first] = {it->second.getCompressionMethod(), content};
    }
    return entries;
  }

  static std::set<std::string> Names(const std::map<std::string, ZipEntry>& entries)
  {
    std::set<std::string> names;
    for (const auto& entry : entries)
      names.insert(entry.first);
    return names;
  }

  /** original.fib with a weight per fiber, an RGBA colour per point (alpha 100, which TRX does not
   *  keep) and the oblique reference grid of gs.trx */
  mitk::FiberBundle::Pointer AnnotatedBundle()
  {
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("FiberProcessing/original.fib"));
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      fib->SetFiberWeight(i, 0.25f + 0.125f * static_cast<float>(i));
    auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(4);
    colors->SetNumberOfTuples(fib->GetFiberPolyData()->GetNumberOfPoints());
    for (vtkIdType k = 0; k < colors->GetNumberOfTuples(); ++k)
    {
      const unsigned char rgba[4] = {static_cast<unsigned char>(k % 251), static_cast<unsigned char>((3 * k) % 256), static_cast<unsigned char>((7 * k) % 256), 100};
      colors->SetTypedTuple(k, rgba);
    }
    fib->SetFiberColors(colors);
    fib->SetTrackVisHeader(mitk::IOUtil::Load<mitk::FiberBundle>(Data("gs.trx"))->GetTrackVisHeader());
    return fib;
  }

  static mitk::IFileWriter::Options WriterOptions(const std::string& positionType, const std::string& compression, bool weights, bool colors)
  {
    mitk::IFileWriter::Options options;
    options["Position type"] = positionType;
    options["Compression"] = compression;
    options["Save fiber weights"] = weights;
    options["Save color information"] = colors;
    return options;
  }

  static void AssertSameGrid(mitk::FiberBundle* expected, mitk::FiberBundle* actual)
  {
    const auto e = expected->GetTrackVisHeader();
    const auto a = actual->GetTrackVisHeader();
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c)
        CPPUNIT_ASSERT_EQUAL(e.vox_to_ras[r][c], a.vox_to_ras[r][c]);
    for (int c = 0; c < 3; ++c)
      CPPUNIT_ASSERT_EQUAL(e.dim[c], a.dim[c]);
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("mitkFiberBundleTrxTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void HalfConversion_AllHalves_RoundTrip()
  {
    for (std::uint32_t bits = 0; bits <= 0xffff; ++bits)
    {
      const auto half = static_cast<std::uint16_t>(bits);
      const float value = mitk::TrxFile::HalfToFloat(half);
      const bool isNan = (half & 0x7c00) == 0x7c00 && (half & 0x3ff) != 0;
      CPPUNIT_ASSERT_EQUAL_MESSAGE(std::to_string(bits), isNan, static_cast<bool>(std::isnan(value)));
      if (isNan)
        CPPUNIT_ASSERT_MESSAGE(std::to_string(bits), std::isnan(mitk::TrxFile::HalfToFloat(mitk::TrxFile::FloatToHalf(value))));
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE(std::to_string(bits), half, mitk::TrxFile::FloatToHalf(value));
    }
    // a few values by definition
    CPPUNIT_ASSERT_EQUAL(1.0f, mitk::TrxFile::HalfToFloat(0x3c00));
    CPPUNIT_ASSERT_EQUAL(-2.0f, mitk::TrxFile::HalfToFloat(0xc000));
    CPPUNIT_ASSERT_EQUAL(65504.0f, mitk::TrxFile::HalfToFloat(0x7bff));
    CPPUNIT_ASSERT_EQUAL(std::ldexp(1.0f, -24), mitk::TrxFile::HalfToFloat(0x0001));
    CPPUNIT_ASSERT_EQUAL(std::ldexp(1.0f, -14), mitk::TrxFile::HalfToFloat(0x0400));
    CPPUNIT_ASSERT(std::isinf(mitk::TrxFile::HalfToFloat(0x7c00)));
  }

  void HalfConversion_Rounding_TiesToEvenAndOverflow()
  {
    using mitk::TrxFile::FloatToHalf;
    // between 1 and the next half (1 + 2^-10): halfway goes to the even one (1), above to the next
    const float ulp = std::ldexp(1.0f, -10);
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x3c00), FloatToHalf(1.0f + ulp / 2));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x3c01), FloatToHalf(std::nextafter(1.0f + ulp / 2, 2.0f)));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x3c00), FloatToHalf(std::nextafter(1.0f + ulp / 2, 0.0f)));
    // halfway between 1 + 2^-10 (odd) and 1 + 2^-9 (even): up
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x3c02), FloatToHalf(1.0f + 1.5f * ulp));
    // subnormals: halfway between 0 and 2^-24 goes to 0, between 2^-24 and 2^-23 to 2^-23
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x0000), FloatToHalf(std::ldexp(1.0f, -25)));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x0001), FloatToHalf(std::nextafter(std::ldexp(1.0f, -25), 1.0f)));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x0002), FloatToHalf(1.5f * std::ldexp(1.0f, -24)));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x8000), FloatToHalf(-std::ldexp(1.0f, -26)));
    // the largest subnormal rounds up to the smallest normal half
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x0400), FloatToHalf(std::nextafter(std::ldexp(1.0f, -14), 0.0f)));
    // overflow: 65519 rounds to 65504, 65520 (halfway, odd mantissa) and beyond to infinity
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x7bff), FloatToHalf(65519.0f));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x7c00), FloatToHalf(65520.0f));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0xfc00), FloatToHalf(-1e10f));
    CPPUNIT_ASSERT_EQUAL(std::uint16_t(0x7c00), FloatToHalf(std::numeric_limits<float>::infinity()));
    // every float within half an ulp of the half it gives
    for (float v = -70000.0f; v < 70000.0f; v += 0.37f)
    {
      const float rounded = mitk::TrxFile::HalfToFloat(FloatToHalf(v));
      if (std::isinf(rounded))
        continue;
      const float next = mitk::TrxFile::HalfToFloat(static_cast<std::uint16_t>(FloatToHalf(std::fabs(v)) + 1));
      CPPUNIT_ASSERT_MESSAGE(std::to_string(v), std::fabs(rounded - v) <= (next - std::fabs(rounded)) / 2);
    }
  }

  void Reader_GoldStandard_PointsAndColours()
  {
    // gs_rasmm_space.txt holds the points of gs.trx in RAS millimetres, one per line.
    std::vector<Point> ras;
    std::ifstream text(Data("gs_rasmm_space.txt"));
    Point p;
    while (text >> p[0] >> p[1] >> p[2])
      ras.push_back(p);
    CPPUNIT_ASSERT_EQUAL(std::size_t(104), ras.size());
    const auto red = ReadRaw<float>(Data("gs_fldr.trx/dpv/color_x.float32"));
    const auto green = ReadRaw<float>(Data("gs_fldr.trx/dpv/color_y.float32"));
    const auto blue = ReadRaw<float>(Data("gs_fldr.trx/dpv/color_z.float32"));
    CPPUNIT_ASSERT_EQUAL(std::size_t(104), red.size());

    for (const std::string name : {"gs.trx", "gs_fldr.trx"})
    {
      auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(Data(name));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 13u, fib->GetNumFibers());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 104u, fib->GetNumberOfPoints());
      std::size_t k = 0;
      for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      {
        const Fiber points = FiberPoints(fib, i);
        for (std::size_t j = 0; j < points.size(); ++j, ++k)
        {
          const std::string where = name + " point " + std::to_string(k);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(where, -ras[k][0], points[j][0], 1e-5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(where, -ras[k][1], points[j][1], 1e-5);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(where, ras[k][2], points[j][2], 1e-5);
          const std::array<int, 4> expected = {static_cast<int>(red[k]), static_cast<int>(green[k]), static_cast<int>(blue[k]), 255};
          CPPUNIT_ASSERT_MESSAGE(where, expected == FiberColor(fib, i, static_cast<vtkIdType>(j)));
        }
      }
      // weights: the file has none (dps/random_coord is dropped)
      for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 1.0f, fib->GetFiberWeight(i));
    }
  }

  void Reader_GoldStandard_HeaderDescribesReferenceImage()
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(Data("gs.nii"));
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(Data("gs.trx"));
    CPPUNIT_ASSERT(fib->HasTrackVisHeader());
    const auto header = fib->GetTrackVisHeader();
    for (int c = 0; c < 3; ++c)
      CPPUNIT_ASSERT_EQUAL(static_cast<int>(image->GetDimension(static_cast<unsigned int>(c))), static_cast<int>(header.dim[c]));
    // index -> LPS of the header (vox_to_ras with x and y negated) and of the image, at the
    // corners of the grid; the header holds float32 values (as trx-python reads them).
    for (int corner = 0; corner < 8; ++corner)
    {
      mitk::Point3D index;
      for (int c = 0; c < 3; ++c)
        index[c] = (corner >> c & 1) ? header.dim[c] - 1 : 0;
      mitk::Point3D world;
      image->GetGeometry()->IndexToWorld(index, world);
      for (int r = 0; r < 3; ++r)
      {
        double ras = header.vox_to_ras[r][3];
        for (int c = 0; c < 3; ++c)
          ras += header.vox_to_ras[r][c] * index[c];
        CPPUNIT_ASSERT_DOUBLES_EQUAL(world[r], r < 2 ? -ras : ras, 1e-4);
      }
    }
  }

  void Reader_StoredDeflatedAndFolder_GiveTheSameBundle()
  {
    // float16 positions, uint64 offsets; small_compressed.trx is deflated
    auto stored = mitk::IOUtil::Load<mitk::FiberBundle>(Data("small.trx"));
    CPPUNIT_ASSERT_EQUAL(1000u, stored->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(33886u, stored->GetNumberOfPoints());
    const Fiber first = FiberPoints(stored, 0);
    // file: (20.28, -33.44, -41.9) in RAS, as float16
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-20.28, first[0][0], 0.01);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(33.44, first[0][1], 0.01);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-41.9, first[0][2], 0.01);
    CPPUNIT_ASSERT(118 == FiberColor(stored, 0, 0)[0]);
    // dps/commit_weights is no fiber weight
    for (unsigned int i = 0; i < stored->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(1.0f, stored->GetFiberWeight(i));

    // The header keeps the LAS axes of VOXEL_TO_RASMM.
    const auto header = stored->GetTrackVisHeader();
    CPPUNIT_ASSERT_EQUAL(std::string("LAS"), std::string(header.voxel_order, 3));
    const float matrix[4][4] = {{-1.25f, 0, 0, 72.5f}, {0, 1.25f, 0, -109.75f}, {0, 0, 1.25f, -64.5f}, {0, 0, 0, 1}};
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c)
        CPPUNIT_ASSERT_EQUAL(matrix[r][c], header.vox_to_ras[r][c]);
    CPPUNIT_ASSERT_EQUAL(short(117), header.dim[0]);
    CPPUNIT_ASSERT_EQUAL(short(151), header.dim[1]);
    CPPUNIT_ASSERT_EQUAL(short(115), header.dim[2]);

    for (const std::string name : {"small_compressed.trx", "small_fldr.trx"})
    {
      auto other = mitk::IOUtil::Load<mitk::FiberBundle>(Data(name));
      AssertSamePoints(stored, other, name);
      for (unsigned int i = 0; i < stored->GetNumFibers(); ++i)
      {
        const auto numPoints = static_cast<vtkIdType>(FiberPoints(stored, i).size());
        for (vtkIdType j = 0; j < numPoints; ++j)
          CPPUNIT_ASSERT_MESSAGE(name, FiberColor(stored, i, j) == FiberColor(other, i, j));
      }
      const auto otherHeader = other->GetTrackVisHeader();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 0, std::memcmp(&header, &otherHeader, sizeof(header)));
    }
  }

  void Reader_Groups_OptionOffAndOn()
  {
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), LoadAll(Data("small.trx"), false).size());

    const auto data = LoadAll(Data("small.trx"), true);
    const std::vector<std::string> names = {"AF_L", "AF_R", "CC", "CST_L", "CST_R"};
    const std::vector<unsigned int> counts = {2, 5, 34, 14, 9};
    CPPUNIT_ASSERT_EQUAL(std::size_t(6), data.size());
    auto whole = dynamic_cast<mitk::FiberBundle*>(data[0].GetPointer());
    CPPUNIT_ASSERT(whole != nullptr);
    CPPUNIT_ASSERT_EQUAL(1000u, whole->GetNumFibers());
    for (std::size_t g = 0; g < names.size(); ++g)
    {
      auto group = dynamic_cast<mitk::FiberBundle*>(data[g + 1].GetPointer());
      CPPUNIT_ASSERT(group != nullptr);
      CPPUNIT_ASSERT_EQUAL(names[g], group->GetProperty("name")->GetValueAsString());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(names[g], counts[g], group->GetNumFibers());
      const auto wholeHeader = whole->GetTrackVisHeader();
      const auto groupHeader = group->GetTrackVisHeader();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(names[g], 0, std::memcmp(&wholeHeader, &groupHeader, sizeof(wholeHeader)));
      // the fibers of the group file, in its order (with repeated indices)
      const auto indices = ReadRaw<std::int64_t>(Data("small_fldr.trx/groups/" + names[g] + ".int64"));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(names[g], indices.size(), std::size_t(group->GetNumFibers()));
      for (std::size_t j = 0; j < indices.size(); ++j)
      {
        const auto i = static_cast<unsigned int>(indices[j]);
        CPPUNIT_ASSERT_MESSAGE(names[g], FiberPoints(whole, i) == FiberPoints(group, static_cast<unsigned int>(j)));
        CPPUNIT_ASSERT_MESSAGE(names[g], FiberColor(whole, i, 0) == FiberColor(group, static_cast<unsigned int>(j), 0));
      }
    }
  }

  void Reader_PositionAndOffsetTypes_GiveTheSamePoints()
  {
    const auto values = FlatPoints(RasFibers);
    const std::map<std::string, std::string> positions = {{"positions.3.float16", Raw(Halves(values))},
                                                          {"positions.3.float32", Raw(Cast<float>(values))},
                                                          {"positions.3.float64", Raw(values)}};
    const std::map<std::string, std::string> offsets = {{"offsets.uint32", Raw(Offsets<std::uint32_t>(RasFibers))},
                                                        {"offsets.uint64", Raw(Offsets<std::uint64_t>(RasFibers))}};
    for (const auto& [positionsName, positionsBytes] : positions)
      for (const auto& [offsetsName, offsetsBytes] : offsets)
      {
        const std::string name = positionsName + "_" + offsetsName;
        const auto folder = WriteFolder(name, {{"header.json", Header(3, 9)}, {positionsName, positionsBytes}, {offsetsName, offsetsBytes}});
        AssertRasFibers(mitk::IOUtil::Load<mitk::FiberBundle>(folder), name);
      }
  }

  void Reader_GeneratedWeightsAndColours()
  {
    // trx_python.trx: trx-python's defaults; weights_colors.trx: float64 positions, uint64
    // offsets, dps/weights, dpv/color and a group of streamlines 0 and 2
    AssertRasFibers(mitk::IOUtil::Load<mitk::FiberBundle>(Data("trx_python.trx")), "trx_python.trx");
    const auto data = LoadAll(Data("weights_colors.trx"), true);
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), data.size());
    auto fib = dynamic_cast<mitk::FiberBundle*>(data[0].GetPointer());
    AssertRasFibers(fib, "weights_colors.trx");
    const std::vector<float> weights = {0.5f, 1.0f, 2.0f};
    int k = 0;
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(weights[i], fib->GetFiberWeight(i));
      for (vtkIdType j = 0; j < static_cast<vtkIdType>(RasFibers[i].size()); ++j, ++k)
      {
        const std::array<int, 4> expected = {20 * k, 100 + 10 * k, 250 - 25 * k, 255};
        CPPUNIT_ASSERT(expected == FiberColor(fib, i, j));
      }
    }
    auto group = dynamic_cast<mitk::FiberBundle*>(data[1].GetPointer());
    CPPUNIT_ASSERT_EQUAL(std::string("first_and_last"), group->GetProperty("name")->GetValueAsString());
    CPPUNIT_ASSERT_EQUAL(2u, group->GetNumFibers());
    CPPUNIT_ASSERT(FiberPoints(fib, 0) == FiberPoints(group, 0));
    CPPUNIT_ASSERT(FiberPoints(fib, 2) == FiberPoints(group, 1));
    CPPUNIT_ASSERT_EQUAL(0.5f, group->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(2.0f, group->GetFiberWeight(1));
    CPPUNIT_ASSERT(FiberColor(fib, 2, 3) == FiberColor(group, 1, 3));
  }

  void Reader_StreamlinesWithoutPoints_DroppedWithTheirData()
  {
    // streamlines of 3, 0, 2 and 1 points; weights 1 to 4; point k has colour (k, 2k, 3k); a group
    // of streamlines 1, 2 and 3
    std::vector<float> positions;
    std::vector<std::uint8_t> colors;
    for (int k = 0; k < 6; ++k)
    {
      positions.insert(positions.end(), {static_cast<float>(k), 0.5f * k, 1.0f});
      colors.insert(colors.end(), {static_cast<std::uint8_t>(k), static_cast<std::uint8_t>(2 * k), static_cast<std::uint8_t>(3 * k)});
    }
    const auto folder = WriteFolder("empty_streamline", {{"header.json", Header(4, 6)},
                                                         {"positions.3.float32", Raw(positions)},
                                                         {"offsets.uint64", Raw(std::vector<std::uint64_t>{0, 3, 3, 5, 6})},
                                                         {"dps/weights.float32", Raw(std::vector<float>{1, 2, 3, 4})},
                                                         {"dpv/color.3.uint8", Raw(colors)},
                                                         {"groups/later.uint32", Raw(std::vector<std::uint32_t>{1, 2, 3})}});
    const auto data = LoadAll(folder, true);
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), data.size());
    auto fib = dynamic_cast<mitk::FiberBundle*>(data[0].GetPointer());
    auto group = dynamic_cast<mitk::FiberBundle*>(data[1].GetPointer());
    CPPUNIT_ASSERT_EQUAL(3u, fib->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(2u, group->GetNumFibers());
    const std::vector<float> weights = {1, 3, 4};
    const std::vector<std::vector<int>> pointsOfFibers = {{0, 1, 2}, {3, 4}, {5}};
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(weights[i], fib->GetFiberWeight(i));
      const Fiber points = FiberPoints(fib, i);
      CPPUNIT_ASSERT_EQUAL(pointsOfFibers[i].size(), points.size());
      for (std::size_t j = 0; j < points.size(); ++j)
      {
        const int k = pointsOfFibers[i][j];
        CPPUNIT_ASSERT_EQUAL(-static_cast<double>(k), points[j][0]);
        const std::array<int, 4> expected = {k, 2 * k, 3 * k, 255};
        CPPUNIT_ASSERT(expected == FiberColor(fib, i, static_cast<vtkIdType>(j)));
      }
    }
    CPPUNIT_ASSERT(FiberPoints(fib, 1) == FiberPoints(group, 0));
    CPPUNIT_ASSERT(FiberPoints(fib, 2) == FiberPoints(group, 1));
    CPPUNIT_ASSERT_EQUAL(3.0f, group->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(4.0f, group->GetFiberWeight(1));
    CPPUNIT_ASSERT(FiberColor(fib, 2, 0) == FiberColor(group, 1, 0));
  }

  void Reader_MalformedFiles_ThrowNamingTheFile()
  {
    const auto valid = ValidEntries();
    auto with = [&valid](const std::string& name, const std::string& bytes)
    {
      auto entries = valid;
      entries[name] = bytes;
      return entries;
    };
    auto without = [&valid](const std::string& name)
    {
      auto entries = valid;
      entries.erase(name);
      return entries;
    };
    auto replaced = [&valid](const std::string& oldName, const std::string& name, const std::string& bytes)
    {
      auto entries = valid;
      entries.erase(oldName);
      entries[name] = bytes;
      return entries;
    };
    std::vector<float> nan = Cast<float>(FlatPoints(RasFibers));
    nan[4] = std::numeric_limits<float>::quiet_NaN();
    std::vector<double> beyondFloat = FlatPoints(RasFibers);
    beyondFloat[7] = 1e300;

    struct Case
    {
      std::string name;
      std::map<std::string, std::string> entries;
      /** whether groups are loaded */
      bool groups;
      /** a part of the expected message */
      std::string problem;
    };
    const std::string positions = "positions.3.float32";
    const std::string offsets = "offsets.uint32";
    const std::vector<Case> cases = {
      {"no_header", without("header.json"), false, "header.json is missing"},
      {"invalid_json", with("header.json", "{\"NB_STREAMLINES\": 3,"), false, "no valid JSON"},
      {"no_count", with("header.json", ReferenceHeader + "\"NB_STREAMLINES\": 3}"), false, "header.json has no NB_VERTICES"},
      {"negative_count", with("header.json", ReferenceHeader + "\"NB_STREAMLINES\": -3, \"NB_VERTICES\": 9}"), false, "NB_STREAMLINES in header.json is no count"},
      {"no_positions", without(positions), false, "positions array is missing"},
      {"no_offsets", without(offsets), false, "offsets array is missing"},
      {"positions_int", replaced(positions, "positions.3.int16", Raw(std::vector<std::int16_t>(27, 1))), false, "positions must be"},
      {"positions_components", replaced(positions, "positions.float32", Raw(Cast<float>(FlatPoints(RasFibers)))), false, "positions must be"},
      {"positions_dtype", replaced(positions, "positions.3.float128", ""), false, "unsupported name or data type"},
      {"two_positions", with("positions.3.float64", Raw(FlatPoints(RasFibers))), false, "two positions arrays"},
      {"offsets_signed", replaced(offsets, "offsets.int32", Raw(Offsets<std::int32_t>(RasFibers))), false, "uint32 or uint64"},
      {"offsets_count", with(offsets, Raw(std::vector<std::uint32_t>{0, 3, 9})), false, "NB_STREAMLINES + 1"},
      {"offsets_first", with(offsets, Raw(std::vector<std::uint32_t>{1, 3, 5, 9})), false, "first offset"},
      {"offsets_decrease", with(offsets, Raw(std::vector<std::uint32_t>{0, 5, 3, 9})), false, "decrease"},
      {"offsets_last", with(offsets, Raw(std::vector<std::uint32_t>{0, 3, 5, 8})), false, "last offset"},
      {"positions_size", with(positions, Raw(std::vector<float>(24, 1.0f))), false, "NB_VERTICES"},
      {"positions_nan", with(positions, Raw(nan)), false, "not finite"},
      {"positions_beyond_float", replaced(positions, "positions.3.float64", Raw(beyondFloat)), false, "not finite"},
      {"dps_rows", with("dps/length.float32", Raw(std::vector<float>{1, 2})), false, "NB_STREAMLINES"},
      {"dpv_rows", with("dpv/color.3.uint8", Raw(std::vector<std::uint8_t>(24, 1))), false, "NB_VERTICES"},
      {"group_range", with("groups/g.uint32", Raw(std::vector<std::uint32_t>{0, 3})), true, "streamline index 3"},
      {"group_negative", with("groups/g.int64", Raw(std::vector<std::int64_t>{-1})), true, "streamline index -1"},
    };
    for (const auto& c : cases)
      AssertThrows(WriteFolder(c.name, c.entries), c.groups, c.problem);

    // a group index out of range does no harm while groups are not loaded
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), LoadAll(WriteFolder("group_not_loaded", with("groups/g.uint32", Raw(std::vector<std::uint32_t>{0, 3}))), false).size());

    // no ZIP archive, and a truncated one
    const std::string garbage = m_TempDir + "/garbage.trx";
    std::ofstream(garbage, std::ios::binary) << "This is no ZIP archive.";
    AssertThrows(garbage, false, "is no TRX file");
    std::ifstream in(Data("small.trx"), std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string truncated = m_TempDir + "/truncated.trx";
    std::ofstream(truncated, std::ios::binary) << bytes.substr(0, bytes.size() / 2);
    AssertThrows(truncated, false, "");
  }

  void Writer_Float32RoundTrip_KeepsFibersWeightsColoursAndGrid()
  {
    auto fib = AnnotatedBundle();
    const std::string path = m_TempDir + "/roundtrip.trx";
    mitk::IOUtil::Save(fib, path);
    auto read = mitk::IOUtil::Load<mitk::FiberBundle>(path);
    AssertSamePoints(fib, read, "float32");
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberWeight(i), read->GetFiberWeight(i));
      const auto numPoints = static_cast<vtkIdType>(FiberPoints(fib, i).size());
      for (vtkIdType j = 0; j < numPoints; ++j)
      {
        auto expected = FiberColor(fib, i, j);
        expected[3] = 255;
        CPPUNIT_ASSERT(expected == FiberColor(read, i, j));
      }
    }
    AssertSameGrid(fib, read);
  }

  void Writer_TrxFromTrx_KeepsVoxelToRasmmAndDimensions()
  {
    // the LAS grid of small.trx and the oblique one of gs.trx
    for (const std::string name : {"small.trx", "gs.trx"})
    {
      auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(Data(name));
      const std::string path = m_TempDir + "/copy_of_" + name;
      mitk::IOUtil::Save(fib, path);
      const auto original = nlohmann::json::parse(ZipEntries(Data(name)).at("header.json").content);
      const auto written = nlohmann::json::parse(ZipEntries(path).at("header.json").content);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, original["VOXEL_TO_RASMM"].dump(), written["VOXEL_TO_RASMM"].dump());
      CPPUNIT_ASSERT_MESSAGE(name, original["VOXEL_TO_RASMM"] == written["VOXEL_TO_RASMM"]);
      CPPUNIT_ASSERT_MESSAGE(name, original["DIMENSIONS"] == written["DIMENSIONS"]);
      CPPUNIT_ASSERT_MESSAGE(name, original["NB_STREAMLINES"] == written["NB_STREAMLINES"]);
      CPPUNIT_ASSERT_MESSAGE(name, original["NB_VERTICES"] == written["NB_VERTICES"]);
      AssertSamePoints(fib, mitk::IOUtil::Load<mitk::FiberBundle>(path), name);
    }
  }

  void Writer_Float16_WithinHalfAStep()
  {
    auto fib = AnnotatedBundle();
    const std::string path = m_TempDir + "/half.trx";
    mitk::IOUtil::Save(fib, path, WriterOptions("float16", "stored", true, true));
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), ZipEntries(path).count("positions.3.float16"));
    auto read = mitk::IOUtil::Load<mitk::FiberBundle>(path);
    CPPUNIT_ASSERT_EQUAL(fib->GetNumFibers(), read->GetNumFibers());
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      const Fiber expected = FiberPoints(fib, i);
      const Fiber actual = FiberPoints(read, i);
      CPPUNIT_ASSERT_EQUAL(expected.size(), actual.size());
      for (std::size_t j = 0; j < expected.size(); ++j)
        for (std::size_t c = 0; c < 3; ++c)
        {
          // half of the float16 step at the value: 2^(exponent - 10) / 2, at least the subnormal 2^-25
          const double v = std::fabs(expected[j][c]);
          const double halfStep = v < std::ldexp(1.0, -14) ? std::ldexp(1.0, -25) : std::ldexp(1.0, std::ilogb(v) - 11);
          CPPUNIT_ASSERT(std::fabs(expected[j][c] - actual[j][c]) <= halfStep);
          CPPUNIT_ASSERT_EQUAL(mitk::TrxFile::HalfToFloat(mitk::TrxFile::FloatToHalf(static_cast<float>(actual[j][c]))), static_cast<float>(actual[j][c]));
        }
    }
  }

  void Writer_Options_EntriesAndCompression()
  {
    auto fib = AnnotatedBundle();
    const std::string stored = m_TempDir + "/stored.trx";
    mitk::IOUtil::Save(fib, stored);
    auto entries = ZipEntries(stored);
    const std::set<std::string> all = {"header.json", "positions.3.float32", "offsets.uint32", "dps/weights.float32", "dpv/color.3.uint8"};
    CPPUNIT_ASSERT(all == Names(entries));
    for (const auto& entry : entries)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(entry.first, int(Poco::Zip::ZipCommon::CM_STORE), entry.second.method);
    // offsets with the last one, NB_VERTICES
    std::vector<std::uint32_t> offsets = {0};
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      offsets.push_back(offsets.back() + static_cast<std::uint32_t>(FiberPoints(fib, i).size()));
    CPPUNIT_ASSERT(Raw(offsets) == entries["offsets.uint32"].content);
    std::vector<float> weights;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      weights.push_back(fib->GetFiberWeight(i));
    CPPUNIT_ASSERT(Raw(weights) == entries["dps/weights.float32"].content);

    const std::string deflated = m_TempDir + "/deflated.trx";
    mitk::IOUtil::Save(fib, deflated, WriterOptions("float16", "deflated", false, false));
    entries = ZipEntries(deflated);
    const std::set<std::string> some = {"header.json", "positions.3.float16", "offsets.uint32"};
    CPPUNIT_ASSERT(some == Names(entries));
    for (const auto& entry : entries)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(entry.first, int(Poco::Zip::ZipCommon::CM_DEFLATE), entry.second.method);
    auto read = mitk::IOUtil::Load<mitk::FiberBundle>(deflated);
    CPPUNIT_ASSERT_EQUAL(fib->GetNumFibers(), read->GetNumFibers());
    for (unsigned int i = 0; i < read->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(1.0f, read->GetFiberWeight(i));
    AssertSameGrid(fib, read);

    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(fib, m_TempDir + "/wrong.trx", WriterOptions("float64", "stored", true, true)), mitk::Exception);
  }

  void Writer_ForcedZip64_ReadBack()
  {
    auto fib = AnnotatedBundle();
    const std::string path = m_TempDir + "/zip64.trx";
    mitk::TrxFile::WriteOptions options;
    options.zip64Threshold = 0;
    mitk::TrxFile::Write(path, fib, options);
    std::ifstream in(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    // the Zip64 end of central directory record and its locator
    CPPUNIT_ASSERT(bytes.find(std::string("PK\x06\x06", 4)) != std::string::npos);
    CPPUNIT_ASSERT(bytes.find(std::string("PK\x06\x07", 4)) != std::string::npos);
    auto read = mitk::IOUtil::Load<mitk::FiberBundle>(path);
    AssertSamePoints(fib, read, "zip64");
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberWeight(i), read->GetFiberWeight(i));
  }

  void Writer_EmptyBundle_HeaderOnly()
  {
    // as trx-python writes a file without streamlines
    const std::string path = m_TempDir + "/empty.trx";
    mitk::IOUtil::Save(mitk::FiberBundle::New(), path);
    CPPUNIT_ASSERT(std::set<std::string>{"header.json"} == Names(ZipEntries(path)));
    CPPUNIT_ASSERT_EQUAL(0u, mitk::IOUtil::Load<mitk::FiberBundle>(path)->GetNumFibers());
  }

  void Writer_WriteErrors_Throw()
  {
    auto fib = AnnotatedBundle();
    CPPUNIT_ASSERT_THROW(mitk::IOUtil::Save(fib, m_TempDir + "/no/such/folder/fibers.trx"), mitk::Exception);
    mitk::TrxFile::WriteOptions options;
    CPPUNIT_ASSERT_THROW(mitk::TrxFile::Write(m_TempDir, fib, options), mitk::Exception);
  }

  /** Loading path throws; the reader's message names path (IOUtil names it once more) and
   *  contains problem. */
  void AssertThrows(const std::string& path, bool groups, const std::string& problem)
  {
    try
    {
      LoadAll(path, groups);
    }
    catch (const mitk::Exception& e)
    {
      const std::string message = e.what();
      std::size_t count = 0;
      for (auto pos = message.find(path); pos != std::string::npos; pos = message.find(path, pos + 1))
        ++count;
      CPPUNIT_ASSERT_MESSAGE(path + ": " + message, count >= 2);
      CPPUNIT_ASSERT_MESSAGE(path + ": " + message, message.find(problem) != std::string::npos);
      return;
    }
    CPPUNIT_FAIL(path + " loaded without an exception");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleTrx)
