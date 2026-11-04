/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkTrxFile.h"
#include "mitkFiberPolyDataBuilder.h"
#include "mitkTrackVisFile.h"

#include <mitkExceptionMacro.h>
#include <mitkLog.h>

#include <nlohmann/json.hpp>
#include <Poco/DateTime.h>
#include <Poco/Exception.h>
#include <Poco/Path.h>
#include <Poco/Zip/Compress.h>
#include <Poco/Zip/ZipArchive.h>
#include <Poco/Zip/ZipLocalFileHeader.h>
#include <Poco/Zip/ZipStream.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <streambuf>
#include <type_traits>

namespace
{
  enum class DataType { Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64, Float16, Float32, Float64, Bit };

  struct DataTypeName
  {
    const char* name;
    DataType type;
    std::size_t size;
  };

  /** The data types of TRX file names; trx-python stores "bit" (numpy bool) as one byte. */
  constexpr std::array<DataTypeName, 12> DataTypeNames = {{
    {"int8", DataType::Int8, 1}, {"int16", DataType::Int16, 2}, {"int32", DataType::Int32, 4}, {"int64", DataType::Int64, 8},
    {"uint8", DataType::UInt8, 1}, {"uint16", DataType::UInt16, 2}, {"uint32", DataType::UInt32, 4}, {"uint64", DataType::UInt64, 8},
    {"float16", DataType::Float16, 2}, {"float32", DataType::Float32, 4}, {"float64", DataType::Float64, 8}, {"bit", DataType::Bit, 1}}};

  std::size_t SizeOf(DataType type)
  {
    for (const auto& entry : DataTypeNames)
      if (entry.type == type)
        return entry.size;
    return 1;
  }

  bool IsFloat(DataType type)
  {
    return type == DataType::Float16 || type == DataType::Float32 || type == DataType::Float64;
  }

  bool IsIndexType(DataType type)
  {
    return !IsFloat(type) && type != DataType::Bit;
  }

  /** A half-precision value as the file stores it. */
  struct Half
  {
    std::uint16_t bits;
  };

  /** The value of type T stored little-endian at bytes, as TRX stores all arrays. */
  template <typename T>
  T LittleEndian(const char* bytes)
  {
    std::array<char, sizeof(T)> raw;
    std::memcpy(raw.data(), bytes, sizeof(T));
    if constexpr (std::endian::native == std::endian::big)
      std::reverse(raw.begin(), raw.end());
    T value;
    std::memcpy(&value, raw.data(), sizeof(T));
    return value;
  }

  /** The value of file type S at bytes; halves as float. */
  template <typename S>
  auto Decode(const char* bytes)
  {
    if constexpr (std::is_same_v<S, Half>)
      return mitk::TrxFile::HalfToFloat(LittleEndian<std::uint16_t>(bytes));
    else
      return LittleEndian<S>(bytes);
  }

  /** Calls f with a value of the C++ type of the data type (Half for float16, one byte for bit). */
  template <typename F>
  void WithType(DataType type, F&& f)
  {
    switch (type)
    {
      case DataType::Int8: f(std::int8_t()); break;
      case DataType::Int16: f(std::int16_t()); break;
      case DataType::Int32: f(std::int32_t()); break;
      case DataType::Int64: f(std::int64_t()); break;
      case DataType::UInt8: f(std::uint8_t()); break;
      case DataType::UInt16: f(std::uint16_t()); break;
      case DataType::UInt32: f(std::uint32_t()); break;
      case DataType::UInt64: f(std::uint64_t()); break;
      case DataType::Float16: f(Half()); break;
      case DataType::Float32: f(float()); break;
      case DataType::Float64: f(double()); break;
      case DataType::Bit: f(std::uint8_t()); break;
    }
  }

  /** An array entry, named as trx-python names it: <folder>/<name>.<dtype> or
   *  <folder>/<name>.<components>.<dtype>; trx-python splits at every dot, so names hold none. */
  struct ArrayName
  {
    std::string folder;
    std::string name;
    std::uint64_t components = 1;
    DataType type = DataType::Float32;
  };

  std::optional<ArrayName> ParseArrayName(const std::string& entryName)
  {
    ArrayName result;
    const auto slash = entryName.rfind('/');
    if (slash != std::string::npos)
      result.folder = entryName.substr(0, slash);
    const std::string base = slash == std::string::npos ? entryName : entryName.substr(slash + 1);

    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true)
    {
      const auto dot = base.find('.', start);
      parts.push_back(base.substr(start, dot == std::string::npos ? std::string::npos : dot - start));
      if (dot == std::string::npos)
        break;
      start = dot + 1;
    }
    if ((parts.size() != 2 && parts.size() != 3) || parts.front().empty())
      return std::nullopt;
    result.name = parts.front();
    if (parts.size() == 3)
    {
      const std::string& count = parts[1];
      if (count.empty() || count.size() > 9 || !std::all_of(count.begin(), count.end(), [](char c) { return c >= '0' && c <= '9'; }))
        return std::nullopt;
      result.components = std::stoull(count);
      if (result.components == 0)
        return std::nullopt;
    }
    const auto type = std::find_if(DataTypeNames.begin(), DataTypeNames.end(), [&parts](const DataTypeName& t) { return parts.back() == t.name; });
    if (type == DataTypeNames.end())
      return std::nullopt;
    result.type = type->type;
    return result;
  }

  std::string FolderOf(const std::string& entryName)
  {
    const auto slash = entryName.rfind('/');
    return slash == std::string::npos ? std::string() : entryName.substr(0, slash);
  }

  std::string BaseNameOf(const std::string& entryName)
  {
    const auto slash = entryName.rfind('/');
    return slash == std::string::npos ? entryName : entryName.substr(slash + 1);
  }

  /** The bytes of rows x components values of size bytes; throws when they exceed 64 bits. */
  std::uint64_t ArrayBytes(std::uint64_t rows, std::uint64_t components, std::size_t size, const std::string& path, const std::string& entryName)
  {
    const std::uint64_t rowBytes = components * size;
    if (rows > std::numeric_limits<std::uint64_t>::max() / rowBytes)
      mitkThrow() << path << ": " << entryName << " would hold more than 2^64 bytes.";
    return rows * rowBytes;
  }

  struct SourceEntry
  {
    /** path relative to the archive or folder, with '/' */
    std::string name;
    /** uncompressed size in bytes */
    std::uint64_t size = 0;
    // ZIP entries only
    std::streamoff dataStart = 0;
    int compressionMethod = Poco::Zip::ZipCommon::CM_STORE;
    bool encrypted = false;
    const Poco::Zip::ZipLocalFileHeader* zipHeader = nullptr;
  };

  /** The entries of a TRX ZIP archive or folder. */
  class Source
  {
  public:
    explicit Source(const std::string& path)
      : m_Path(path)
    {
      std::error_code error;
      if (std::filesystem::is_directory(path, error))
      {
        m_IsFolder = true;
        try
        {
          const std::filesystem::path root(path);
          for (const auto& item : std::filesystem::recursive_directory_iterator(root))
            if (item.is_regular_file())
              m_Entries.push_back({item.path().lexically_relative(root).generic_string(), item.file_size()});
        }
        catch (const std::filesystem::filesystem_error& e)
        {
          mitkThrow() << "Could not list the TRX folder " << path << ": " << e.what();
        }
      }
      else if (std::filesystem::is_regular_file(path, error))
      {
        std::ifstream in(path, std::ios::binary);
        if (!in)
          mitkThrow() << "Could not open " << path;
        try
        {
          m_Archive = std::make_unique<Poco::Zip::ZipArchive>(in);
        }
        catch (const Poco::Exception& e)
        {
          mitkThrow() << path << " is no TRX file: no valid ZIP archive (" << e.displayText() << ").";
        }
        for (auto it = m_Archive->headerBegin(); it != m_Archive->headerEnd(); ++it)
        {
          const Poco::Zip::ZipLocalFileHeader& header = it->second;
          if (header.isDirectory())
            continue;
          SourceEntry entry;
          entry.name = header.getFileName();
          entry.size = header.getUncompressedSize();
          entry.dataStart = header.getDataStartPos();
          entry.compressionMethod = header.getCompressionMethod();
          entry.encrypted = header.isEncrypted();
          entry.zipHeader = &header;
          m_Entries.push_back(entry);
        }
      }
      else
        mitkThrow() << path << " does not exist or is no file or folder.";
      std::sort(m_Entries.begin(), m_Entries.end(), [](const SourceEntry& a, const SourceEntry& b) { return a.name < b.name; });
    }

    const std::vector<SourceEntry>& GetEntries() const { return m_Entries; }

    const SourceEntry* Find(const std::string& name) const
    {
      for (const auto& entry : m_Entries)
        if (entry.name == name)
          return &entry;
      return nullptr;
    }

    /** Passes the bytes of the entry to consume in pieces of 1 MiB (the last one shorter), a
     *  multiple of every value size. Throws when the entry holds fewer bytes than its size. */
    void Read(const SourceEntry& entry, const std::function<void(const char*, std::size_t)>& consume) const
    {
      std::vector<char> buffer(static_cast<std::size_t>(std::min<std::uint64_t>(std::uint64_t(1) << 20, std::max<std::uint64_t>(entry.size, 1))));
      auto readAll = [&](std::istream& in)
      {
        std::uint64_t remaining = entry.size;
        while (remaining > 0)
        {
          const std::size_t n = static_cast<std::size_t>(std::min<std::uint64_t>(buffer.size(), remaining));
          in.read(buffer.data(), static_cast<std::streamsize>(n));
          if (static_cast<std::size_t>(in.gcount()) != n)
            mitkThrow() << m_Path << ": " << entry.name << " is truncated.";
          consume(buffer.data(), n);
          remaining -= n;
        }
      };

      if (m_IsFolder)
      {
        std::ifstream in(std::filesystem::path(m_Path) / entry.name, std::ios::binary);
        if (!in)
          mitkThrow() << "Could not open " << entry.name << " in " << m_Path;
        readAll(in);
        return;
      }

      if (entry.encrypted)
        mitkThrow() << m_Path << ": " << entry.name << " is encrypted.";
      std::ifstream in(m_Path, std::ios::binary);
      if (!in)
        mitkThrow() << "Could not open " << m_Path;
      if (entry.compressionMethod == Poco::Zip::ZipCommon::CM_STORE)
      {
        // Stored entries are read straight from the file, without a stream per value.
        in.seekg(entry.dataStart);
        if (!in)
          mitkThrow() << m_Path << ": " << entry.name << " lies outside the file.";
        readAll(in);
      }
      else if (entry.compressionMethod == Poco::Zip::ZipCommon::CM_DEFLATE)
      {
        try
        {
          Poco::Zip::ZipInputStream zipIn(in, *entry.zipHeader, true);
          // Poco's errors (CRC, inflate) reach the catch below instead of a short read.
          zipIn.exceptions(std::ios::badbit);
          readAll(zipIn);
        }
        catch (const Poco::Exception& e)
        {
          mitkThrow() << m_Path << ": could not decompress " << entry.name << " (" << e.displayText() << ").";
        }
      }
      else
        mitkThrow() << m_Path << ": " << entry.name << " has the unsupported ZIP compression method " << entry.compressionMethod
                    << " (TRX files are stored or deflated).";
    }

  private:
    std::string m_Path;
    bool m_IsFolder = false;
    std::unique_ptr<Poco::Zip::ZipArchive> m_Archive;
    std::vector<SourceEntry> m_Entries;
  };

  /** Calls store(index, value) for each value of the entry, decoded from the data type (halves as
   *  float); the entry holds a whole number of values. */
  template <typename F>
  void ForEachValue(const Source& source, const SourceEntry& entry, DataType type, F&& store)
  {
    WithType(type, [&](auto tag)
    {
      using S = decltype(tag);
      std::uint64_t index = 0;
      source.Read(entry, [&](const char* data, std::size_t bytes)
      {
        for (std::size_t k = 0; k < bytes; k += sizeof(S))
          store(index++, Decode<S>(data + k));
      });
    });
  }

  std::uint64_t Count(const nlohmann::json& header, const char* key, const std::string& path)
  {
    const auto value = header.find(key);
    if (value == header.end())
      mitkThrow() << path << ": header.json has no " << key << ".";
    if (!value->is_number_unsigned())
      mitkThrow() << path << ": " << key << " in header.json is no count (" << value->dump() << ").";
    return value->get<std::uint64_t>();
  }

  /** VOXEL_TO_RASMM (4 rows of 4 numbers, or 16 numbers, as trx-python reshapes them) and
   *  DIMENSIONS (3 numbers); an empty string, or what is missing. */
  std::string ReadGrid(const nlohmann::json& header, double voxelToRasmm[3][4], double dimensions[3])
  {
    const auto matrix = header.find("VOXEL_TO_RASMM");
    const auto dims = header.find("DIMENSIONS");
    if (matrix == header.end() || dims == header.end())
      return "header.json has no VOXEL_TO_RASMM or DIMENSIONS";

    std::vector<nlohmann::json> values;
    if (matrix->is_array() && matrix->size() == 4 && std::all_of(matrix->begin(), matrix->end(), [](const nlohmann::json& row) { return row.is_array() && row.size() == 4; }))
    {
      for (const auto& row : *matrix)
        for (const auto& value : row)
          values.push_back(value);
    }
    else if (matrix->is_array() && matrix->size() == 16)
      values.assign(matrix->begin(), matrix->end());
    if (values.size() != 16 || !std::all_of(values.begin(), values.end(), [](const nlohmann::json& v) { return v.is_number() && std::isfinite(v.get<double>()); }))
      return "VOXEL_TO_RASMM is no 4 x 4 matrix";
    if (!dims->is_array() || dims->size() != 3 || !std::all_of(dims->begin(), dims->end(), [](const nlohmann::json& v) { return v.is_number(); }))
      return "DIMENSIONS are no 3 numbers";

    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 4; ++c)
        voxelToRasmm[r][c] = values[static_cast<std::size_t>(4 * r + c)].get<double>();
    for (int i = 0; i < 3; ++i)
      dimensions[i] = (*dims)[static_cast<std::size_t>(i)].get<double>();
    return {};
  }

  unsigned char ToColorValue(double value)
  {
    if (!(value > 0))
      return 0;
    if (value >= 255)
      return 255;
    return static_cast<unsigned char>(std::lround(value));
  }

  void CheckRows(const SourceEntry& entry, const ArrayName& name, std::uint64_t rows, const char* rowName, const std::string& path)
  {
    const std::size_t size = SizeOf(name.type);
    if (entry.size != ArrayBytes(rows, name.components, size, path, entry.name))
      mitkThrow() << path << ": " << entry.name << " holds " << entry.size << " bytes, but " << rowName << " = " << rows << " rows of "
                  << name.components << " values of " << size << " bytes are " << rows * name.components * size << ".";
  }

  mitk::TrxFile::Contents ReadContents(const std::string& path, bool loadGroups)
  {
    const Source source(path);
    mitk::TrxFile::Contents contents;

    const SourceEntry* headerEntry = source.Find("header.json");
    if (headerEntry == nullptr)
      mitkThrow() << path << " is no TRX file: header.json is missing.";
    if (headerEntry->size > (std::uint64_t(1) << 24))
      mitkThrow() << path << ": header.json holds " << headerEntry->size << " bytes, too many for a TRX header.";
    std::string text;
    source.Read(*headerEntry, [&text](const char* data, std::size_t bytes) { text.append(data, bytes); });
    nlohmann::json header;
    try
    {
      header = nlohmann::json::parse(text);
    }
    catch (const nlohmann::json::exception& e)
    {
      mitkThrow() << path << ": header.json is no valid JSON (" << e.what() << ").";
    }
    if (!header.is_object())
      mitkThrow() << path << ": header.json holds no JSON object.";
    contents.numStreamlines = Count(header, "NB_STREAMLINES", path);
    contents.numVertices = Count(header, "NB_VERTICES", path);
    const std::string gridProblem = ReadGrid(header, contents.voxelToRasmm, contents.dimensions);
    contents.hasGrid = gridProblem.empty();
    if (!contents.hasGrid)
      MITK_INFO << path << ": no reference geometry (" << gridProblem << ").";

    // Entries by kind
    const SourceEntry* positionsEntry = nullptr;
    const SourceEntry* offsetsEntry = nullptr;
    ArrayName positionsName;
    ArrayName offsetsName;
    std::vector<std::pair<const SourceEntry*, ArrayName>> dps;
    std::vector<std::pair<const SourceEntry*, ArrayName>> dpv;
    std::vector<std::pair<const SourceEntry*, ArrayName>> groups;
    for (const auto& entry : source.GetEntries())
    {
      if (entry.name == "header.json")
        continue;
      const auto name = ParseArrayName(entry.name);
      const std::string folder = FolderOf(entry.name);
      const std::string base = BaseNameOf(entry.name);
      const bool isPositions = folder.empty() && base.starts_with("positions.");
      const bool isOffsets = folder.empty() && base.starts_with("offsets.");
      if (isPositions || isOffsets)
      {
        if (!name)
          mitkThrow() << path << ": " << entry.name << " has an unsupported name or data type.";
        const SourceEntry*& slot = isPositions ? positionsEntry : offsetsEntry;
        if (slot != nullptr)
          mitkThrow() << path << " holds two " << (isPositions ? "positions" : "offsets") << " arrays: " << slot->name << " and " << entry.name << ".";
        slot = &entry;
        (isPositions ? positionsName : offsetsName) = *name;
      }
      else if (name && name->folder == "dps")
        dps.emplace_back(&entry, *name);
      else if (name && name->folder == "dpv")
        dpv.emplace_back(&entry, *name);
      else if (name && name->folder == "groups")
        groups.emplace_back(&entry, *name);
      else
        contents.droppedEntries.push_back(entry.name);
    }

    // Offsets: NB_STREAMLINES + 1 values from 0 to NB_VERTICES, as trx-python requires them.
    if (offsetsEntry == nullptr)
    {
      if (contents.numStreamlines > 0)
        mitkThrow() << path << " is no TRX file: the offsets array is missing.";
      contents.offsets.assign(1, 0);
    }
    else
    {
      if (offsetsName.components != 1 || (offsetsName.type != DataType::UInt32 && offsetsName.type != DataType::UInt64))
        mitkThrow() << path << ": " << offsetsEntry->name << ": offsets must be uint32 or uint64 values.";
      if (contents.numStreamlines == std::numeric_limits<std::uint64_t>::max())
        mitkThrow() << path << ": NB_STREAMLINES is too large.";
      const std::uint64_t count = contents.numStreamlines + 1;
      const std::size_t size = SizeOf(offsetsName.type);
      if (offsetsEntry->size != ArrayBytes(count, 1, size, path, offsetsEntry->name))
        mitkThrow() << path << ": " << offsetsEntry->name << " holds " << offsetsEntry->size / size << " offsets, but NB_STREAMLINES + 1 = " << count << ".";
      contents.offsets.resize(static_cast<std::size_t>(count));
      ForEachValue(source, *offsetsEntry, offsetsName.type, [&contents](std::uint64_t i, auto value) { contents.offsets[i] = static_cast<std::uint64_t>(value); });
    }
    if (contents.offsets.front() != 0)
      mitkThrow() << path << ": the first offset is " << contents.offsets.front() << ", not 0.";
    for (std::size_t i = 1; i < contents.offsets.size(); ++i)
      if (contents.offsets[i] < contents.offsets[i - 1])
        mitkThrow() << path << ": the offsets decrease at streamline " << i << " (" << contents.offsets[i - 1] << ", then " << contents.offsets[i] << ").";
    if (contents.offsets.back() != contents.numVertices)
      mitkThrow() << path << ": the last offset is " << contents.offsets.back() << ", not NB_VERTICES = " << contents.numVertices << ".";

    // Positions: RAS millimetres of the file to MITK's LPS, in one pass.
    contents.positions = vtkSmartPointer<vtkFloatArray>::New();
    contents.positions->SetNumberOfComponents(3);
    if (positionsEntry == nullptr)
    {
      if (contents.numVertices > 0)
        mitkThrow() << path << " is no TRX file: the positions array is missing.";
    }
    else
    {
      if (positionsName.components != 3 || !IsFloat(positionsName.type))
        mitkThrow() << path << ": " << positionsEntry->name << ": positions must be float16, float32 or float64 triples.";
      CheckRows(*positionsEntry, positionsName, contents.numVertices, "NB_VERTICES", path);
      contents.positions->SetNumberOfTuples(static_cast<vtkIdType>(contents.numVertices));
      if (contents.numVertices > 0)
      {
        float* out = contents.positions->GetPointer(0);
        ForEachValue(source, *positionsEntry, positionsName.type, [&](std::uint64_t i, auto value)
        {
          const float v = static_cast<float>(value);
          if (!std::isfinite(v))
            mitkThrow() << path << ": point " << i / 3 << " of " << positionsEntry->name << " is not finite.";
          out[i] = i % 3 == 2 ? v : -v;
        });
      }
    }

    for (const auto& [entry, name] : dps)
    {
      CheckRows(*entry, name, contents.numStreamlines, "NB_STREAMLINES", path);
      if (name.name == "weights" && name.components == 1 && contents.weights == nullptr)
      {
        contents.weights = vtkSmartPointer<vtkFloatArray>::New();
        contents.weights->SetNumberOfValues(static_cast<vtkIdType>(contents.numStreamlines));
        float* out = contents.weights->GetPointer(0);
        ForEachValue(source, *entry, name.type, [out](std::uint64_t i, auto value) { out[i] = static_cast<float>(value); });
      }
      else
        contents.droppedEntries.push_back(entry->name);
    }

    // Colours: the "color" keyword of scilpy and MI-Brain (3 components), or trx-python's color_x,
    // color_y and color_z of its test files.
    const SourceEntry* colorEntry = nullptr;
    std::array<const SourceEntry*, 3> colorComponents = {nullptr, nullptr, nullptr};
    std::map<const SourceEntry*, DataType> types;
    for (const auto& [entry, name] : dpv)
    {
      CheckRows(*entry, name, contents.numVertices, "NB_VERTICES", path);
      types[entry] = name.type;
      if (name.name == "color" && name.components == 3 && colorEntry == nullptr)
        colorEntry = entry;
      else if (name.components == 1 && (name.name == "color_x" || name.name == "color_y" || name.name == "color_z") && colorComponents[static_cast<std::size_t>(name.name.back() - 'x')] == nullptr)
        colorComponents[static_cast<std::size_t>(name.name.back() - 'x')] = entry;
      else
        contents.droppedEntries.push_back(entry->name);
    }
    const bool allComponents = std::all_of(colorComponents.begin(), colorComponents.end(), [](const SourceEntry* e) { return e != nullptr; });
    if (colorEntry != nullptr || allComponents)
    {
      contents.colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
      contents.colors->SetNumberOfComponents(3);
      contents.colors->SetNumberOfTuples(static_cast<vtkIdType>(contents.numVertices));
    }
    if (colorEntry != nullptr && contents.numVertices > 0)
    {
      unsigned char* out = contents.colors->GetPointer(0);
      ForEachValue(source, *colorEntry, types[colorEntry], [out](std::uint64_t i, auto value) { out[i] = ToColorValue(static_cast<double>(value)); });
    }
    else if (allComponents && contents.numVertices > 0)
    {
      unsigned char* out = contents.colors->GetPointer(0);
      for (std::size_t c = 0; c < 3; ++c)
        ForEachValue(source, *colorComponents[c], types[colorComponents[c]], [out, c](std::uint64_t i, auto value) { out[3 * i + c] = ToColorValue(static_cast<double>(value)); });
    }
    for (const SourceEntry* entry : colorComponents)
      if (entry != nullptr && (colorEntry != nullptr || !allComponents))
        contents.droppedEntries.push_back(entry->name);

    for (const auto& [entry, name] : groups)
    {
      if (name.components != 1 || !IsIndexType(name.type))
      {
        contents.droppedEntries.push_back(entry->name);
        continue;
      }
      const std::size_t size = SizeOf(name.type);
      if (entry->size % size != 0)
        mitkThrow() << path << ": " << entry->name << " holds " << entry->size << " bytes, no whole number of " << size << "-byte indices.";
      mitk::TrxFile::Contents::Group group{name.name, {}};
      if (loadGroups)
      {
        group.streamlines.resize(static_cast<std::size_t>(entry->size / size));
        ForEachValue(source, *entry, name.type, [&](std::uint64_t i, auto value)
        {
          bool inRange = true;
          if constexpr (std::is_signed_v<decltype(value)>)
            inRange = value >= 0;
          if (!inRange || static_cast<std::uint64_t>(value) >= contents.numStreamlines)
            mitkThrow() << path << ": group " << name.name << " holds the streamline index " << std::to_string(value)
                        << ", not below NB_STREAMLINES = " << contents.numStreamlines << ".";
          group.streamlines[i] = static_cast<std::uint64_t>(value);
        });
      }
      contents.groups.push_back(std::move(group));
    }
    std::stable_sort(contents.groups.begin(), contents.groups.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
    return contents;
  }

  /** Stores value little-endian at out; the position after it. */
  template <typename T>
  char* StoreLittleEndian(T value, char* out)
  {
    std::memcpy(out, &value, sizeof(T));
    if constexpr (std::endian::native == std::endian::big)
      std::reverse(out, out + sizeof(T));
    return out + sizeof(T);
  }

  /** An input stream buffer of numRows rows of rowBytes bytes that fill produces chunk by chunk
   *  when they are read, so that an archive entry need not be held in memory. Poco's Compress
   *  reads the first byte, seeks to the end for the length and back, then reads to the end. */
  class GeneratedStreamBuf : public std::streambuf
  {
  public:
    /** Writes numRows rows from row firstRow on to out. */
    using Fill = std::function<void(std::uint64_t firstRow, std::size_t numRows, char* out)>;

    GeneratedStreamBuf(std::uint64_t numRows, std::size_t rowBytes, Fill fill)
      : m_Size(numRows * rowBytes)
      , m_RowBytes(rowBytes)
      , m_Fill(std::move(fill))
      , m_Buffer(rowBytes * std::max<std::size_t>(1, (std::size_t(1) << 20) / rowBytes))
    {
    }

  protected:
    int_type underflow() override
    {
      if (this->gptr() != nullptr && this->gptr() < this->egptr())
        return traits_type::to_int_type(*this->gptr());
      const std::uint64_t position = this->Position();
      if (position >= m_Size)
        return traits_type::eof();
      const std::uint64_t chunkStart = position / m_Buffer.size() * m_Buffer.size();
      const auto chunkBytes = static_cast<std::size_t>(std::min<std::uint64_t>(m_Buffer.size(), m_Size - chunkStart));
      m_Fill(chunkStart / m_RowBytes, chunkBytes / m_RowBytes, m_Buffer.data());
      m_ChunkStart = chunkStart;
      this->setg(m_Buffer.data(), m_Buffer.data() + (position - chunkStart), m_Buffer.data() + chunkBytes);
      return traits_type::to_int_type(*this->gptr());
    }

    pos_type seekoff(off_type offset, std::ios_base::seekdir direction, std::ios_base::openmode which) override
    {
      if ((which & std::ios_base::in) == 0)
        return pos_type(off_type(-1));
      const std::uint64_t base = direction == std::ios_base::beg ? 0 : direction == std::ios_base::cur ? this->Position() : m_Size;
      const auto target = static_cast<off_type>(base) + offset;
      if (target < 0 || static_cast<std::uint64_t>(target) > m_Size)
        return pos_type(off_type(-1));
      const auto position = static_cast<std::uint64_t>(target);
      if (this->eback() != nullptr && position >= m_ChunkStart && position <= m_ChunkStart + static_cast<std::uint64_t>(this->egptr() - this->eback()))
        this->setg(this->eback(), this->eback() + (position - m_ChunkStart), this->egptr());
      else
      {
        this->setg(nullptr, nullptr, nullptr);
        m_Position = position;
      }
      return pos_type(target);
    }

    pos_type seekpos(pos_type position, std::ios_base::openmode which) override
    {
      return this->seekoff(off_type(position), std::ios_base::beg, which);
    }

  private:
    std::uint64_t Position() const
    {
      return this->eback() != nullptr ? m_ChunkStart + static_cast<std::uint64_t>(this->gptr() - this->eback()) : m_Position;
    }

    const std::uint64_t m_Size;
    const std::size_t m_RowBytes;
    Fill m_Fill;
    std::vector<char> m_Buffer;
    std::uint64_t m_ChunkStart = 0;
    /** the position while there is no chunk */
    std::uint64_t m_Position = 0;
  };

  /** Walks the point ids of the fibers of a polydata in fiber order, the vertices of a TRX file;
   *  a visit before the current fiber starts over at the first one. */
  class VertexCursor
  {
  public:
    explicit VertexCursor(vtkPolyData* polyData)
      : m_PolyData(polyData)
      , m_NumFibers(polyData->GetNumberOfCells())
      , m_Buffer(vtkSmartPointer<vtkIdList>::New())
    {
      this->Load();
    }

    /** Calls f(pointId) for count vertices from vertex first on. */
    template <typename F>
    void Visit(std::uint64_t first, std::size_t count, F&& f)
    {
      if (first < m_FiberStart)
      {
        m_Fiber = 0;
        m_FiberStart = 0;
        this->Load();
      }
      for (std::uint64_t v = first; v < first + count; ++v)
      {
        while (v >= m_FiberStart + m_Ids.size())
        {
          m_FiberStart += m_Ids.size();
          ++m_Fiber;
          if (m_Fiber >= m_NumFibers)
            mitkThrow() << "Vertex " << v << " lies beyond the last fiber.";
          this->Load();
        }
        f(m_Ids[static_cast<std::size_t>(v - m_FiberStart)]);
      }
    }

  private:
    void Load()
    {
      m_Ids = m_Fiber < m_NumFibers ? mitk::FiberBundleDetail::FiberPointIds(m_PolyData, m_Fiber, m_Buffer) : std::span<const vtkIdType>();
    }

    vtkPolyData* m_PolyData;
    const vtkIdType m_NumFibers;
    vtkSmartPointer<vtkIdList> m_Buffer;
    vtkIdType m_Fiber = 0;
    std::uint64_t m_FiberStart = 0;
    std::span<const vtkIdType> m_Ids;
  };
}

float mitk::TrxFile::HalfToFloat(std::uint16_t half)
{
  const std::uint32_t sign = static_cast<std::uint32_t>(half & 0x8000u) << 16;
  const std::uint32_t exponent = (half >> 10) & 0x1fu;
  const std::uint32_t mantissa = half & 0x3ffu;
  if (exponent == 0)
  {
    // zero and subnormals: mantissa x 2^-24, exact in float
    const float value = std::ldexp(static_cast<float>(mantissa), -24);
    return sign != 0 ? -value : value;
  }
  if (exponent == 0x1f) // infinity and NaN, payload kept
    return std::bit_cast<float>(sign | 0x7f800000u | (mantissa << 13));
  return std::bit_cast<float>(sign | ((exponent + 112) << 23) | (mantissa << 13));
}

std::uint16_t mitk::TrxFile::FloatToHalf(float value)
{
  const std::uint32_t bits = std::bit_cast<std::uint32_t>(value);
  const auto sign = static_cast<std::uint16_t>((bits >> 16) & 0x8000u);
  const std::uint32_t magnitude = bits & 0x7fffffffu;
  if (magnitude > 0x7f800000u) // NaN: quiet, with the upper payload bits
    return static_cast<std::uint16_t>(sign | 0x7e00u | ((magnitude >> 13) & 0x3ffu));
  // 65520 lies halfway between the largest half, 65504 (odd mantissa), and 65536: to infinity
  if (magnitude >= 0x477ff000u)
    return static_cast<std::uint16_t>(sign | 0x7c00u);
  if (magnitude >= 0x38800000u) // normal halves from 2^-14 on: exponent rebiased by 127 - 15
  {
    std::uint32_t half = (magnitude - 0x38000000u) >> 13;
    const std::uint32_t rest = magnitude & 0x1fffu;
    if (rest > 0x1000u || (rest == 0x1000u && (half & 1u) != 0))
      ++half; // a carry into the exponent is the next half
    return static_cast<std::uint16_t>(sign | half);
  }
  // subnormal halves: value x 2^24 rounded to an integer; below 2^-25 (exponent 102) rounds to 0
  const std::uint32_t exponent = magnitude >> 23;
  if (exponent < 102)
    return sign;
  const std::uint32_t mantissa = (magnitude & 0x7fffffu) | 0x800000u;
  const std::uint32_t shift = 126 - exponent;
  std::uint32_t half = mantissa >> shift;
  const std::uint32_t rest = mantissa & ((1u << shift) - 1);
  const std::uint32_t halfway = 1u << (shift - 1);
  if (rest > halfway || (rest == halfway && (half & 1u) != 0))
    ++half;
  return static_cast<std::uint16_t>(sign | half);
}

mitk::TrxFile::Contents mitk::TrxFile::Read(const std::string& path, bool loadGroups)
{
  try
  {
    return ReadContents(path, loadGroups);
  }
  catch (const Poco::Exception& e)
  {
    mitkThrow() << path << ": " << e.displayText();
  }
}

void mitk::TrxFile::Write(const std::string& filename, const FiberBundle* fib, const WriteOptions& options)
{
  vtkPolyData* polyData = fib->GetFiberPolyData();
  const auto numFibers = static_cast<std::uint64_t>(polyData->GetNumberOfCells());
  std::vector<std::uint64_t> offsets(static_cast<std::size_t>(numFibers + 1), 0);
  vtkSmartPointer<vtkIdList> idBuffer = vtkSmartPointer<vtkIdList>::New();
  for (std::uint64_t i = 0; i < numFibers; ++i)
    offsets[i + 1] = offsets[i] + FiberBundleDetail::FiberPointIds(polyData, static_cast<vtkIdType>(i), idBuffer).size();
  const std::uint64_t numVertices = offsets.back();

  // The reference grid as MITK reads it; a grid read from a TRX file comes back unchanged.
  const TrackVisFile::Header grid = TrackVisFile::ConsistentHeader(fib->GetTrackVisHeader(), false);
  nlohmann::json matrix = nlohmann::json::array();
  for (int r = 0; r < 3; ++r)
    matrix.push_back({static_cast<double>(grid.vox_to_ras[r][0]), static_cast<double>(grid.vox_to_ras[r][1]),
                      static_cast<double>(grid.vox_to_ras[r][2]), static_cast<double>(grid.vox_to_ras[r][3])});
  matrix.push_back({0.0, 0.0, 0.0, 1.0});
  nlohmann::json header;
  header["DIMENSIONS"] = nlohmann::json::array({grid.dim[0], grid.dim[1], grid.dim[2]});
  header["VOXEL_TO_RASMM"] = matrix;
  header["NB_STREAMLINES"] = numFibers;
  header["NB_VERTICES"] = numVertices;
  std::istringstream headerStream(header.dump());

  const bool half = options.positionType == PositionType::Float16;
  const std::size_t positionBytes = half ? 2 : 4;
  const bool offsets64 = numVertices > std::numeric_limits<std::uint32_t>::max();
  vtkFloatArray* weights = options.weights ? fib->GetFiberWeights().GetPointer() : nullptr;
  if (weights != nullptr && static_cast<std::uint64_t>(weights->GetNumberOfValues()) != numFibers)
    mitkThrow() << "Cannot write " << filename << ": " << weights->GetNumberOfValues() << " fiber weights for " << numFibers << " fibers.";
  vtkUnsignedCharArray* colors = options.colors ? fib->GetFiberColors().GetPointer() : nullptr;
  if (colors != nullptr && (colors->GetNumberOfComponents() != 4 || colors->GetNumberOfTuples() < polyData->GetNumberOfPoints()))
    mitkThrow() << "Cannot write " << filename << ": the fiber colours are no RGBA values per point.";
  const bool entries = numFibers > 0;
  const std::uint64_t bytes = (entries ? numVertices * (3 * positionBytes + (colors != nullptr ? 3 : 0)) + (numFibers + 1) * (offsets64 ? 8 : 4)
                                         + (weights != nullptr ? 4 * numFibers : 0) : 0) + (std::uint64_t(1) << 20);

  std::ofstream out(filename, std::ios::binary | std::ios::trunc);
  if (!out)
    mitkThrow() << "Could not create " << filename;
  try
  {
    Poco::Zip::Compress compress(out, true, bytes >= options.zip64Threshold);
    const auto method = options.compression == Compression::Deflated ? Poco::Zip::ZipCommon::CM_DEFLATE : Poco::Zip::ZipCommon::CM_STORE;
    auto add = [&](const std::string& name, std::istream& in)
    {
      // Errors of the fill functions leave the stream instead of ending the entry early.
      in.exceptions(std::ios::badbit);
      compress.addFile(in, Poco::DateTime(), Poco::Path(name, Poco::Path::PATH_UNIX), method, Poco::Zip::ZipCommon::CL_NORMAL);
    };
    add("header.json", headerStream);

    // trx-python writes no positions and offsets for a file without streamlines.
    if (entries)
    {
      vtkPoints* points = polyData->GetPoints();
      VertexCursor positionCursor(polyData);
      GeneratedStreamBuf positionBuffer(numVertices, 3 * positionBytes, [&](std::uint64_t first, std::size_t count, char* data)
      {
        positionCursor.Visit(first, count, [&](vtkIdType id)
        {
          double p[3];
          points->GetPoint(id, p);
          // LPS to RAS; exact for MITK's float points, so that float32 files keep them bit for bit
          const float ras[3] = {static_cast<float>(-p[0]), static_cast<float>(-p[1]), static_cast<float>(p[2])};
          for (const float v : ras)
            data = half ? StoreLittleEndian(FloatToHalf(v), data) : StoreLittleEndian(v, data);
        });
      });
      std::istream positionStream(&positionBuffer);
      add(half ? "positions.3.float16" : "positions.3.float32", positionStream);

      GeneratedStreamBuf offsetBuffer(numFibers + 1, offsets64 ? 8 : 4, [&](std::uint64_t first, std::size_t count, char* data)
      {
        for (std::size_t k = 0; k < count; ++k)
          data = offsets64 ? StoreLittleEndian(offsets[first + k], data) : StoreLittleEndian(static_cast<std::uint32_t>(offsets[first + k]), data);
      });
      std::istream offsetStream(&offsetBuffer);
      add(offsets64 ? "offsets.uint64" : "offsets.uint32", offsetStream);

      if (colors != nullptr)
      {
        const unsigned char* rgba = colors->GetPointer(0);
        VertexCursor colorCursor(polyData);
        GeneratedStreamBuf colorBuffer(numVertices, 3, [&](std::uint64_t first, std::size_t count, char* data)
        {
          colorCursor.Visit(first, count, [&](vtkIdType id)
          {
            std::memcpy(data, rgba + 4 * id, 3);
            data += 3;
          });
        });
        std::istream colorStream(&colorBuffer);
        add("dpv/color.3.uint8", colorStream);
      }

      if (weights != nullptr)
      {
        GeneratedStreamBuf weightBuffer(numFibers, 4, [&](std::uint64_t first, std::size_t count, char* data)
        {
          for (std::size_t k = 0; k < count; ++k)
            data = StoreLittleEndian(weights->GetValue(static_cast<vtkIdType>(first + k)), data);
        });
        std::istream weightStream(&weightBuffer);
        add("dps/weights.float32", weightStream);
      }
    }
    compress.close();
  }
  catch (const Poco::Exception& e)
  {
    mitkThrow() << "Could not write " << filename << ": " << e.displayText();
  }
  catch (const std::ios_base::failure& e)
  {
    mitkThrow() << "Could not write " << filename << ": " << e.what();
  }
  // Buffered write errors (disk full) show up only here.
  out.close();
  if (!out)
    mitkThrow() << "Could not finish writing " << filename;
}
