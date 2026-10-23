/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleTckReader.h"
#include "mitkFiberBundleMimeTypes.h"

#include <itksys/SystemTools.hxx>
#include <mitkExceptionMacro.h>
#include <mitkFiberBundle.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>


mitk::FiberBundleTckReader::FiberBundleTckReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_TCK_MIMETYPE_NAME(), "tck Fiber Bundle Reader (MRtrix format)" )
{
  this->RegisterService();
}

mitk::FiberBundleTckReader * mitk::FiberBundleTckReader::Clone() const
{
  return new FiberBundleTckReader(*this);
}

namespace
{
  /** Reads the MRtrix header up to the line "END": offset of the data ("file: . <offset>"),
   *  "datatype" and "count". Throws when the header is incomplete. */
  struct TckHeader
  {
    std::streamoff offset = -1;
    std::string datatype = "Float32LE";
    long long count = -1;
  };

  TckHeader ReadTckHeader(std::istream& in, const std::string& filename)
  {
    TckHeader header;
    std::string line;
    bool first = true;
    bool end = false;
    while (std::getline(in, line))
    {
      // MRtrix pads header lines with spaces (room for later updates of count etc.)
      const auto begin = line.find_first_not_of(" \t\r");
      line = begin==std::string::npos ? std::string() : line.substr(begin, line.find_last_not_of(" \t\r") - begin + 1);
      if (first && line!="mrtrix tracks")
        mitkThrow() << filename << " is no MRtrix tracks file (first line \"" << line << "\").";
      first = false;
      if (line=="END")
      {
        end = true;
        break;
      }

      const auto colon = line.find(':');
      if (colon==std::string::npos)
        continue;
      std::string key = line.substr(0, colon);
      key.erase(key.find_last_not_of(" \t") + 1);
      std::istringstream value(line.substr(colon+1));
      if (key=="file")
      {
        std::string dot;
        value >> dot >> header.offset;
        if (dot!="." || value.fail())
          mitkThrow() << "Unsupported \"file:\" entry in " << filename << ": \"" << line << "\"";
      }
      else if (key=="datatype")
        value >> header.datatype;
      else if (key=="count")
        value >> header.count;
    }
    if (!end)
      mitkThrow() << "MRtrix header of " << filename << " has no END line.";
    if (header.offset<0)
      mitkThrow() << "MRtrix header of " << filename << " has no \"file: . <offset>\" entry.";
    return header;
  }

  /** One coordinate triple in the file's data type and byte order, converted to double. */
  template <typename T>
  bool ReadTriple(std::istream& in, bool swapBytes, double* triple)
  {
    for (int c=0; c<3; ++c)
    {
      std::array<char, sizeof(T)> bytes;
      if (!in.read(bytes.data(), bytes.size()))
        return false;
      if (swapBytes)
        std::reverse(bytes.begin(), bytes.end());
      T value;
      std::memcpy(&value, bytes.data(), sizeof(T));
      triple[c] = static_cast<double>(value);
    }
    return true;
  }
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleTckReader::DoRead()
{
  std::vector<itk::SmartPointer<mitk::BaseData> > result;
  const std::string filename = this->GetInputLocation();
  MITK_INFO << "Loading tractogram (MRtrix format): " << itksys::SystemTools::GetFilenameName(filename);

  // Opened for reading only, so that read-only files and shares load.
  std::ifstream in(filename, std::ios::binary);
  if (!in)
    mitkThrow() << "Could not open " << filename;

  const TckHeader header = ReadTckHeader(in, filename);

  bool isDouble = false;
  bool isBigEndian = false;
  if (header.datatype=="Float32LE" || header.datatype=="Float32BE")
    isDouble = false;
  else if (header.datatype=="Float64LE" || header.datatype=="Float64BE")
    isDouble = true;
  else
    mitkThrow() << "Unsupported MRtrix data type " << header.datatype << " in " << filename;
  isBigEndian = header.datatype.compare(header.datatype.size()-2, 2, "BE")==0;
  const bool swapBytes = isBigEndian != (std::endian::native==std::endian::big);

  in.clear();
  in.seekg(header.offset);
  if (!in)
    mitkThrow() << "Data offset " << header.offset << " lies outside " << filename;

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

  // Tracks end with a NaN triple, the file with an Inf triple.
  bool terminated = false;
  double p[3];
  while (isDouble ? ReadTriple<double>(in, swapBytes, p) : ReadTriple<float>(in, swapBytes, p))
  {
    if (std::isinf(p[0]) || std::isinf(p[1]) || std::isinf(p[2]))
    {
      terminated = true;
      break;
    }
    else if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2]))
    {
      vtkNewCells->InsertNextCell(container);
      container = vtkSmartPointer<vtkPolyLine>::New();
    }
    else
    {
      // RAS (MRtrix) to LPS (MITK)
      p[0] = -p[0];
      p[1] = -p[1];
      vtkIdType id = vtkNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
  }

  // MRtrix itself reads truncated files up to the last complete track, so they load with a warning.
  if (!terminated)
    MITK_WARN << filename << " ends without the end-of-file marker; " << container->GetNumberOfPoints() << " points of an incomplete track are dropped.";
  if (header.count>=0 && header.count!=vtkNewCells->GetNumberOfCells())
    MITK_WARN << filename << " declares " << header.count << " tracks but holds " << vtkNewCells->GetNumberOfCells() << ".";

  vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  fiberPolyData->SetPoints(vtkNewPoints);
  fiberPolyData->SetLines(vtkNewCells);

  FiberBundle::Pointer fib = FiberBundle::New(fiberPolyData);
  result.push_back(fib.GetPointer());
  MITK_INFO << "Fiber bundle read";
  return result;
}
