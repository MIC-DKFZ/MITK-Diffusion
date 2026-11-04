/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundleTrxReader.h"
#include "mitkFiberBundleMimeTypes.h"
#include "mitkFiberPolyDataBuilder.h"
#include "mitkTrackVisFile.h"
#include "mitkTrxFile.h"

#include <itksys/SystemTools.hxx>
#include <mitkExceptionMacro.h>
#include <mitkFiberBundle.h>
#include <mitkStringProperty.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <algorithm>
#include <numeric>
#include <optional>
#include <sstream>

namespace
{
  const char* const GROUPS_OPTION = "Load groups as separate bundles";

  /** A bundle of the given streamlines of the file, in their order (all of them for nullptr), with
   *  their weights and colours. The polydata is built from the arrays as they are (points,
   *  offsets and point ids), not point by point. */
  mitk::FiberBundle::Pointer CreateBundle(const mitk::TrxFile::Contents& contents, const std::vector<std::uint64_t>* streamlines, const std::string& filename)
  {
    vtkSmartPointer<vtkFloatArray> positions = contents.positions;
    vtkSmartPointer<vtkUnsignedCharArray> colors = contents.colors;
    vtkSmartPointer<vtkFloatArray> weights = contents.weights;
    vtkSmartPointer<vtkIdTypeArray> offsets = vtkSmartPointer<vtkIdTypeArray>::New();
    if (streamlines == nullptr)
    {
      offsets->SetNumberOfValues(static_cast<vtkIdType>(contents.offsets.size()));
      for (std::size_t i = 0; i < contents.offsets.size(); ++i)
        offsets->SetValue(static_cast<vtkIdType>(i), static_cast<vtkIdType>(contents.offsets[i]));
    }
    else
    {
      std::uint64_t numPoints = 0;
      for (const auto s : *streamlines)
        numPoints += contents.offsets[s + 1] - contents.offsets[s];
      positions = vtkSmartPointer<vtkFloatArray>::New();
      positions->SetNumberOfComponents(3);
      positions->SetNumberOfTuples(static_cast<vtkIdType>(numPoints));
      if (contents.colors != nullptr)
      {
        colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        colors->SetNumberOfComponents(3);
        colors->SetNumberOfTuples(static_cast<vtkIdType>(numPoints));
      }
      if (contents.weights != nullptr)
      {
        weights = vtkSmartPointer<vtkFloatArray>::New();
        weights->SetNumberOfValues(static_cast<vtkIdType>(streamlines->size()));
      }
      offsets->SetNumberOfValues(static_cast<vtkIdType>(streamlines->size() + 1));
      offsets->SetValue(0, 0);
      vtkIdType point = 0;
      for (std::size_t i = 0; i < streamlines->size(); ++i)
      {
        const auto s = (*streamlines)[i];
        const auto first = static_cast<vtkIdType>(contents.offsets[s]);
        const auto count = static_cast<vtkIdType>(contents.offsets[s + 1]) - first;
        for (vtkIdType k = 0; k < 3 * count; ++k)
          positions->SetValue(3 * point + k, contents.positions->GetValue(3 * first + k));
        if (colors != nullptr)
          for (vtkIdType k = 0; k < 3 * count; ++k)
            colors->SetValue(3 * point + k, contents.colors->GetValue(3 * first + k));
        if (weights != nullptr)
          weights->SetValue(static_cast<vtkIdType>(i), contents.weights->GetValue(static_cast<vtkIdType>(s)));
        point += count;
        offsets->SetValue(static_cast<vtkIdType>(i + 1), point);
      }
    }

    vtkSmartPointer<vtkIdTypeArray> pointIds = vtkSmartPointer<vtkIdTypeArray>::New();
    pointIds->SetNumberOfValues(positions->GetNumberOfTuples());
    for (vtkIdType i = 0; i < positions->GetNumberOfTuples(); ++i)
      pointIds->SetValue(i, i);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->SetData(offsets, pointIds);
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetData(positions);

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    if (weights != nullptr)
    {
      weights->SetName("FIBER_WEIGHTS");
      polyData->GetCellData()->AddArray(weights);
    }
    if (colors != nullptr)
    {
      colors->SetName("FIBER_COLORS");
      polyData->GetPointData()->AddArray(colors);
    }

    // The constructor drops streamlines without points; weights and colours go along.
    mitk::FiberBundle::Pointer fib = mitk::FiberBundle::New(polyData);
    mitk::FiberBundleDetail::SetFiberArraysOfReader(fib, filename);
    return fib;
  }

  std::string Join(const std::vector<std::string>& names)
  {
    std::ostringstream text;
    for (std::size_t i = 0; i < names.size(); ++i)
      text << (i > 0 ? ", " : "") << names[i];
    return text.str();
  }
}

mitk::FiberBundleTrxReader::FiberBundleTrxReader()
  : mitk::AbstractFileReader(mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRX_MIMETYPE_NAME(), "TRX Fiber Bundle Reader")
{
  Options defaultOptions;
  defaultOptions[GROUPS_OPTION] = false;
  this->SetDefaultOptions(defaultOptions);

  this->RegisterService();
}

mitk::FiberBundleTrxReader * mitk::FiberBundleTrxReader::Clone() const
{
  return new FiberBundleTrxReader(*this);
}

mitk::IFileReader::ConfidenceLevel mitk::FiberBundleTrxReader::GetConfidenceLevel() const
{
  const std::string path = this->GetInputLocation();
  // The extension matched already; a folder without header.json is refused by DoRead, which
  // names what is missing.
  if (this->GetInputStream() == nullptr && itksys::SystemTools::FileIsDirectory(path))
    return Supported;
  return AbstractFileReader::GetConfidenceLevel();
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleTrxReader::DoRead()
{
  std::vector<itk::SmartPointer<mitk::BaseData> > result;
  const std::string filename = this->GetInputLocation();
  Options options = this->GetOptions();
  const bool loadGroups = us::any_cast<bool>(options[GROUPS_OPTION]);
  MITK_INFO << "Loading tractogram (TRX format): " << itksys::SystemTools::GetFilenameName(filename);

  const TrxFile::Contents contents = TrxFile::Read(filename, loadGroups);
  if (!contents.droppedEntries.empty())
    MITK_INFO << "Not loaded from " << filename << ": " << Join(contents.droppedEntries);

  // The reference grid keeps the axis order of the file, so that saving as TRX gives back
  // VOXEL_TO_RASMM and DIMENSIONS unchanged.
  std::optional<TrackVisFile::Header> header;
  if (contents.hasGrid)
  {
    try
    {
      header = TrackVisFile::HeaderOfGrid(contents.voxelToRasmm, contents.dimensions);
    }
    catch (const mitk::Exception& e)
    {
      MITK_WARN << filename << ": no reference geometry; " << e.GetDescription();
    }
  }

  // The whole bundle comes first: MitkFileFormatConverter and IOUtil::Load<FiberBundle> take
  // the first output.
  FiberBundle::Pointer fib = CreateBundle(contents, nullptr, filename);
  if (header)
    fib->SetTrackVisHeader(*header);
  result.push_back(fib.GetPointer());

  std::vector<std::string> groupNames;
  for (const auto& group : contents.groups)
    groupNames.push_back(group.name);
  if (!loadGroups && !groupNames.empty())
    MITK_INFO << groupNames.size() << " groups of " << filename << " (" << Join(groupNames) << ") not loaded; the reader option \""
              << GROUPS_OPTION << "\" loads them as bundles of their own.";
  if (loadGroups)
  {
    for (const auto& group : contents.groups)
    {
      FiberBundle::Pointer groupFib = CreateBundle(contents, &group.streamlines, filename);
      if (groupFib->GetNumFibers() == 0)
      {
        MITK_INFO << "Group " << group.name << " of " << filename << " holds no streamline with points; not loaded.";
        continue;
      }
      if (header)
        groupFib->SetTrackVisHeader(*header);
      groupFib->SetProperty("name", mitk::StringProperty::New(group.name));
      result.push_back(groupFib.GetPointer());
    }
  }

  MITK_INFO << "Fiber bundle read";
  return result;
}
