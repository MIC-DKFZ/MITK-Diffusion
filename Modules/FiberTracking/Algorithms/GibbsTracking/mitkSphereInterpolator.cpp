/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkSphereInterpolator.h"

#include <usGetModuleContext.h>
#include <usModuleContext.h>
#include <usModule.h>
#include <usModuleResource.h>
#include <usModuleResourceStream.h>

#include <mitkOdfImage.h>

#include <fstream>
#include <exception>

static const std::string BaryCoordsFileName = "FiberTrackingLUTBaryCoords.bin";
static const std::string IndicesFileName = "FiberTrackingLUTIndices.bin";

SphereInterpolator::SphereInterpolator(const std::string& lutPath)
{
  // set before loading, so that no member is left uninitialized when loading fails
  size = 301;
  sN = (size-1)/2;
  nverts = ODF_SAMPLING_SIZE;
  beta = 0.5;

  inva = (sqrt(1+beta)-sqrt(beta));
  b = 1/(1-sqrt(1/beta + 1));

  m_ValidState = lutPath.length()==0 ? LoadLookuptables() : LoadLookuptables(lutPath);
}

SphereInterpolator::~SphereInterpolator()
{

}

bool SphereInterpolator::LoadLookuptables(const std::string& lutPath)
{
  MITK_INFO << "SphereInterpolator: loading lookuptables from custom path: " << lutPath;

  std::string path = lutPath; path.append(BaryCoordsFileName);
  std::ifstream BaryCoordsStream;
  BaryCoordsStream.open(path.c_str(), ios::in | ios::binary);
  if (!BaryCoordsStream.is_open())
  {
    MITK_INFO << "SphereInterpolator: could not load FiberTrackingLUTBaryCoords.bin from " << path;
    return false;
  }

  std::ifstream IndicesStream;
  path = lutPath; path.append("FiberTrackingLUTIndices.bin");
  IndicesStream.open(path.c_str(), ios::in | ios::binary);
  if (!IndicesStream.is_open())
  {
    MITK_INFO << "SphereInterpolator: could not load FiberTrackingLUTIndices.bin from " << path;
    return false;
  }

  if (LoadLookuptables(BaryCoordsStream, IndicesStream))
  {
    MITK_INFO << "SphereInterpolator: first and second lut loaded successfully";
    return true;
  }

  return false;
}

bool SphereInterpolator::LoadLookuptables()
{
  MITK_INFO << "SphereInterpolator: loading lookuptables";

  us::Module* module = us::GetModuleContext()->GetModule();
  us::ModuleResource BaryCoordsRes = module->GetResource(BaryCoordsFileName);
  if (!BaryCoordsRes.IsValid())
  {
    MITK_INFO << "Could not retrieve resource " << BaryCoordsFileName;
    return false;
  }

  us::ModuleResource IndicesRes = module->GetResource(IndicesFileName);
  if (!IndicesRes)
  {
    MITK_INFO << "Could not retrieve resource " << IndicesFileName;
    return false;
  }

  us::ModuleResourceStream BaryCoordsStream(BaryCoordsRes, std::ios_base::binary);
  us::ModuleResourceStream IndicesStream(IndicesRes, std::ios_base::binary);
  return LoadLookuptables(BaryCoordsStream, IndicesStream);
}

bool SphereInterpolator::LoadLookuptables(std::istream& BaryCoordsStream, std::istream& IndicesStream)
{
  // getInterpolation reads 3 values for each of 6 cube faces at size x size positions
  const std::size_t expected = static_cast<std::size_t>(3*6*size*size);

  barycoords.clear();
  float coord;
  BaryCoordsStream.seekg (0, ios::beg);
  while (BaryCoordsStream.read(reinterpret_cast<char *>(&coord), sizeof(coord)))
    barycoords.push_back(coord);

  indices.clear();
  int index;
  IndicesStream.seekg (0, ios::beg);
  while (IndicesStream.read(reinterpret_cast<char *>(&index), sizeof(index)))
    indices.push_back(index);

  if (barycoords.size()!=expected || indices.size()!=expected)
  {
    MITK_ERROR << "SphereInterpolator: lookup tables have " << barycoords.size() << " and " << indices.size()
               << " entries, expected " << expected << " each.";
    return false;
  }
  return true;
}
