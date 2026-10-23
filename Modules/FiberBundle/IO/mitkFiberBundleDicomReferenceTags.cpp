/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundleDicomReferenceTags.h"

#include <mitkExceptionMacro.h>

#include <array>
#include <utility>

std::vector<std::string> mitk::FiberBundleDicomTags::GetEditableTags()
{
  return {PatientId, PatientName, StudyInstanceUid, SeriesInstanceUid, SopInstanceUid, SopClassUid,
          FrameOfReferenceUid, AlgorithmFamilyCodeValue, AlgorithmFamilyCodeMeaning, ModelCodeValue,
          ModelCodeMeaning, AnatomyValue, AnatomyMeaning};
}

std::vector<std::string> mitk::CopyDicomReferenceTags(const IPropertyProvider* source, FiberBundle* target)
{
  if (source == nullptr || target == nullptr)
    mitkThrow() << "Copying DICOM reference tags needs a source and a fiber bundle.";

  // bundle property (and name used by the diffusion DICOM import), MITK's DICOM tag property
  using namespace FiberBundleDicomTags;
  const std::array<std::pair<std::string, std::string>, 7> tags = {{
    {PatientId, "DICOM.0010.0020"},
    {PatientName, "DICOM.0010.0010"},
    {StudyInstanceUid, "DICOM.0020.000D"},
    {SeriesInstanceUid, "DICOM.0020.000E"},
    {SopInstanceUid, "DICOM.0008.0018"},
    {SopClassUid, "DICOM.0008.0016"},
    {FrameOfReferenceUid, "DICOM.0020.0052"}}};

  std::vector<std::string> missing;
  for (const auto& [name, dicomTag] : tags)
  {
    auto property = source->GetConstProperty(name);
    if (property.IsNull())
      property = source->GetConstProperty(dicomTag);
    const std::string value = property.IsNotNull() ? property->GetValueAsString() : std::string();
    if (value.empty())
      missing.push_back(name);
    else
      target->GetPropertyList()->SetStringProperty(name.c_str(), value.c_str());
  }
  return missing;
}
