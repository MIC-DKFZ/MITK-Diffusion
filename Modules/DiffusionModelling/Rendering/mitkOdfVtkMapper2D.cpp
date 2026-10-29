/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkOdfVtkMapper2D.h"

#include <mitkBaseRenderer.h>
#include <mitkOdfNormalizationMethodProperty.h>
#include <mitkOdfScaleByProperty.h>
#include <mitkProperties.h>
#include <mitkTensorImage.h>

#include <vtkLookupTable.h>
#include <vtkProperty.h>

#include <algorithm>

mitk::OdfVtkMapper2D::LocalStorage::LocalStorage()
  : m_Actor(vtkSmartPointer<vtkActor>::New()), m_Mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
{
  m_Mapper->SetLookupTable(vtkSmartPointer<vtkLookupTable>::New());
  m_Mapper->ScalarVisibilityOn();
  m_Mapper->SetScalarModeToUsePointFieldData();
  m_Mapper->SelectColorArray("ODF_COLORS");
  m_Mapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  m_Actor->SetMapper(m_Mapper);
  m_Actor->GetProperty()->SetInterpolationToGouraud();
  m_Actor->VisibilityOff();
}

vtkProp *mitk::OdfVtkMapper2D::GetVtkProp(BaseRenderer *renderer)
{
  return m_LSH.GetLocalStorage(renderer)->m_Actor;
}

const char *mitk::OdfVtkMapper2D::GetVisibilityPropertyName(OdfGlyphGenerator::Orientation orientation)
{
  switch (orientation)
  {
    case OdfGlyphGenerator::Orientation::Sagittal:
      return "VisibleOdfs_S";
    case OdfGlyphGenerator::Orientation::Coronal:
      return "VisibleOdfs_C";
    default:
      return "VisibleOdfs_T";
  }
}

mitk::OdfGlyphGenerator::Settings mitk::OdfVtkMapper2D::GetSettings() const
{
  const DataNode *node = this->GetDataNode();
  OdfGlyphGenerator::Settings settings;
  node->GetFloatProperty("Scaling", settings.scaling);
  node->GetIntProperty("ShowMaxNumber", settings.maxGlyphs);
  if (auto normalization = dynamic_cast<const OdfNormalizationMethodProperty *>(node->GetProperty("Normalization")))
    settings.normalization = normalization->GetNormalization();
  if (auto scaleBy = dynamic_cast<const OdfScaleByProperty *>(node->GetProperty("ScaleBy")))
    settings.scaleBy = scaleBy->GetScaleBy();
  node->GetFloatProperty("IndexParam1", settings.indexParam1);
  node->GetFloatProperty("IndexParam2", settings.indexParam2);
  node->GetBoolProperty("DiffusionCore.Rendering.OdfVtkMapper.SwitchTensorView", settings.tensorEllipsoids);
  node->GetBoolProperty("DiffusionCore.Rendering.OdfVtkMapper.ColourisationModeBit", settings.colourByDirection);
  node->GetBoolProperty("DiffusionCore.Rendering.OdfVtkMapper.RandomModeBit", settings.randomPlacement);
  return settings;
}

void mitk::OdfVtkMapper2D::GenerateDataForRenderer(BaseRenderer *renderer)
{
  LocalStorage *localStorage = m_LSH.GetLocalStorage(renderer);
  const DataNode *node = this->GetDataNode();
  const auto *image = dynamic_cast<const Image *>(node->GetData());
  const PlaneGeometry *plane = renderer->GetCurrentWorldPlaneGeometry();
  if (image == nullptr || plane == nullptr)
  {
    this->ResetMapper(renderer);
    return;
  }

  const Vector2D viewportSize = renderer->GetViewportSizeInMM().GetVectorFromOrigin();
  const Vector2D viewportOrigin = renderer->GetOriginInMM().GetVectorFromOrigin();
  const bool viewChanged = localStorage->m_LastPlaneOrigin != plane->GetOrigin() ||
                           localStorage->m_LastPlaneAxis0 != plane->GetAxisVector(0) ||
                           localStorage->m_LastPlaneAxis1 != plane->GetAxisVector(1) ||
                           localStorage->m_LastViewportSize != viewportSize ||
                           localStorage->m_LastViewportOrigin != viewportOrigin;
  // GetPropertyList(renderer) creates the renderer's list on its first call, which counts as a
  // change; it is taken before the generate time, or the next render would regenerate once.
  const auto propertiesTime = std::max(node->GetPropertyList()->GetMTime(), node->GetPropertyList(renderer)->GetMTime());
  const auto lastTime = localStorage->GetLastGenerateDataTime().GetMTime();
  if (!viewChanged && !localStorage->IsGenerateDataRequired(renderer, this, this->GetDataNode()) && lastTime >= propertiesTime)
    return;

  localStorage->UpdateGenerateDataTime();
  localStorage->m_LastPlaneOrigin = plane->GetOrigin();
  localStorage->m_LastPlaneAxis0 = plane->GetAxisVector(0);
  localStorage->m_LastPlaneAxis1 = plane->GetAxisVector(1);
  localStorage->m_LastViewportSize = viewportSize;
  localStorage->m_LastViewportOrigin = viewportOrigin;

  bool visible = false;
  node->GetVisibility(visible, renderer, "visible");
  const auto slice = OdfGlyphGenerator::GetSlice(*image, this->GetTimestep(), *plane);
  // Glyphs are shown only where their orientation is switched on; a missing property hides them.
  bool orientationVisible = false;
  if (visible && slice)
    node->GetBoolProperty(GetVisibilityPropertyName(slice->orientation), orientationVisible, renderer);
  if (!orientationVisible)
  {
    this->ResetMapper(renderer);
    return;
  }

  // The rectangle of the plane in view.
  OdfGlyphGenerator::Window window;
  Point2D center2D;
  center2D[0] = viewportOrigin[0] + viewportSize[0] / 2;
  center2D[1] = viewportOrigin[1] + viewportSize[1] / 2;
  plane->Map(center2D, window.center);
  window.axis1 = plane->GetAxisVector(0);
  window.axis1.Normalize();
  window.axis2 = plane->GetAxisVector(1);
  window.axis2.Normalize();
  window.halfSize1 = viewportSize[0] / 2;
  window.halfSize2 = viewportSize[1] / 2;

  try
  {
    localStorage->m_Mapper->SetInputData(OdfGlyphGenerator::Generate(*image, this->GetTimestep(), *slice, &window, this->GetSettings()));
    localStorage->m_Actor->VisibilityOn();
  }
  catch (const Exception &e)
  {
    MITK_ERROR << e.GetDescription();
    this->ResetMapper(renderer);
  }
}

void mitk::OdfVtkMapper2D::ResetMapper(BaseRenderer *renderer)
{
  LocalStorage *localStorage = m_LSH.GetLocalStorage(renderer);
  localStorage->m_Mapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  localStorage->m_Actor->VisibilityOff();
}

void mitk::OdfVtkMapper2D::SetDefaultProperties(DataNode *node, BaseRenderer *renderer, bool overwrite)
{
  node->AddProperty("ShowMaxNumber", IntProperty::New(500), renderer, overwrite);
  node->AddProperty("Normalization", OdfNormalizationMethodProperty::New(), renderer, overwrite);

  auto scaleBy = OdfScaleByProperty::New();
  scaleBy->SetScaleByGFA();
  node->AddProperty("ScaleBy", scaleBy, renderer, overwrite);

  // Tensors are drawn larger and coloured by their principal direction.
  const bool isTensorImage = dynamic_cast<TensorImage *>(node->GetData()) != nullptr;
  node->AddProperty("DiffusionCore.Rendering.OdfVtkMapper.ColourisationModeBit", BoolProperty::New(isTensorImage), renderer, overwrite);
  node->AddProperty("Scaling", FloatProperty::New(isTensorImage ? 3.0 : 1.5), renderer, overwrite);

  node->AddProperty("IndexParam1", FloatProperty::New(2), renderer, overwrite);
  node->AddProperty("IndexParam2", FloatProperty::New(1), renderer, overwrite);
  node->AddProperty("visible", BoolProperty::New(true), renderer, overwrite);
  node->AddProperty("VisibleOdfs_T", BoolProperty::New(false), renderer, overwrite);
  node->AddProperty("VisibleOdfs_C", BoolProperty::New(false), renderer, overwrite);
  node->AddProperty("VisibleOdfs_S", BoolProperty::New(false), renderer, overwrite);
  node->AddProperty("layer", IntProperty::New(100), renderer, overwrite);
  node->AddProperty("DiffusionCore.Rendering.OdfVtkMapper.SwitchTensorView", BoolProperty::New(true), renderer, overwrite);
  node->AddProperty("DiffusionCore.Rendering.OdfVtkMapper.RandomModeBit", BoolProperty::New(true), renderer, overwrite);
}
