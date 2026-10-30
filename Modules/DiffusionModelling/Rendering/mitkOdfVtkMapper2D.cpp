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
#include <mitkRenderingManager.h>
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

namespace
{
  bool SameWindow(const mitk::OdfGlyphGenerator::Window &a, const mitk::OdfGlyphGenerator::Window &b)
  {
    return a.center == b.center && a.axis1 == b.axis1 && a.axis2 == b.axis2 && a.halfSize1 == b.halfSize1 &&
           a.halfSize2 == b.halfSize2;
  }
}

bool mitk::OdfVtkMapper2D::IsGenerationRequired(const std::optional<GlyphView> &generated, const GlyphView &current, bool contentChanged, int lod)
{
  if (!generated || contentChanged)
    return true;
  const bool sliceChanged = generated->slice.normalAxis != current.slice.normalAxis ||
                            generated->slice.index != current.slice.index ||
                            generated->slice.planeIndex != current.slice.planeIndex ||
                            generated->planeAxis0 != current.planeAxis0 || generated->planeAxis1 != current.planeAxis1;
  if (sliceChanged)
    return true;
  // Generating and uploading thousands of glyphs at every mouse step of a pan or zoom made the
  // view lag; the glyphs in place are still right, only the area they cover is out of date.
  return lod > 0 && !SameWindow(generated->window, current.window);
}

bool mitk::OdfVtkMapper2D::IsLODEnabled(BaseRenderer *renderer) const
{
  return this->GetVisibleView(renderer).has_value();
}

std::optional<mitk::OdfVtkMapper2D::GlyphView> mitk::OdfVtkMapper2D::GetVisibleView(BaseRenderer *renderer) const
{
  const DataNode *node = this->GetDataNode();
  const auto *image = dynamic_cast<const Image *>(node->GetData());
  const PlaneGeometry *plane = renderer->GetCurrentWorldPlaneGeometry();
  if (image == nullptr || plane == nullptr)
    return std::nullopt;

  bool visible = false;
  node->GetVisibility(visible, renderer, "visible");
  if (!visible)
    return std::nullopt;
  const auto slice = OdfGlyphGenerator::GetSlice(*image, this->GetTimestep(), *plane);
  // Glyphs are shown only where their orientation is switched on; a missing property hides them.
  bool orientationVisible = false;
  if (slice)
    node->GetBoolProperty(GetVisibilityPropertyName(slice->orientation), orientationVisible, renderer);
  if (!orientationVisible)
    return std::nullopt;

  GlyphView view;
  view.slice = *slice;
  view.planeAxis0 = plane->GetAxisVector(0);
  view.planeAxis1 = plane->GetAxisVector(1);

  // The rectangle of the plane in view.
  const Vector2D viewportSize = renderer->GetViewportSizeInMM().GetVectorFromOrigin();
  const Vector2D viewportOrigin = renderer->GetOriginInMM().GetVectorFromOrigin();
  Point2D center2D;
  center2D[0] = viewportOrigin[0] + viewportSize[0] / 2;
  center2D[1] = viewportOrigin[1] + viewportSize[1] / 2;
  plane->Map(center2D, view.window.center);
  view.window.axis1 = view.planeAxis0;
  view.window.axis1.Normalize();
  view.window.axis2 = view.planeAxis1;
  view.window.axis2.Normalize();
  view.window.halfSize1 = viewportSize[0] / 2;
  view.window.halfSize2 = viewportSize[1] / 2;
  return view;
}

void mitk::OdfVtkMapper2D::GenerateDataForRenderer(BaseRenderer *renderer)
{
  LocalStorage *localStorage = m_LSH.GetLocalStorage(renderer);
  const DataNode *node = this->GetDataNode();
  // GetPropertyList(renderer) creates the renderer's list on its first call, which counts as a
  // change; it is taken before the generate time, or the next render would regenerate once.
  const auto propertiesTime = std::max(node->GetPropertyList()->GetMTime(), node->GetPropertyList(renderer)->GetMTime());
  const bool contentChanged = localStorage->IsGenerateDataRequired(renderer, this, this->GetDataNode()) ||
                              localStorage->GetLastGenerateDataTime().GetMTime() < propertiesTime;

  const auto view = this->GetVisibleView(renderer);
  if (!view)
  {
    if (localStorage->m_GeneratedView)
    {
      localStorage->m_GeneratedView.reset();
      this->ResetMapper(renderer);
    }
    return;
  }
  if (!IsGenerationRequired(localStorage->m_GeneratedView, *view, contentChanged, RenderingManager::GetInstance()->GetNextLOD(renderer)))
    return;

  localStorage->UpdateGenerateDataTime();
  // Also set when the generation fails, so that the error is not repeated at every render.
  localStorage->m_GeneratedView = view;
  try
  {
    const auto *image = static_cast<const Image *>(node->GetData());
    localStorage->m_Mapper->SetInputData(OdfGlyphGenerator::Generate(*image, this->GetTimestep(), view->slice, &view->window, this->GetSettings()));
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
