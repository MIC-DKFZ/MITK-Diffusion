/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkCompositeMapper.h"

#include <mitkLevelWindowProperty.h>

mitk::CompositeMapper::CompositeMapper()
  : m_OdfMapper(OdfVtkMapper2D::New()), m_ImgMapper(ImageVtkMapper2D::New())
{
}

int mitk::CompositeMapper::MitkRenderOverlay(BaseRenderer *renderer)
{
  return m_ImgMapper->MitkRenderOverlay(renderer) + m_OdfMapper->MitkRenderOverlay(renderer);
}

int mitk::CompositeMapper::MitkRenderOpaqueGeometry(BaseRenderer *renderer)
{
  return m_ImgMapper->MitkRenderOpaqueGeometry(renderer) + m_OdfMapper->MitkRenderOpaqueGeometry(renderer);
}

int mitk::CompositeMapper::MitkRenderTranslucentGeometry(BaseRenderer *renderer)
{
  return m_ImgMapper->MitkRenderTranslucentGeometry(renderer) + m_OdfMapper->MitkRenderTranslucentGeometry(renderer);
}

int mitk::CompositeMapper::MitkRenderVolumetricGeometry(BaseRenderer *renderer)
{
  return m_ImgMapper->MitkRenderVolumetricGeometry(renderer) + m_OdfMapper->MitkRenderVolumetricGeometry(renderer);
}

void mitk::CompositeMapper::SetDataNode(DataNode *node)
{
  Superclass::SetDataNode(node);
  m_ImgMapper->SetDataNode(node);
  m_OdfMapper->SetDataNode(node);
}

void mitk::CompositeMapper::Update(BaseRenderer *renderer)
{
  // Each sub-mapper checks itself whether its output is out of date.
  m_OdfMapper->Update(renderer);
  m_ImgMapper->Update(renderer);
}

mitk::ImageVtkMapper2D::Pointer mitk::CompositeMapper::GetImageMapper()
{
  return m_ImgMapper;
}

vtkProp *mitk::CompositeMapper::GetVtkProp(BaseRenderer *renderer)
{
  // The image mapper switches its prop between its actors and an empty assembly, so the parts are
  // checked on every call; the assembly itself stays the same per renderer.
  vtkPropAssembly *assembly = m_LSH.GetLocalStorage(renderer)->m_PropAssembly;
  vtkProp *odfProp = m_OdfMapper->GetVtkProp(renderer);
  vtkProp *imageProp = m_ImgMapper->GetVtkProp(renderer);
  vtkPropCollection *parts = assembly->GetParts();
  if (parts->GetNumberOfItems() != 2 || parts->GetItemAsObject(0) != odfProp || parts->GetItemAsObject(1) != imageProp)
  {
    while (parts->GetNumberOfItems() > 0)
      assembly->RemovePart(vtkProp::SafeDownCast(parts->GetItemAsObject(0)));
    assembly->AddPart(odfProp);
    assembly->AddPart(imageProp);
  }
  return assembly;
}

bool mitk::CompositeMapper::HasVtkProp(const vtkProp *prop, BaseRenderer *renderer)
{
  return prop == m_LSH.GetLocalStorage(renderer)->m_PropAssembly || m_OdfMapper->HasVtkProp(prop, renderer) ||
         m_ImgMapper->HasVtkProp(prop, renderer);
}

void mitk::CompositeMapper::ReleaseGraphicsResources(BaseRenderer *renderer)
{
  m_ImgMapper->ReleaseGraphicsResources(renderer);
  m_OdfMapper->ReleaseGraphicsResources(renderer);
}

bool mitk::CompositeMapper::IsLODEnabled(BaseRenderer *renderer) const
{
  return m_ImgMapper->IsLODEnabled(renderer) || m_OdfMapper->IsLODEnabled(renderer);
}

void mitk::CompositeMapper::SetDefaultProperties(DataNode *node, BaseRenderer *renderer, bool overwrite)
{
  OdfVtkMapper2D::SetDefaultProperties(node, renderer, overwrite);
  ImageVtkMapper2D::SetDefaultProperties(node, renderer, overwrite);

  LevelWindow opacityLevelWindow;
  opacityLevelWindow.SetRangeMinMax(0, 255);
  opacityLevelWindow.SetWindowBounds(0, 0);
  node->AddProperty("opaclevelwindow", LevelWindowProperty::New(opacityLevelWindow));
}
