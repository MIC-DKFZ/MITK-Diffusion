/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkCompositeMapper_h
#define mitkCompositeMapper_h

#include <MitkDiffusionModellingExports.h>

#include <mitkImageVtkMapper2D.h>
#include <mitkLocalStorageHandler.h>
#include <mitkOdfVtkMapper2D.h>
#include <mitkVtkMapper.h>

#include <vtkPropAssembly.h>
#include <vtkSmartPointer.h>

namespace mitk
{
  /**
   * \brief 2D mapper of tensor, ODF and SH images: the colour map of ImageVtkMapper2D with the
   *        glyphs of OdfVtkMapper2D on top.
   *
   * The class name and the image mapper in mapper slot 3 (set by the service activator) are a
   * contract with MITK: PlaneGeometryDataVtkMapper3D takes the image mapper from slot 3 of nodes
   * whose 2D mapper is named "CompositeMapper" to texture the slices in the 3D view.
   */
  class MITKDIFFUSIONMODELLING_EXPORT CompositeMapper : public VtkMapper
  {
  public:
    mitkClassMacro(CompositeMapper, VtkMapper);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

    int MitkRenderOverlay(BaseRenderer *renderer) override;
    int MitkRenderOpaqueGeometry(BaseRenderer *renderer) override;
    int MitkRenderTranslucentGeometry(BaseRenderer *renderer) override;
    int MitkRenderVolumetricGeometry(BaseRenderer *renderer) override;

    void SetDataNode(DataNode *node) override;
    void Update(BaseRenderer *renderer) override;

    ImageVtkMapper2D::Pointer GetImageMapper();

    vtkProp *GetVtkProp(BaseRenderer *renderer) override;
    bool HasVtkProp(const vtkProp *prop, BaseRenderer *renderer) override;
    void ReleaseGraphicsResources(BaseRenderer *renderer) override;
    bool IsLODEnabled(BaseRenderer *renderer) const override;

    static void SetDefaultProperties(DataNode *node, BaseRenderer *renderer = nullptr, bool overwrite = false);

    class LocalStorage : public Mapper::BaseLocalStorage
    {
    public:
      vtkSmartPointer<vtkPropAssembly> m_PropAssembly = vtkSmartPointer<vtkPropAssembly>::New();
    };

  protected:
    CompositeMapper();
    ~CompositeMapper() override = default;

  private:
    OdfVtkMapper2D::Pointer m_OdfMapper;
    ImageVtkMapper2D::Pointer m_ImgMapper;
    LocalStorageHandler<LocalStorage> m_LSH;
  };
}

#endif
