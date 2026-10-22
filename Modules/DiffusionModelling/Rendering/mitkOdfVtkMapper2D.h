/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkOdfVtkMapper2D_h
#define mitkOdfVtkMapper2D_h

#include <MitkDiffusionModellingExports.h>

#include <mitkLocalStorageHandler.h>
#include <mitkOdfGlyphGenerator.h>
#include <mitkVtkMapper.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

namespace mitk
{
  /**
   * \brief 2D mapper for the glyphs of tensor, ODF and SH images.
   *
   * Draws the glyphs of the image slice the renderer's plane shows, for the voxels in view, if
   * the node property "VisibleOdfs_T", "_S" or "_C" of the plane's orientation is on. The glyphs
   * themselves come from OdfGlyphGenerator; the node properties set in SetDefaultProperties
   * control them.
   */
  class MITKDIFFUSIONMODELLING_EXPORT OdfVtkMapper2D : public VtkMapper
  {
  public:
    mitkClassMacro(OdfVtkMapper2D, VtkMapper);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

    vtkProp *GetVtkProp(BaseRenderer *renderer) override;
    bool IsLODEnabled(BaseRenderer * /*renderer*/) const override { return true; }

    static void SetDefaultProperties(DataNode *node, BaseRenderer *renderer = nullptr, bool overwrite = false);

    /** Name of the node property that switches the glyphs of an orientation on. */
    static const char *GetVisibilityPropertyName(OdfGlyphGenerator::Orientation orientation);

    class LocalStorage : public Mapper::BaseLocalStorage
    {
    public:
      LocalStorage();

      vtkSmartPointer<vtkActor> m_Actor;
      vtkSmartPointer<vtkPolyDataMapper> m_Mapper;

      /** The view the glyphs were generated for. */
      Point3D m_LastPlaneOrigin;
      Vector3D m_LastPlaneAxis0;
      Vector3D m_LastPlaneAxis1;
      Vector2D m_LastViewportSize;
      Vector2D m_LastViewportOrigin;
    };

  protected:
    OdfVtkMapper2D() = default;
    ~OdfVtkMapper2D() override = default;

    void GenerateDataForRenderer(BaseRenderer *renderer) override;
    void ResetMapper(BaseRenderer *renderer) override;

  private:
    OdfGlyphGenerator::Settings GetSettings() const;

    LocalStorageHandler<LocalStorage> m_LSH;
  };
}

#endif
