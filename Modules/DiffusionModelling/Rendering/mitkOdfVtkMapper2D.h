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
#include <vtkGlyph3DMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

#include <optional>

namespace mitk
{
  /**
   * \brief 2D mapper for the glyphs of tensor, ODF and SH images.
   *
   * Draws the glyphs of the image slice the renderer's plane shows, for the voxels in view, if
   * the node property "VisibleOdfs_T", "_S" or "_C" of the plane's orientation is on. The glyphs
   * themselves come from OdfGlyphGenerator; the node properties set in SetDefaultProperties
   * control them.
   *
   * The glyphs are in world coordinates, so a pan or zoom leaves them in place: while the view
   * moves, the mapper keeps them, and fills the new view at the render of level of detail 1 that
   * MITK's RenderingManager runs about 200 ms after the last render of a renderer with LOD
   * enabled (IsLODEnabled). A new slice, new data or new properties are drawn at once, and so is
   * a zoom whose glyphs need a finer mesh, for the voxels the glyphs cover (GetViewToGenerate).
   *
   * Tensor ellipsoids coloured by direction (OdfGlyphGenerator::CanDrawEllipsoidInstances) are
   * drawn as GPU instances of one sphere, all other glyphs as one mesh.
   */
  class MITKDIFFUSIONMODELLING_EXPORT OdfVtkMapper2D : public VtkMapper
  {
  public:
    mitkClassMacro(OdfVtkMapper2D, VtkMapper);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

    vtkProp *GetVtkProp(BaseRenderer *renderer) override;

    /** The glyphs are already in world coordinates; the image geometry must not move them again. */
    void UpdateVtkTransform(BaseRenderer *) override {}

    /** Both VTK mappers hold GPU buffers, also the one the actor does not use at the moment. */
    void ReleaseGraphicsResources(BaseRenderer *renderer) override;

    /**
     * True while the renderer shows glyphs of the node, so that the RenderingManager renders it
     * once more at LOD 1 after an interaction, and only then.
     */
    bool IsLODEnabled(BaseRenderer *renderer) const override;

    static void SetDefaultProperties(DataNode *node, BaseRenderer *renderer = nullptr, bool overwrite = false);

    /** Name of the node property that switches the glyphs of an orientation on. */
    static const char *GetVisibilityPropertyName(OdfGlyphGenerator::Orientation orientation);

    /** What the glyphs of a renderer are generated for. */
    struct GlyphView
    {
      OdfGlyphGenerator::Slice slice;
      Vector3D planeAxis0;
      Vector3D planeAxis1;
      /** The part of the plane in view. */
      OdfGlyphGenerator::Window window;
      /** Zoom of the view, which chooses the glyph mesh (OdfGlyphGenerator::Settings). */
      double pixelsPerMm = 0;
      /** Points of each glyph's mesh at this zoom (OdfGlyphGenerator::GetPointsPerGlyph). */
      unsigned int pointsPerGlyph = 0;
    };

    /**
     * The view to generate the glyphs of a renderer for, or none if the current glyphs stay. The
     * current view for new data, new node properties or another time step (contentChanged), for
     * another slice or plane at any level of detail, and for a pan or zoom (another window or
     * pixelsPerMm) at LOD 1 and higher. A zoom whose glyphs need a finer mesh than the current
     * ones (more pointsPerGlyph) at LOD 0: the generated window at the current zoom, so the same
     * voxels get finer glyphs at once; a magnified screenshot renders at LOD 0 with its
     * magnification's zoom and so gets the mesh it shows.
     * \param generated the view the current glyphs were generated for; none if there are none
     * \param lod the level of detail of the render, RenderingManager::GetNextLOD
     */
    static std::optional<GlyphView> GetViewToGenerate(const std::optional<GlyphView> &generated, const GlyphView &current, bool contentChanged, int lod);

    class LocalStorage : public Mapper::BaseLocalStorage
    {
    public:
      LocalStorage();

      vtkSmartPointer<vtkActor> m_Actor;
      /** Draws the glyph mesh of OdfGlyphGenerator::Generate. */
      vtkSmartPointer<vtkPolyDataMapper> m_Mapper;
      /** Draws the tensor ellipsoids of OdfGlyphGenerator::GenerateEllipsoidInstances. */
      vtkSmartPointer<vtkGlyph3DMapper> m_InstanceMapper;

      /** The view the glyphs were generated (or failed) for; none while no glyphs are shown. */
      std::optional<GlyphView> m_GeneratedView;
    };

  protected:
    OdfVtkMapper2D() = default;
    ~OdfVtkMapper2D() override = default;

    void GenerateDataForRenderer(BaseRenderer *renderer) override;
    void ResetMapper(BaseRenderer *renderer) override;

  private:
    OdfGlyphGenerator::Settings GetSettings() const;

    /** The view to draw glyphs for; none if the renderer shows no glyphs of the node. */
    std::optional<GlyphView> GetVisibleView(BaseRenderer *renderer) const;

    LocalStorageHandler<LocalStorage> m_LSH;
  };
}

#endif
