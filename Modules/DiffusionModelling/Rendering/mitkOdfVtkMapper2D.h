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
    };

    /**
     * Whether the glyphs of a renderer have to be generated again: for new data, new node
     * properties or another time step (contentChanged), and for any change of the view.
     * \param generated the view the current glyphs were generated for; none if there are none
     */
    static bool IsGenerationRequired(const std::optional<GlyphView> &generated, const GlyphView &current, bool contentChanged);

    class LocalStorage : public Mapper::BaseLocalStorage
    {
    public:
      LocalStorage();

      vtkSmartPointer<vtkActor> m_Actor;
      vtkSmartPointer<vtkPolyDataMapper> m_Mapper;

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
