/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkOdfGlyphGenerator_h
#define mitkOdfGlyphGenerator_h

#include <MitkDiffusionModellingExports.h>

#include <mitkImage.h>
#include <mitkOdfNormalizationMethodProperty.h>
#include <mitkOdfScaleByProperty.h>
#include <mitkPlaneGeometry.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <optional>

namespace mitk
{
  /**
   * \brief Builds the glyphs of one slice of a tensor, ODF or SH image as one vtkPolyData,
   *        independent of any renderer.
   *
   * Each glyph is the ODF base mesh, its points scaled by the voxel's normalized ODF values and
   * centred on the voxel in world coordinates. The type of the image decides how its voxels are
   * turned into ODFs: TensorImage (6 values), OdfImage (ODF_SAMPLING_SIZE values) or ShImage
   * (coefficients in the image's SH convention). The output has the point arrays "ODF_COLORS"
   * (RGBA) and "Normals".
   *
   * Glyphs that are small on screen (Settings::pixelsPerMm) are geodesic spheres of fewer
   * directions, with the voxel's values at those directions; ODF images give each the value of the
   * nearest ODF direction. Their size, colour by direction and normalization still come from the
   * ODF at the ODF_SAMPLING_SIZE directions, and their values are limited to its range.
   */
  class MITKDIFFUSIONMODELLING_EXPORT OdfGlyphGenerator
  {
  public:
    /** View a slice belongs to; selects the node property "VisibleOdfs_T", "_S" or "_C". */
    enum class Orientation
    {
      Axial,
      Sagittal,
      Coronal
    };

    /** Display settings, read from the data node's properties by the mapper. */
    struct Settings
    {
      float scaling = 1;
      int normalization = ODFN_MAX;
      int scaleBy = ODFSB_NONE;
      float indexParam1 = 2;
      float indexParam2 = 1;
      /** Tensors as ellipsoids instead of tensor ODFs. */
      bool tensorEllipsoids = false;
      /** One colour per glyph from its principal direction instead of colours by ODF value. */
      bool colourByDirection = false;
      int maxGlyphs = 500;
      /**
       * With more voxels than maxGlyphs: a fixed pseudo-random subset of exactly maxGlyphs voxels,
       * in which a voxel keeps its glyph while it stays in view; else every n-th voxel.
       */
      bool randomPlacement = true;
      /**
       * Screen pixels per mm in the plane, for the mesh of the glyphs. A glyph is at most
       * scaling x glyph spacing x pixelsPerMm pixels wide; below 16 pixels it has 42 points (80
       * triangles), below 40 pixels 162 points (320 triangles), else the ODF_SAMPLING_SIZE points
       * of the ODF base mesh. 0 gives every glyph the base mesh.
       */
      double pixelsPerMm = 0;
    };

    /** The voxel layer of the image a plane shows. */
    struct Slice
    {
      /** Image axis perpendicular to the plane. */
      unsigned int normalAxis = 2;
      /** Voxel layer along normalAxis whose values are drawn. */
      itk::IndexValueType index = 0;
      /** Continuous index of the plane along normalAxis; the glyph centres lie on the plane. */
      double planeIndex = 0;
      Orientation orientation = Orientation::Axial;
      /** Smaller spacing of the two in-plane image axes; scales the glyphs. */
      double glyphSpacing = 1;
    };

    /** Rectangle in the plane, in world coordinates; voxels outside it get no glyph. */
    struct Window
    {
      Point3D center;
      /** Unit vectors along the two edges of the rectangle. */
      Vector3D axis1;
      Vector3D axis2;
      double halfSize1 = 0;
      double halfSize2 = 0;
    };

    /**
     * The slice of the image at time step t that the plane shows. None if the plane is not
     * parallel to two image axes, or lies more than half a voxel outside the image.
     */
    static std::optional<Slice> GetSlice(const Image &image, TimeStepType t, const PlaneGeometry &plane);

    /**
     * Glyphs of the slice at time step t, for the voxels inside the window (all voxels of the
     * slice without one), at most settings.maxGlyphs of them.
     * \throw mitk::Exception for images other than tensor, ODF and SH images of float values.
     */
    static vtkSmartPointer<vtkPolyData> Generate(const Image &image,
                                                 TimeStepType t,
                                                 const Slice &slice,
                                                 const Window *window,
                                                 const Settings &settings);

    /**
     * Whether the glyphs of the image with these settings can be drawn as instances of one sphere
     * (GenerateEllipsoidInstances): tensor images as ellipsoids, coloured by direction and
     * normalized by their maximum. Colours by value change along a glyph, and the other
     * normalizations change its shape, so those glyphs need a mesh of their own (Generate).
     */
    static bool CanDrawEllipsoidInstances(const Image &image, const Settings &settings);

    /**
     * The tensor ellipsoids of Generate as instances of a unit sphere (GetInstanceSphere), for the
     * same voxels: one point per glyph at its centre, with the point arrays "Orientation" (unit
     * quaternion w, x, y, z of the rotation that takes the sphere's x, y and z axes to the
     * eigenvectors in ascending order of their eigenvalues), "Scale" (the semi-axes along them)
     * and "ODF_COLORS" (RGBA, the colour of the principal direction, as Generate's). A semi-axis is
     * scaling x (GFA or principal curvature of the ellipsoid's ODF, or 1) x glyph spacing / 2 x
     * eigenvalue / largest eigenvalue, so the glyph is the exact ellipsoid as long as Generate's
     * largest radius; all are zero for a null tensor or an eigenvalue <= 0.
     * \throw mitk::Exception for images other than tensor images of float values.
     */
    static vtkSmartPointer<vtkPolyData> GenerateEllipsoidInstances(const Image &image,
                                                                   TimeStepType t,
                                                                   const Slice &slice,
                                                                   const Window *window,
                                                                   const Settings &settings);

    /**
     * The unit sphere of the instanced ellipsoids, with its normals: the mesh of Generate's glyphs
     * for this slice and settings (Settings::pixelsPerMm), so that small instances are as cheap to
     * draw as small glyph meshes. A new copy at each call; the meshes of different sizes have
     * different numbers of points.
     */
    static vtkSmartPointer<vtkPolyData> GetInstanceSphere(const Slice &slice, const Settings &settings);
  };
}

#endif
