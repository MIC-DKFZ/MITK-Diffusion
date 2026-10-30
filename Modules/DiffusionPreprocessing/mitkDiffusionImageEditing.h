/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionImageEditing_h
#define mitkDiffusionImageEditing_h

#include <MitkDiffusionPreprocessingExports.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkImage.h>

#include <array>
#include <map>
#include <optional>
#include <vector>

namespace mitk
{
  /**
   * \brief Edits of images and DWIs offered by the preprocessing view.
   *
   * Each function returns a new image and leaves its input unchanged. For DWIs the result follows
   * the gradient convention of MITK Diffusion: the modified gradients are world (LPS) vectors.
   * Invalid input throws mitk::Exception (or itk::ExceptionObject from the ITK filters used).
   */
  namespace DiffusionImageEditing
  {
    using MatrixType = DiffusionPropertyHelper::MeasurementFrameType;

    /**
     * \brief Mirrors the voxel data along the selected voxel axes and keeps the geometry.
     *
     * The image content is mirrored along the world direction of each flipped axis (a column of
     * the direction matrix), and so are the gradients of a DWI. Supports DWIs and 3D images with
     * one component.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer FlipAxes(const Image* image, const std::array<bool, 3>& axes);

    /** \brief The header fields to change; fields left empty keep their exact current value. */
    struct HeaderChange
    {
      std::optional<Vector3D> spacing;
      std::optional<Point3D> origin;
      /** Columns: world direction of each image axis; they are normalized and must be orthogonal. */
      std::optional<MatrixType> direction;
      /** Replaces the stored measurement frame (frame axes as rows); DWIs only. */
      std::optional<MatrixType> measurementFrame;
      /**
       * Whether the gradients of a DWI turn with a new direction matrix, i.e. keep their relation to
       * the voxel grid (as FSL bvecs do). Otherwise they stay fixed in world.
       */
      bool rotateGradientsWithImage = true;
    };

    /**
     * \brief Changes spacing, origin, direction matrix and measurement frame.
     *
     * When the direction matrix or the frame of a DWI changes, the new modified gradients are
     * computed explicitly and stored as the measurement frame applied to the kept original
     * gradients, with the image rotation switched off; they then no longer depend on whether the
     * input applied the image rotation.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ApplyHeader(const Image* image, const HeaderChange& change);

    /** \brief Whether the modified gradients of the DWI include the image rotation (legacy NRRD, DICOM). */
    MITKDIFFUSIONPREPROCESSING_EXPORT bool AppliesImageRotationToGradients(const Image* dwi);

    /** \brief Recomputes the modified gradients of the DWI without the image rotation. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ClearGradientRotation(const Image* dwi);

    /**
     * \brief Concatenates the volumes of DWIs with the same geometry.
     *
     * The gradients are rescaled to the largest reference b-value and stored as world vectors with
     * an identity frame and without the image rotation.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer MergeDwis(const std::vector<const Image*>& dwis);

    /** \brief Gradients with rounded b-values and the reference b-value they refer to. */
    struct RoundedBValues
    {
      DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer gradients;
      double referenceBValue = 0;
    };

    /**
     * \brief Rounds the b-value of each gradient to a multiple of the rounding value.
     *
     * The new reference b-value is the smallest rounded b-value above 0; the gradients are rescaled
     * to it, and a gradient whose b-value rounds to 0 becomes the zero vector (a b0 volume). Throws
     * for a reference b-value of 0, a rounding value below 1 and gradients without a weighted
     * volume.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT RoundedBValues RoundBValues(const DiffusionPropertyHelper::GradientDirectionsContainerType* gradients, double referenceBValue, int roundingValue);

    /** \brief The DWI with rounded b-values (see above). */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer RoundBValues(const Image* dwi, int roundingValue);

    /**
     * \brief Removes the given numbers of voxels at the lower and upper end of each axis.
     *
     * The remaining voxels keep their world position. Supports DWIs and 3D images with one
     * component; throws if no voxel would remain along an axis.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer Crop(const Image* image, const std::array<unsigned int, 3>& lower, const std::array<unsigned int, 3>& upper);

    enum class ResampleGrid { SamplingFactor, Spacing, Size };
    enum class Interpolation { NearestNeighbour, Linear, BSpline, WindowedSinc };

    /**
     * \brief Resamples onto a grid given relative to the input (itk::ResampledImageGrid) that keeps
     * the outer corner of the image and its direction.
     *
     * Supports DWIs (volume by volume) and 3D images with one component.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer Resample(const Image* image, ResampleGrid grid, const Vector3D& value, Interpolation interpolation);

    /** \brief One volume of a DWI as a 3D image with the DWI's geometry and no properties. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ExtractVolume(const Image* dwi, unsigned int index);

    /** \brief The DWI without the given volumes; an index listed twice counts once. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer RemoveVolumes(const Image* dwi, const std::vector<unsigned int>& indices);

    /** \brief The average of the b0 volumes (see mitk::DiffusionPropertyHelper::IsB0) as a 3D image without properties. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ExtractB0(const Image* dwi);

    /** \brief The b0 volumes (see mitk::DiffusionPropertyHelper::IsB0) as the time steps of a 3D+t image without properties. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ExtractB0Volumes(const Image* dwi);

    /** \brief Negates the selected world (LPS) components of every gradient. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer FlipGradients(const Image* dwi, const std::array<bool, 3>& components);

    /**
     * \brief Reorders the world (LPS) components of every gradient: component i of the result is
     * component order[i] of the input. Throws if order is not a permutation of 0, 1 and 2.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer PermuteGradientComponents(const Image* dwi, const std::array<unsigned int, 3>& order);

    /** \brief Negates every gradient with a negative x (world) component. */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer MirrorGradientsToHalfSphere(const Image* dwi);

    /**
     * \brief Keeps per selected shell (b-value of the DWI's b-value map) the given number of
     * evenly spread directions (itk::ElectrostaticRepulsionDiffusionGradientReductionFilter); with
     * useFirstN the first ones. Volumes of shells not listed are dropped.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer ReduceGradients(const Image* dwi, const std::map<unsigned int, unsigned int>& directionsPerShell, bool useFirstN);

    /** \brief Averages volumes whose gradients agree within the precision (DiffusionPropertyHelper::AverageRedundantGradients). */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer AverageRepeatedGradients(const Image* dwi, double precision);
  }
}

#endif
