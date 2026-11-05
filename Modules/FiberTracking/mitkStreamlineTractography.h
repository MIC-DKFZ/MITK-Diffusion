/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkStreamlineTractography_h
#define mitkStreamlineTractography_h

#include <MitkFiberTrackingExports.h>

#include <itkStreamlineTrackingFilter.h>
#include <mitkImage.h>
#include <mitkStreamlineTractographyParameters.h>

#include <memory>
#include <vector>

namespace mitk
{
  class TrackingDataHandler;

  /**
   * \brief The tracking handler of a streamline tractography run, built from the model images
   *        as the Streamline Tractography view and MitkStreamlineTractography build it.
   *
   * \param modelImages One SH, ODF or peak image, or tensor images: several only in
   *        deterministic mode, as a multi-tensor model; in probabilistic mode a tensor image is
   *        tracked as its ODFs.
   * \param mode The tracking mode, which decides how tensors are tracked.
   * \param faImage Optional FA (tensors) or GFA (SH, ODF) image for the cutoff; without it the
   *        handler takes the values of its own model. Not used for peaks.
   * \return The handler, which the tracker only borrows
   *         (itk::StreamlineTrackingFilter::SetTrackingHandler).
   * \throw mitk::Exception without an image, for an image that is no tensor, ODF, SH or peak
   *        image, for images of several types, and for several images that are not tensors in
   *        deterministic mode.
   */
  MITKFIBERTRACKING_EXPORT std::unique_ptr<TrackingDataHandler> CreateTrackingHandler(
    const std::vector<const Image*>& modelImages,
    StreamlineTractographyParameters::MODE mode,
    const Image* faImage = nullptr);

  /**
   * \brief The prior handler of a peak image, for
   *        itk::StreamlineTrackingFilter::SetTrackingPriorHandler.
   *
   * The handler works on a copy of the image, so it can be kept across runs while the image in
   * the data storage changes. Give it the parameters of each run with
   * SetTrackingPriorParameters().
   *
   * \throw mitk::Exception if the image is no peak image.
   */
  MITKFIBERTRACKING_EXPORT std::unique_ptr<TrackingDataHandler> CreateTrackingPriorHandler(const Image* priorImage);

  /**
   * \brief Gives a prior handler the parameters of a run: the tracker's, with the prior's own
   *        flips (m_PriorFlipX, ...) and without the data cutoff, which applies to the model
   *        and with "restrict to prior" would otherwise stop streamlines where the prior's
   *        peaks are short.
   */
  MITKFIBERTRACKING_EXPORT void SetTrackingPriorParameters(TrackingDataHandler& priorHandler,
                                                           const StreamlineTractographyParameters& parameters);

  /**
   * \brief The ROI images of a run: masks of any pixel type, non-zero inside, or for the target
   *        of the "label difference" endpoint constraint a label image. Null where the run has
   *        none.
   */
  struct StreamlineTractographyRois
  {
    const Image* seed = nullptr;
    const Image* mask = nullptr;
    const Image* stop = nullptr;
    const Image* target = nullptr;
    const Image* exclusion = nullptr;
  };

  /**
   * \brief Casts the ROI images to the float images the tracker reads and sets all five of
   *        them; a null image leaves the tracker without that region.
   */
  MITKFIBERTRACKING_EXPORT void SetTrackingRois(itk::StreamlineTrackingFilter& tracker, const StreamlineTractographyRois& rois);
}

#endif
