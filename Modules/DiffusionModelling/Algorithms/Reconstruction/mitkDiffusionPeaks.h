/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionPeaks_h
#define mitkDiffusionPeaks_h

#include <MitkDiffusionModellingExports.h>
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <atomic>

/** Peak extraction from SH and tensor images, as the peak extraction view runs it. */
namespace mitk
{
  /** \brief Length of the extracted ODF peaks. */
  enum class PeakNormalization
  {
    None,        ///< the ODF value in the peak direction
    MaxVector,   ///< relative to the largest peak of the voxel
    SingleVector ///< unit vectors
  };

  struct PeakParameters
  {
    unsigned int maxPeaks = 3;
    /** Peaks below this fraction of the largest peak of their voxel are dropped. */
    double relativeThreshold = 0.4;
    /** Peaks whose ODF value is below are dropped. */
    double absoluteThreshold = 0.03;
    /** In degrees: of two peaks closer than this, the smaller one is dropped. */
    double angularThreshold = 15;
    bool scaleByGfa = false;
    PeakNormalization normalization = PeakNormalization::MaxVector;
  };

  struct PeakResult
  {
    /** Peaks in world coordinates, three components per peak. */
    PeakImage::Pointer peaks;
    /** unsigned char image with the number of peaks per voxel. */
    Image::Pointer numberOfPeaks;
  };

  /**
   * \brief Maxima of the ODFs of an SH image (itk::OdfMaximaExtractionFilter on the ODF directions),
   * evaluated in the SH convention stored in the image.
   * \param mask optional; voxels where it is 0 get no peaks. Any pixel type, cast to unsigned char.
   * \param stop optional; once set, the extraction starts no further block of voxels and ends with
   * itk::ProcessAborted after the blocks already running.
   * \throw mitk::Exception without SH image or for an SH order other than 2 to 12;
   * itk::ExceptionObject for a mask of another geometry.
   */
  MITKDIFFUSIONMODELLING_EXPORT PeakResult ExtractPeaks(const ShImage *shImage,
                                                        const PeakParameters &parameters,
                                                        const Image *mask = nullptr,
                                                        const std::atomic<bool> *stop = nullptr);

  /**
   * \brief Principal eigenvectors of a tensor image (itk::DiffusionTensorPrincipalDirectionImageFilter).
   * \param faThreshold voxels with a lower FA get no peak.
   * \param normalize unit vectors; otherwise the vectors are scaled by the FA.
   * \param mask optional; voxels where it is 0 get no peak.
   * \throw mitk::Exception without tensor image; itk::ExceptionObject for a mask of another geometry.
   */
  MITKDIFFUSIONMODELLING_EXPORT PeakResult ExtractTensorPeaks(const TensorImage *tensors,
                                                              double faThreshold,
                                                              bool normalize,
                                                              const Image *mask = nullptr);
}

#endif
