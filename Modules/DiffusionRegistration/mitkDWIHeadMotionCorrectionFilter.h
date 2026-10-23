/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDWIHEADMOTIONCORRECTIONFILTER_H
#define MITKDWIHEADMOTIONCORRECTIONFILTER_H

#include <mitkImageToImageFilter.h>

#include <itkVectorImage.h>
#include <MitkDiffusionRegistrationExports.h>

namespace mitk
{

/**
 * @class DWIHeadMotionCorrectionFilter
 *
 * @brief Corrects head motion in a DWI by an affine registration of every volume to a reference.
 *
 * The reference is the first volume of the lowest b-value shell (usually a b0 volume). Each other
 * volume is registered to it (MatchPoint's multi-modal affine default algorithm) and resampled with
 * B-spline interpolation; the output keeps the volume order of the input, and each gradient is
 * turned by the rotation part of its volume's transform. Throws for an input that is not a DWI.
 */
class MITKDIFFUSIONREGISTRATION_EXPORT DWIHeadMotionCorrectionFilter
    : public ImageToImageFilter
{
public:

  mitkClassMacro( DWIHeadMotionCorrectionFilter,
                  ImageToImageFilter )

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  using DiffusionPixelType = short;
  using InputImageType = Superclass::InputImageType;
  using OutputImageType = Superclass::OutputImageType;
  using ITKDiffusionImageType = itk::VectorImage<DiffusionPixelType, 3>;
  using ITKDiffusionVolumeType = itk::Image<DiffusionPixelType, 3>;

  /** The corrected DWI, the same image as GetOutput(). */
  mitk::Image::Pointer GetCorrectedImage();

protected:
  DWIHeadMotionCorrectionFilter() = default;
  ~DWIHeadMotionCorrectionFilter() override = default;

  void GenerateData() override;

};

} //end namespace mitk


#endif // MITKDWIHEADMOTIONCORRECTIONFILTER_H
