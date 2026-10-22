/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef DIFFUSIONKURTOSISRECONSTRUCTIONIMAGEFILTER_H
#define DIFFUSIONKURTOSISRECONSTRUCTIONIMAGEFILTER_H

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <mitkKurtosisFit.h>

namespace itk
{
  /**
    @class DiffusionKurtosisReconstructionImageFilter
    @brief Maps of the kurtosis model parameters D (output 0) and K (output 1), fitted voxel by voxel
    with mitk::kurtosis::FitVoxel.
    */
  template <class TInputPixelType, class TOutputPixelType>
  class DiffusionKurtosisReconstructionImageFilter
    : public ImageToImageFilter<VectorImage<TInputPixelType, 3>, Image<TOutputPixelType, 3>>
  {
  public:
    ITK_DISALLOW_COPY_AND_MOVE(DiffusionKurtosisReconstructionImageFilter);

    using Self = DiffusionKurtosisReconstructionImageFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<VectorImage<TInputPixelType, 3>, Image<TOutputPixelType, 3>>;

    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);
    itkOverrideGetNameOfClassMacro(DiffusionKurtosisReconstructionImageFilter);

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;
    using MaskImageType = Image<short, 3>;
    using GradientDirectionContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

    /** The reference b-value of the gradients; no default, the filter throws when it is not set. */
    itkSetMacro(ReferenceBValue, double);

    /** The gradient of each volume of the input, scaled so that b = reference b-value * |g|^2. */
    itkSetConstObjectMacro(GradientDirections, GradientDirectionContainerType);

    /** Restricts the fit to voxels where the mask is not 0, others get 0; the mask must have the
        geometry of the input. */
    itkSetConstObjectMacro(ImageMask, MaskImageType);

    /** The variance in mm^2, not the sigma, of a Gaussian smoothing of each volume before the fit;
        0 (the default) means no smoothing. */
    itkSetMacro(SmoothingVariance, double);

    void SetFitConfiguration(const mitk::kurtosis::FitConfiguration& configuration)
    {
      m_FitConfiguration = configuration;
      this->Modified();
    }

    const mitk::kurtosis::FitConfiguration& GetFitConfiguration() const { return m_FitConfiguration; }

  protected:
    DiffusionKurtosisReconstructionImageFilter();
    ~DiffusionKurtosisReconstructionImageFilter() override = default;

    void VerifyPreconditions() const override;
    void VerifyInputInformation() const override;
    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) override;

  private:
    double m_ReferenceBValue{0};
    GradientDirectionContainerType::ConstPointer m_GradientDirections;
    MaskImageType::ConstPointer m_ImageMask;
    double m_SmoothingVariance{0};
    mitk::kurtosis::FitConfiguration m_FitConfiguration;

    /** Made for each run and only read by the threads. */
    vnl_vector<double> m_BValues;
    typename InputImageType::ConstPointer m_ProcessedInputImage;
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionKurtosisReconstructionImageFilter.cpp"
#endif

#endif // DIFFUSIONKURTOSISRECONSTRUCTIONIMAGEFILTER_H
