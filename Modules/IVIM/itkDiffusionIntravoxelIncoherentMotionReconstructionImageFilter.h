/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter_h
#define itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <mitkIvimFit.h>
#include <mitkMaskGeometry.h>

namespace itk
{
  /** \class DiffusionIntravoxelIncoherentMotionReconstructionImageFilter
      \brief Maps of the IVIM parameters f (output 0), D (output 1) and D* (output 2), fitted voxel by
      voxel with a mitk::ivim::VoxelFitter per thread; the voxels are independent, so the filter runs
      multithreaded. */
  template <class TInputPixelType, class TOutputPixelType>
  class DiffusionIntravoxelIncoherentMotionReconstructionImageFilter
    : public ImageToImageFilter<VectorImage<TInputPixelType, 3>, Image<TOutputPixelType, 3>>
  {
  public:
    ITK_DISALLOW_COPY_AND_MOVE(DiffusionIntravoxelIncoherentMotionReconstructionImageFilter);

    using Self = DiffusionIntravoxelIncoherentMotionReconstructionImageFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<VectorImage<TInputPixelType, 3>, Image<TOutputPixelType, 3>>;

    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);
    itkOverrideGetNameOfClassMacro(DiffusionIntravoxelIncoherentMotionReconstructionImageFilter);

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;
    using MaskImageType = Image<short, 3>;
    using GradientDirectionContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

    /** The gradient of each volume of the input, scaled so that b = reference b-value * |g|^2. */
    itkSetConstObjectMacro(GradientDirections, GradientDirectionContainerType);

    /** The reference b-value of the gradients; no default, the filter throws when it is not set. */
    itkSetMacro(BValue, double);

    /** Restricts the fit to voxels where the mask is not 0, others get 0; the mask must have the
        voxel grid of the input (mitk::HasVoxelGridOf). */
    itkSetConstObjectMacro(ImageMask, MaskImageType);

    void SetFitParameters(const mitk::ivim::FitParameters& parameters)
    {
      m_FitParameters = parameters;
      this->Modified();
    }

    const mitk::ivim::FitParameters& GetFitParameters() const { return m_FitParameters; }

  protected:
    DiffusionIntravoxelIncoherentMotionReconstructionImageFilter();
    ~DiffusionIntravoxelIncoherentMotionReconstructionImageFilter() override = default;

    void PrintSelf(std::ostream& os, Indent indent) const override;
    void VerifyPreconditions() const override;
    void VerifyInputInformation() const override;
    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) override;

  private:
    GradientDirectionContainerType::ConstPointer m_GradientDirections;
    double m_BValue{0};
    MaskImageType::ConstPointer m_ImageMask;
    mitk::ivim::FitParameters m_FitParameters;

    /** Made from the gradients for each run and only read by the threads. */
    mitk::ivim::Protocol m_Protocol;
  };
}

#include "itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.tpp"

#endif
