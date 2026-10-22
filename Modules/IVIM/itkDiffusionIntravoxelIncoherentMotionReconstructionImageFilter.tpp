/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace itk
{
  template <class TIn, class TOut>
  DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::DiffusionIntravoxelIncoherentMotionReconstructionImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(3);
    for (unsigned int i = 0; i < 3; ++i)
      this->SetNthOutput(i, OutputImageType::New().GetPointer());
  }

  template <class TIn, class TOut>
  void DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::VerifyPreconditions() const
  {
    Superclass::VerifyPreconditions();

    if (m_GradientDirections.IsNull())
      itkExceptionMacro("Gradient directions not set.");
    if (m_BValue <= 0)
      itkExceptionMacro("Reference b-value not set.");
  }

  template <class TIn, class TOut>
  void DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::VerifyInputInformation() const
  {
    Superclass::VerifyInputInformation();

    if (m_GradientDirections->Size() != this->GetInput()->GetNumberOfComponentsPerPixel())
      itkExceptionMacro("The input has " << this->GetInput()->GetNumberOfComponentsPerPixel() << " volumes but "
                                         << m_GradientDirections->Size() << " gradient directions.");
  }

  template <class TIn, class TOut>
  void DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::BeforeThreadedGenerateData()
  {
    m_Protocol = mitk::ivim::MakeProtocol(m_GradientDirections, m_BValue);
  }

  template <class TIn, class TOut>
  void DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::DynamicThreadedGenerateData(
    const OutputImageRegionType& outputRegionForThread)
  {
    ImageRegionConstIterator<InputImageType> inputIt(this->GetInput(), outputRegionForThread);
    ImageRegionIterator<OutputImageType> fIt(this->GetOutput(0), outputRegionForThread);
    ImageRegionIterator<OutputImageType> dIt(this->GetOutput(1), outputRegionForThread);
    ImageRegionIterator<OutputImageType> dStarIt(this->GetOutput(2), outputRegionForThread);

    vnl_vector<double> signal(this->GetInput()->GetNumberOfComponentsPerPixel());
    for (; !inputIt.IsAtEnd(); ++inputIt, ++fIt, ++dIt, ++dStarIt)
    {
      const auto pixel = inputIt.Get();
      for (unsigned int i = 0; i < signal.size(); ++i)
        signal[i] = pixel[i];

      const auto fit = mitk::ivim::FitVoxel(signal, m_Protocol, m_FitParameters);
      fIt.Set(fit.f);
      dIt.Set(fit.d);
      dStarIt.Set(fit.dStar);
    }
  }

  template <class TIn, class TOut>
  void DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<TIn, TOut>::PrintSelf(std::ostream& os,
                                                                                          Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "BValue: " << m_BValue << std::endl;
    os << indent << "Method: " << static_cast<int>(m_FitParameters.method) << std::endl;
    if (m_GradientDirections)
      os << indent << "GradientDirections: " << m_GradientDirections << std::endl;
    else
      os << indent << "GradientDirections: (not set)" << std::endl;
  }
}
