/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef DIFFUSIONKURTOSISRECONSTRUCTIONIMAGEFILTER_CXX
#define DIFFUSIONKURTOSISRECONSTRUCTIONIMAGEFILTER_CXX

#include "itkDiffusionKurtosisReconstructionImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkComposeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

/** K bounds with lower >= upper admit no K, and equal bounds give the penalty band zero width. */
static void VerifyKurtosisBounds(const vnl_vector_fixed<double, 2>& bounds)
{
  if (!(bounds[0] < bounds[1]))
    itkGenericExceptionMacro("Lower kurtosis bound " << bounds[0] << " is not below the upper bound " << bounds[1] << ".");
}

namespace itk
{
  /** The measurements a kurtosis fit uses. S0 is the mean of the b=0 measurements, not the first
      volume, so the volume order does not matter; it is fitted instead when asked to or when there
      is no b=0 measurement. The log-scale fit leaves out non-positive measurements, which have no
      logarithm. */
  struct KurtosisFitData
  {
    vnl_vector<double> measurements;
    vnl_vector<double> bvalues;
    double bzero = 0;
    unsigned int bzero_count = 0;
    bool fit_bzero = false;
  };
}

template< class TInputPixelType>
static itk::KurtosisFitData AssembleKurtosisFitData( const itk::VariableLengthVector< TInputPixelType > &input,
                                                     const vnl_vector<double>& bvalues,
                                                     const itk::KurtosisFitConfiguration& kf_config)
{
  std::vector<double> measurements;
  std::vector<double> fit_bvalues;
  itk::KurtosisFitData data;
  for( unsigned int i=0; i<bvalues.size(); ++i )
  {
    if( bvalues[i] < vnl_math::eps )
    {
      data.bzero += input.GetElement(i);
      ++data.bzero_count;
    }
    // skip bvalues higher than the limit provided, if the flag is activated
    else if( ( !kf_config.exclude_high_b || bvalues[i] <= kf_config.b_upper_threshold )
             && ( kf_config.fit_scale != itk::LOGARITHMIC || input.GetElement(i) >= vnl_math::eps ) )
    {
      measurements.push_back( input.GetElement(i) );
      fit_bvalues.push_back( bvalues[i] );
    }
  }

  if( data.bzero_count > 0 )
    data.bzero /= data.bzero_count;
  data.fit_bzero = kf_config.omit_bzero || data.bzero_count == 0;
  data.measurements = vnl_vector<double>( measurements.data(), measurements.size() );
  data.bvalues = vnl_vector<double>( fit_bvalues.data(), fit_bvalues.size() );
  return data;
}

/** True when Levenberg-Marquardt gave up or never started; x then still holds the start values. */
static bool KurtosisFitFailed( const vnl_levenberg_marquardt& fit, const vnl_vector<double>& x )
{
  const auto code = fit.get_failure_code();
  return code == vnl_nonlinear_minimizer::ERROR_FAILURE || code == vnl_nonlinear_minimizer::ERROR_DODGY_INPUT
         || !x.is_finite();
}

/** Fits one voxel; result holds the start values (D, K, S0) on input and the fitted (D, K, S0) on
    output, S0 being the b=0 mean when it is not fitted. A voxel without signal, with fewer usable
    measurements than parameters or whose fit fails gets 0, as voxels outside the mask do. */
template< class TInputPixelType>
static void FitSingleVoxel( const itk::VariableLengthVector< TInputPixelType > &input,
                            const vnl_vector<double>& bvalues,
                            vnl_vector<double>& result,
                            const itk::KurtosisFitConfiguration& kf_config)
{
  // check for length
  assert( input.Size() == bvalues.size() );
  assert( result.size() == 3 );

  const auto data = AssembleKurtosisFitData( input, bvalues, kf_config );

  MITK_DEBUG("KurtosisFilter.FitSingleVoxel.Meas") << data.measurements;
  MITK_DEBUG("KurtosisFilter.FitSingleVoxel.Bval") << data.bvalues;

  const unsigned int num_params = data.fit_bzero ? 3 : 2;
  if( data.measurements.size() < num_params || data.measurements.max_value() < vnl_math::eps
      || ( !data.fit_bzero && data.bzero < vnl_math::eps ) )
  {
    result.fill(0);
    return;
  }

  // perform fit on data vectors
  if( data.fit_bzero )
  {
   itk::kurtosis_fit_omit_unweighted kurtosis_cost_fn( data.measurements.size() );
   // configuration
   kurtosis_cost_fn.set_fit_logscale( static_cast<bool>(kf_config.fit_scale) );
   kurtosis_cost_fn.initialize( data.measurements, data.bvalues );

   if( kf_config.use_K_limits)
   {
    kurtosis_cost_fn.set_K_bounds( kf_config.K_limits );
   }

   const bool logscale = kf_config.fit_scale == itk::LOGARITHMIC;
   if( logscale )
     result[2] = log( std::max( result[2], 1.0 ) );

   vnl_levenberg_marquardt nonlinear_fit( kurtosis_cost_fn );
   nonlinear_fit.minimize( result );
   if( KurtosisFitFailed( nonlinear_fit, result ) )
   {
     result.fill(0);
     return;
   }
   kurtosis_cost_fn.clamp_to_bounds( result );

   if( logscale )
     result[2] = exp( result[2] );
  }
  else
  {
   itk::kurtosis_fit_lsq_function kurtosis_cost_fn( data.measurements.size() );
   // configuration
   kurtosis_cost_fn.set_fit_logscale( static_cast<bool>(kf_config.fit_scale) );
   kurtosis_cost_fn.initialize( data.measurements, data.bvalues );
   kurtosis_cost_fn.set_bzero( data.bzero );

   if( kf_config.use_K_limits)
   {
    kurtosis_cost_fn.set_K_bounds( kf_config.K_limits );
   }

   vnl_vector<double> d_and_k = result.extract( 2 );
   vnl_levenberg_marquardt nonlinear_fit( kurtosis_cost_fn );
   nonlinear_fit.minimize( d_and_k );
   if( KurtosisFitFailed( nonlinear_fit, d_and_k ) )
   {
     result.fill(0);
     return;
   }
   kurtosis_cost_fn.clamp_to_bounds( d_and_k );
   result[0] = d_and_k[0];
   result[1] = d_and_k[1];
   result[2] = data.bzero;
  }

  MITK_DEBUG("KurtosisFilter.FitSingleVoxel.Rslt") << result;
}



template< class TInputPixelType, class TOutputPixelType>
itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::DiffusionKurtosisReconstructionImageFilter()
  : m_ReferenceBValue(-1),
    m_OmitBZero(false),
    m_ApplyPriorSmoothing(false),
    m_SmoothingSigma(1.5),
    m_UseKBounds( false ),
    m_MaxFitBValue( 3000 ),
    m_ScaleForFitting( STRAIGHT )
{
  this->m_InitialPosition = vnl_vector<double>(3, 0);
  this->m_InitialPosition[2] = 1000.0; // S_0
  this->m_InitialPosition[0] = 0.001; // D
  this->m_InitialPosition[1] = 1; // K

  this->m_KurtosisBounds.fill(0);

  // single input image
  this->SetNumberOfRequiredInputs(1);

  // two output images (D, K)
  this->SetNumberOfRequiredOutputs(2);
  typename OutputImageType::Pointer outputPtr1 = OutputImageType::New();
  typename OutputImageType::Pointer outputPtr2 = OutputImageType::New();

  this->SetNthOutput(0, outputPtr1.GetPointer() );
  this->SetNthOutput(1, outputPtr2.GetPointer() );

  
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::VerifyPreconditions() const
{
  Superclass::VerifyPreconditions();

  if (this->m_BValues.empty())
    itkExceptionMacro("Gradient directions not set.");
  if (this->m_UseKBounds)
    VerifyKurtosisBounds(this->m_KurtosisBounds);
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::GenerateOutputInformation()
{
  // first call superclass
  Superclass::GenerateOutputInformation();
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::SetImageMask(MaskImageType::Pointer mask)
{
  this->m_MaskImage = mask;
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::BeforeThreadedGenerateData()
{
  // if we have a region set, convert it to a mask image, which is to be used as default for telling
  // we need to set the image anyway, so by default the mask is overall 1
  if( m_MaskImage.IsNull() )
  {
    m_MaskImage = MaskImageType::New();
    m_MaskImage->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
    m_MaskImage->CopyInformation( this->GetInput() );
    m_MaskImage->Allocate();

    if( this->m_MapOutputRegion.GetNumberOfPixels() > 0 )
    {

      m_MaskImage->FillBuffer(0);

      typedef itk::ImageRegionIteratorWithIndex< MaskImageType > MaskIteratorType;
      MaskIteratorType maskIter( this->m_MaskImage, this->m_MapOutputRegion );
      maskIter.GoToBegin();

      while( !maskIter.IsAtEnd() )
      {
        maskIter.Set( 1 );
        ++maskIter;
      }

    }
    else
    {
      m_MaskImage->FillBuffer(1);
    }


  }


  // apply smoothing to the input image
  if( this->m_ApplyPriorSmoothing )
  {
    // filter typedefs
    typedef itk::DiscreteGaussianImageFilter< itk::Image<TInputPixelType, 3 >, itk::Image<TInputPixelType, 3 > > GaussianFilterType;
    typedef itk::VectorIndexSelectionCastImageFilter< InputImageType, itk::Image<TInputPixelType, 3 > > IndexSelectionType;
    typedef itk::ComposeImageFilter< itk::Image<TInputPixelType, 3>, InputImageType > ComposeFilterType;


    auto vectorImage = this->GetInput();
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetInput( vectorImage );

    typename ComposeFilterType::Pointer vec_composer = ComposeFilterType::New();

    for( unsigned int i=0; i<vectorImage->GetVectorLength(); ++i)
    {
      typename GaussianFilterType::Pointer gaussian_filter = GaussianFilterType::New();

      indexSelectionFilter->SetIndex( i );

      gaussian_filter->SetInput( indexSelectionFilter->GetOutput() );
      gaussian_filter->SetVariance( m_SmoothingSigma );

      vec_composer->SetInput(i, gaussian_filter->GetOutput() );

      gaussian_filter->Update();
    }

    try
    {
      vec_composer->Update();
    }
    catch(const itk::ExceptionObject &e)
    {
      mitkThrow() << "[VectorImage.GaussianSmoothing] !! Failed with ITK Exception: " << e.what();
    }

    this->m_ProcessedInputImage = vec_composer->GetOutput();
  }
  else
  {
    this->m_ProcessedInputImage = const_cast<InputImageType*>( this->GetInput() );
  }
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::AfterThreadedGenerateData()
{
 /* // initialize buffer to zero overall, but don't forget the requested region pointer
  for( unsigned int i=0; i<this->GetNumberOfOutputs(); ++i)
  {
    typename OutputImageType::Pointer output = this->GetOutput(i);

    // after generating, set the buffered region to max, we have taken care about filling it with valid data in
    // UpdateOutputInformation, if not set, a writer (or any filter using the LargestPossibleRegion() during update may
    // complain about not getting the requested region
    output->SetBufferedRegion( output->GetLargestPossibleRegion() );

    std::cout << "[DiffusionKurtosisReconstructionImageFilter.After]" << output->GetLargestPossibleRegion() << std::endl;
  }*/

}

template< class TInputPixelType, class TOutputPixelType>
typename itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::KurtosisSnapshot
itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::GetSnapshot(const itk::VariableLengthVector<TInputPixelType> &input, GradientDirectionContainerType::ConstPointer gradients, float bvalue, KurtosisFitConfiguration kf_conf)
{
  // initialize bvalues from reference value and the gradients provided on input
  this->SetReferenceBValue(bvalue);
  this->SetGradientDirections( gradients );

  // call the other method
  return this->GetSnapshot( input, this->m_BValues, kf_conf );
}


template< class TInputPixelType, class TOutputPixelType>
typename itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::KurtosisSnapshot
itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::GetSnapshot(const itk::VariableLengthVector<TInputPixelType> &input, vnl_vector<double> bvalues, KurtosisFitConfiguration kf_conf)
{
  if (kf_conf.use_K_limits)
    VerifyKurtosisBounds(kf_conf.K_limits);

  vnl_vector<double> fit = this->m_InitialPosition;
  FitSingleVoxel( input, bvalues, fit, kf_conf );

  const auto data = AssembleKurtosisFitData( input, bvalues, kf_conf );

  // assembly result snapshot
  KurtosisSnapshot result;
  result.m_D = fit[0];
  result.m_K = fit[1];
  result.m_Bzero = data.bzero;
  result.m_fittedBZero = data.fit_bzero;

  if( data.fit_bzero )
  {
    result.m_BzeroFit = fit[2];
    if( data.bzero_count > 0 )
      result.m_f = 1 - fit[2] / std::fmax(0.01, data.bzero);
  }
  else
    result.m_f = 1;

  result.fit_bvalues = data.bvalues;
  result.fit_measurements = data.measurements;

  result.bvalues = bvalues;
  result.measurements = vnl_vector<double>( input.Size() );
  for( unsigned int i=0; i<input.Size(); ++i )
    result.measurements[i] = input.GetElement(i);

  return result;
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::SetGradientDirections(GradientDirectionContainerType::ConstPointer gradients)
{
  if( this->m_ReferenceBValue < 0)
  {
    itkExceptionMacro( << "Specify reference b-value prior to setting the gradient directions." );
  }

  if( gradients->Size() == 0 )
  {
    itkExceptionMacro( << "Empty gradient directions container retrieved" );
  }

  vnl_vector<double> vnl_bvalues( gradients->Size(), 0 );

  // initialize the b-values
  auto grIter = gradients->Begin();
  unsigned int index = 0;
  while( grIter != gradients->End() )
  {
    const double twonorm = grIter.Value().two_norm();
    vnl_bvalues( index++ ) = this->m_ReferenceBValue * twonorm * twonorm;

    ++grIter;
  }

  this->m_BValues = vnl_bvalues;


}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::SetInitialSolution(const vnl_vector<double>& x0 )
{
  assert( x0.size() == 2 || x0.size() == 3 );

  this->m_InitialPosition[0] = x0[0];
  this->m_InitialPosition[1] = x0[1];
  if( x0.size() == 3 )
    this->m_InitialPosition[2] = x0[2];
}

template< class TInputPixelType, class TOutputPixelType>
void itk::DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread /*threadId*/)
{
  typename OutputImageType::Pointer dImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  itk::ImageRegionIteratorWithIndex< OutputImageType > dImageIt(dImage, outputRegionForThread);
  dImageIt.GoToBegin();

  typename OutputImageType::Pointer kImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(1));
  itk::ImageRegionIteratorWithIndex< OutputImageType > kImageIt(kImage, outputRegionForThread);
  kImageIt.GoToBegin();

  typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InputIteratorType;
  InputIteratorType inputIter( m_ProcessedInputImage, outputRegionForThread );
  inputIter.GoToBegin();

  typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType > MaskIteratorType;
  MaskIteratorType maskIter( this->m_MaskImage, outputRegionForThread );
  maskIter.GoToBegin();

  KurtosisFitConfiguration fit_config;
  fit_config.omit_bzero = this->m_OmitBZero;
  fit_config.fit_scale = this->m_ScaleForFitting;

  fit_config.use_K_limits = this->m_UseKBounds;
  fit_config.K_limits = this->m_KurtosisBounds;

  while( !inputIter.IsAtEnd() )
  {
    // set (reset) each iteration
    vnl_vector<double> result = this->m_InitialPosition;

    // fit single voxel (if inside mask )
    if( maskIter.Get() > 0 )
    {
      FitSingleVoxel( inputIter.Get(), this->m_BValues, result, fit_config );
    }
    else
    {
      result.fill(0);
    }

    // regardless the fit type, the parameters are always in the first two position
    // of the results vector
    dImageIt.Set( result[0] );
    kImageIt.Set( result[1] );

    //std::cout << "[Kurtosis.Fit]" << inputIter.GetIndex() << " --> " << dImageIt.GetIndex() << " result: " << result << std::endl;

    ++maskIter;
    ++inputIter;
    ++dImageIt;
    ++kImageIt;
  }

}

#endif // guards
