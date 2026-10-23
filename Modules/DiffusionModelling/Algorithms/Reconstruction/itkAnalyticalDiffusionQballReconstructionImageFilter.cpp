/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkAnalyticalDiffusionQballReconstructionImageFilter_cpp
#define itkAnalyticalDiffusionQballReconstructionImageFilter_cpp

#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkArray.h>
#include <vnl/vnl_vector.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionModellingHelperFunctions.h>

#include <locale>
#include "mitkPointShell.h"
#include "mitkShRegularization.h"
#include "mitkSingleShellGradients.h"

namespace itk {

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>::AnalyticalDiffusionQballReconstructionImageFilter() :
  m_GradientDirectionContainer(nullptr),
  m_NumberOfGradientDirections(0),
  m_NumberOfBaselineImages(1),
  m_Threshold(NumericTraits< ReferencePixelType >::NonpositiveMin()),
  m_BValue(-1),
  m_Lambda(0.0),
  m_NormalizationMethod(QBAR_STANDARD),
  m_NumberCoefficients(0),
  m_Delta1(0.001),
  m_Delta2(0.001)
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 );
  
}


template<
    class TReferenceImagePixelType,
    class TGradientImagePixelType,
    class TOdfPixelType,
    int ShOrder,
    int NrOdfDirections>
typename itk::AnalyticalDiffusionQballReconstructionImageFilter< TReferenceImagePixelType,TGradientImagePixelType,TOdfPixelType, ShOrder,NrOdfDirections>::OdfPixelType
itk::AnalyticalDiffusionQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType, ShOrder, NrOdfDirections>::Normalize( OdfPixelType odf, double b0 )
{
  switch( m_NormalizationMethod )
  {
  case QBAR_STANDARD:
  {
    TOdfPixelType sum = 0;

    for(int i=0; i<NrOdfDirections; i++)
    {
      sum += odf[i];
    }
    if(sum>0)
      odf /= sum;

    return odf;
    break;
  }
  case QBAR_ADC_ONLY:
  {
    for(int i=0; i<NrOdfDirections; i++)
    {
      odf[i] = ((TOdfPixelType)log((TOdfPixelType)b0)-odf[i])/m_BValue;
    }
    return odf;
    break;
  }
  case QBAR_RAW_SIGNAL:
  {
    return odf;
    break;
  }
  case QBAR_SOLID_ANGLE:
  {
    for(int i=0; i<NrOdfDirections; i++)
      odf[i] *= itk::Math::pi*4/NrOdfDirections;

    break;
  }
  }

  return odf;
}


template<
    class TReferenceImagePixelType,
    class TGradientImagePixelType,
    class TOdfPixelType,
    int ShOrder,
    int NrOdfDirections>
vnl_vector<TOdfPixelType> itk::AnalyticalDiffusionQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType, ShOrder, NrOdfDirections>::PreNormalize( vnl_vector<TOdfPixelType> vec, double b0 )
{
  switch( m_NormalizationMethod )
  {
  case QBAR_STANDARD:
  {
    int n = vec.size();
    double b0f = (double)b0;
    for(int i=0; i<n; i++)
    {
      vec[i] = vec[i]/b0f;
    }
    return vec;
    break;
  }
  case QBAR_ADC_ONLY:
  {
    int n = vec.size();
    for(int i=0; i<n; i++)
    {
      if (vec[i]<=0)
        vec[i] = 0.001;

      vec[i] = log(vec[i]);
    }
    return vec;
    break;
  }
  case QBAR_RAW_SIGNAL:
  {
    return vec;
    break;
  }
  case QBAR_SOLID_ANGLE:
  {
    int n = vec.size();
    double b0f = (double)b0;
    for(int i=0; i<n; i++)
    {
      vec[i] = vec[i]/b0f;

      if (vec[i]<0)
        vec[i] = m_Delta1;
      else if (vec[i]<m_Delta1)
        vec[i] = m_Delta1/2 + vec[i]*vec[i]/(2*m_Delta1);
      else if (vec[i]>=1)
        vec[i] = 1-m_Delta2/2;
      else if (vec[i]>=1-m_Delta2)
        vec[i] = 1-m_Delta2/2-(1-vec[i])*(1-vec[i])/(2*m_Delta2);

      vec[i] = log(-log(vec[i]));
    }
    return vec;
    break;
  }
  }

  return vec;
}

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
void AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>::BeforeThreadedGenerateData()
{
  // If we have more than 2 inputs, then each input, except the first is a
  // gradient image. The number of gradient images must match the number of
  // gradient directions.
  //const unsigned int numberOfInputs = this->GetNumberOfInputs();

  // Input must be an itk::VectorImage.
  std::string gradientImageClassName(
        this->ProcessObject::GetInput(0)->GetNameOfClass());
  if ( strcmp(gradientImageClassName.c_str(),"VectorImage") != 0 )
  {
    itkExceptionMacro( << "There is only one Gradient image. I expect that to be a VectorImage. " << "But its of type: " << gradientImageClassName );
  }

  this->ComputeReconstructionMatrix();

  typename GradientImagesType::Pointer img = static_cast< GradientImagesType * >( this->ProcessObject::GetInput(0) );

  m_BZeroImage = BZeroImageType::New();
  m_BZeroImage->SetSpacing( img->GetSpacing() );   // Set the image spacing
  m_BZeroImage->SetOrigin( img->GetOrigin() );     // Set the image origin
  m_BZeroImage->SetDirection( img->GetDirection() );  // Set the image direction
  m_BZeroImage->SetLargestPossibleRegion( img->GetLargestPossibleRegion());
  m_BZeroImage->SetBufferedRegion( img->GetLargestPossibleRegion() );
  m_BZeroImage->Allocate();

  m_CoefficientImage = CoefficientImageType::New();
  m_CoefficientImage->SetSpacing( img->GetSpacing() );   // Set the image spacing
  m_CoefficientImage->SetOrigin( img->GetOrigin() );     // Set the image origin
  m_CoefficientImage->SetDirection( img->GetDirection() );  // Set the image direction
  m_CoefficientImage->SetLargestPossibleRegion( img->GetLargestPossibleRegion());
  m_CoefficientImage->SetBufferedRegion( img->GetLargestPossibleRegion() );
  m_CoefficientImage->Allocate();
}

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
void AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread )
{
  typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetPrimaryOutput());

  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  ImageRegionIterator< BZeroImageType > oit2(m_BZeroImage, outputRegionForThread);
  oit2.GoToBegin();

  ImageRegionIterator< CoefficientImageType > oit4(m_CoefficientImage, outputRegionForThread);
  oit4.GoToBegin();

  typedef ImageRegionConstIterator< GradientImagesType > GradientIteratorType;
  typedef typename GradientImagesType::PixelType         GradientVectorType;
  typename GradientImagesType::Pointer gradientImagePointer = nullptr;

  // Would have liked a dynamic_cast here, but seems SGI doesn't like it
  // The enum will ensure that an inappropriate cast is not done
  gradientImagePointer = static_cast< GradientImagesType * >(
                           this->ProcessObject::GetInput(0) );

  GradientIteratorType git(gradientImagePointer, outputRegionForThread );
  git.GoToBegin();

  // Indices of the baseline and the gradient images
  const std::vector<unsigned int> &baselineind = m_BaselineIndices;
  const std::vector<unsigned int> &gradientind = m_GradientIndices;

  while( !git.IsAtEnd() )
  {
    GradientVectorType b = git.Get();

    // Average the baseline image pixels, in double so that the mean of several b0 volumes is not
    // truncated (it normalizes the signal, E = S/b0, and is the b0 output)
    double b0 = 0;
    for(unsigned int i = 0; i < baselineind.size(); ++i)
    {
      b0 += b[baselineind[i]];
    }
    b0 /= this->m_NumberOfBaselineImages;

    OdfPixelType odf(0.0);

    typename CoefficientImageType::PixelType coeffPixel(0.0);
    vnl_vector<TO> B(m_NumberOfGradientDirections);

    if( (b0 != 0) && (b0 >= m_Threshold) )
    {
      for( unsigned int i = 0; i< m_NumberOfGradientDirections; i++ )
      {
        B[i] = static_cast<TO>(b[gradientind[i]]);
      }

      B = PreNormalize(B, b0);
      if(m_NormalizationMethod == QBAR_SOLID_ANGLE)
      {
        vnl_vector<TO> coeffs(m_NumberCoefficients);
        coeffs = ( m_CoeffReconstructionMatrix * B );
        coeffs[0] += 1.0/(2.0*sqrt(itk::Math::pi));
        odf = ( m_SphericalHarmonicBasisMatrix * coeffs ).data_block();
        coeffPixel = coeffs.data_block();
      }
      else
      {
        // The constant 1/(4 pi) of the ODF (c0 = 1/(2 sqrt(pi))) belongs to the solid-angle ODF
        // only; here the SH image holds the same function as the ODF output.
        vnl_vector<TO> coeffs(m_NumberCoefficients);
        coeffs = ( m_CoeffReconstructionMatrix * B );
        coeffPixel = coeffs.data_block();
        odf = ( m_ReconstructionMatrix * B ).data_block();
      }
      odf = Normalize(odf, b0);
    }

    oit.Set( odf );
    oit2.Set( static_cast<BZeroPixelType>(b0) );
    oit4.Set(coeffPixel);
    ++oit;  // odf image iterator
    ++oit2; // b0 image iterator
    ++oit4; // coefficient image iterator
    ++git;  // Gradient  image iterator
  }
}

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
void AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>::ComputeReconstructionMatrix()
{
  // Derived per execution. b0 volumes are those whose b-value rounds to 0, as in the b-value map;
  // the model needs one shell.
  const mitk::SingleShellGradients shell = mitk::GetSingleShellGradients(m_GradientDirectionContainer, m_BValue);
  m_BaselineIndices = shell.baselineIndices;
  m_GradientIndices = shell.diffusionWeightedIndices;
  m_NumberOfBaselineImages = m_BaselineIndices.size();
  m_NumberOfGradientDirections = m_GradientIndices.size();

  m_NumberCoefficients = (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder;
  if( m_NumberOfGradientDirections < m_NumberCoefficients )
  {
    itkExceptionMacro( << "Not enough gradient directions supplied (" << m_NumberOfGradientDirections << "). At least " << (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder << " needed for SH-order " << ShOrder);
  }

  // Gradient preprocessing
  {
    // check for duplicate diffusion gradients
    bool warning = false;
    for(GradientDirectionContainerType::ConstIterator gdcit1 = this->m_GradientDirectionContainer->Begin();
        gdcit1 != this->m_GradientDirectionContainer->End(); ++gdcit1)
    {
      for(GradientDirectionContainerType::ConstIterator gdcit2 = this->m_GradientDirectionContainer->Begin();
          gdcit2 != this->m_GradientDirectionContainer->End(); ++gdcit2)
      {
        if(gdcit1.Value() == gdcit2.Value() && gdcit1.Index() != gdcit2.Index())
        {
          itkWarningMacro( << "Some of the Diffusion Gradients equal each other. Corresponding image data should be averaged before calling this filter." );
          warning = true;
          break;
        }
      }
      if (warning)
        break;
    }
  }

  // Create 3xM matrix Q that contains the gradient vectors in polar coordinates. Half-sphere
  // schemes need no mirrored copies: the even SH basis is antipodally symmetric, and a copy would
  // only double the data term and so halve the effective lambda.
  vnl_matrix<double> Q(3, m_NumberOfGradientDirections);
  for (unsigned int i = 0; i < m_NumberOfGradientDirections; ++i)
  {
    const auto &direction = shell.directions[i];
    double cart[3];
    mitk::gradients::Cart2Sph(direction[0], direction[1], direction[2], cart);
    Q(0,i) = cart[0];
    Q(1,i) = cart[1];
    Q(2,i) = cart[2];
  }

  // SH basis B (MITK's basis, in double; the per-voxel products below run in float), the
  // Laplace-Beltrami eigenvalues L, the regularization and the Funk-Radon factors P per coefficient
  const vnl_matrix<double> B = mitk::sh::ComputeSphericalHarmonicsBasis(Q, ShOrder);
  vnl_matrix<double> L(m_NumberCoefficients,m_NumberCoefficients, 0);
  vnl_matrix<double> penalty(m_NumberCoefficients,m_NumberCoefficients, 0);
  vnl_matrix<double> P(m_NumberCoefficients,m_NumberCoefficients, 0);
  for(int k=0; k<=ShOrder; k+=2)
  {
    for(int m=-k; m<=k; m++)
    {
      int j = (k*k + k + 2)/2 + m - 1;
      L(j,j) = -k*(k+1);
      penalty(j,j) = mitk::LaplaceBeltramiPenalty(k, m);
      P(j,j) = mitk::sh::legendre0(k);
    }
  }

  const vnl_matrix<double> B_transpose(B.transpose());
  vnl_matrix_inverse<double> pseudoInverse( B_transpose * B + penalty*m_Lambda );
  vnl_matrix<double> coeffReconstruction = pseudoInverse.pinverse() * B_transpose;

  switch(m_NormalizationMethod)
  {
  case QBAR_ADC_ONLY:
  case QBAR_RAW_SIGNAL:
    break;
  case QBAR_STANDARD:
  {
    coeffReconstruction = P * coeffReconstruction;
    break;
  }
  case QBAR_SOLID_ANGLE:
  {
    coeffReconstruction = (1.0/(8.0*itk::Math::pi)) * P * L * coeffReconstruction;
    break;
  }
  }

  // needed to calculate the ODF values from the SH coefficients
  const vnl_matrix<double> &odfDirections = mitk::pointshell::Directions(NrOdfDirections);
  vnl_matrix<double> odfSpherical(3, NrOdfDirections);
  for (int i = 0; i < NrOdfDirections; ++i)
  {
    double spherical[3];
    mitk::gradients::Cart2Sph(odfDirections(0,i), odfDirections(1,i), odfDirections(2,i), spherical);
    for (unsigned int r = 0; r < 3; ++r)
      odfSpherical(r,i) = spherical[r];
  }
  const vnl_matrix<double> odfBasis = mitk::sh::ComputeSphericalHarmonicsBasis(odfSpherical, ShOrder);
  const auto toFloat = [](const vnl_matrix<double> &matrix) {
    vnl_matrix<float> result(matrix.rows(), matrix.cols());
    for (unsigned int r = 0; r < matrix.rows(); ++r)
      for (unsigned int c = 0; c < matrix.cols(); ++c)
        result(r,c) = static_cast<float>(matrix(r,c));
    return result;
  };
  m_CoeffReconstructionMatrix = toFloat(coeffReconstruction);
  m_SphericalHarmonicBasisMatrix = toFloat(odfBasis);
  m_ReconstructionMatrix = toFloat(odfBasis * coeffReconstruction);
}

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
void AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>
::SetGradientImage(const GradientDirectionContainerType *gradientDirection,
                   const GradientImagesType *gradientImage )
{
  // Copy Gradient Direction Container
  this->m_GradientDirectionContainer = GradientDirectionContainerType::New();
  for(GradientDirectionContainerType::ConstIterator it = gradientDirection->Begin();
      it != gradientDirection->End(); it++)
  {
    this->m_GradientDirectionContainer->push_back(it.Value());
  }

  if (m_BValue<0)
    mitkThrow() << "B-value needs to be set before gradient image! " << m_BValue;

  // The gradients are classified and normalized per execution (ComputeReconstructionMatrix).
  // ensure that the gradient image we received has as many components as
  // the number of gradient directions
  if( gradientImage->GetVectorLength() != gradientDirection->Size() )
  {
    itkExceptionMacro( << gradientDirection->Size() << " directions specified but image has "
                       << gradientImage->GetVectorLength() << " components.");
  }

  this->ProcessObject::SetNthInput( 0, const_cast< GradientImagesType* >(gradientImage) );
}

template< class T, class TG, class TO, int ShOrder, int NrOdfDirections>
void AnalyticalDiffusionQballReconstructionImageFilter<T,TG,TO,ShOrder,NrOdfDirections>
::PrintSelf(std::ostream& os, Indent indent) const
{
  std::locale C("C");
  std::locale originalLocale = os.getloc();
  os.imbue(C);

  Superclass::PrintSelf(os,indent);

  os << indent << "OdfReconstructionMatrix: " << m_ReconstructionMatrix << std::endl;
  if ( m_GradientDirectionContainer )
    os << indent << "GradientDirectionContainer: " << m_GradientDirectionContainer << std::endl;
  else
    os << indent << "GradientDirectionContainer: (Gradient directions not set)" << std::endl;

  os << indent << "NumberOfGradientDirections: " <<  m_NumberOfGradientDirections << std::endl;
  os << indent << "NumberOfBaselineImages: " << m_NumberOfBaselineImages << std::endl;
  os << indent << "Threshold for reference B0 image: " << m_Threshold << std::endl;
  os << indent << "BValue: " << m_BValue << std::endl;
  os.imbue( originalLocale );
}

}

#endif // itkAnalyticalDiffusionQballReconstructionImageFilter_cpp
