/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkDiffusionMultiShellQballReconstructionImageFilter_cpp
#define itkDiffusionMultiShellQballReconstructionImageFilter_cpp

#include <itkDiffusionMultiShellQballReconstructionImageFilter.h>

#include <limits>
#include <optional>
#include <vnl/vnl_sse.h>
#include <mitkPointShell.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include "mitkShRegularization.h"

namespace itk {

template< class T, class TG, class TO, int L, int NODF>
DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::DiffusionMultiShellQballReconstructionImageFilter() :
  m_Interpolation_Flag(false),
  m_MaxDirections(0),
  m_GradientDirectionContainer(nullptr),
  m_NumberOfGradientDirections(0),
  m_NumberOfBaselineImages(0),
  m_Threshold(0),
  m_BZeroImage(nullptr),
  m_CoefficientImage(nullptr),
  m_BValue(1.0),
  m_Lambda(0.0),
  m_ComputeOdf(true)
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 );

  
}


template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::SetGradientImage(const GradientDirectionContainerType *gradientDirection
                    , const GradientImagesType *gradientImage
                    , float bvalue)
{
  m_BValue = bvalue;
  m_NumberOfBaselineImages = 0;

  this->m_GradientDirectionContainer = GradientDirectionContainerType::New();
  for(GradientDirectionContainerType::ConstIterator it = gradientDirection->Begin();
    it != gradientDirection->End(); it++)
  {
    this->m_GradientDirectionContainer->push_back(it.Value());
  }


  if(m_BValueMap.size() == 0){
    itkWarningMacro(<< "DiffusionMultiShellQballReconstructionImageFilter.cpp : no GradientIndexMapAvalible");

    GradientDirectionContainerType::ConstIterator gdcit;
    for( gdcit = m_GradientDirectionContainer->Begin(); gdcit != m_GradientDirectionContainer->End(); ++gdcit)
    {
      double bValueKey = int(((m_BValue * gdcit.Value().two_norm() * gdcit.Value().two_norm())+7.5)/10)*10;
      m_BValueMap[bValueKey].push_back(gdcit.Index());
    }

  }
  if(m_BValueMap.find(0) == m_BValueMap.end())
  {
    itkExceptionMacro(<< "DiffusionMultiShellQballReconstructionImageFilter.cpp : GradientIndxMap with no b-Zero indecies found: check input BValueMap");
  }


  m_NumberOfBaselineImages = m_BValueMap[0].size();
  m_NumberOfGradientDirections = gradientDirection->Size() - m_NumberOfBaselineImages;
  // ensure that the gradient image we received has as many components as
  // the number of gradient directions
  if( gradientImage->GetVectorLength() != m_NumberOfBaselineImages + m_NumberOfGradientDirections )
  {
    itkExceptionMacro( << m_NumberOfGradientDirections << " gradients + " << m_NumberOfBaselineImages
                       << "baselines = " << m_NumberOfGradientDirections + m_NumberOfBaselineImages
                       << " directions specified but image has " << gradientImage->GetVectorLength()
                       << " components.");
  }


  ProcessObject::SetNthInput( 0, const_cast< GradientImagesType* >(gradientImage) );

  std::string gradientImageClassName(ProcessObject::GetInput(0)->GetNameOfClass());
  if ( strcmp(gradientImageClassName.c_str(),"VectorImage") != 0 )
    itkExceptionMacro( << "There is only one Gradient image. I expect that to be a VectorImage. But its of type: " << gradientImageClassName );


  m_BZeroImage = BZeroImageType::New();
  typename GradientImagesType::Pointer img = static_cast< GradientImagesType * >( ProcessObject::GetInput(0) );
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

template<class TReferenceImagePixelType, class TGradientImagePixelType, class TOdfPixelType, int NOrderL, int NrOdfDirections>
void DiffusionMultiShellQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType,NOrderL, NrOdfDirections>
::Projection1(vnl_vector<double> & vec, double delta)
{
  if (delta==0){   //Clip attenuation values. If att<0 => att=0, if att>1 => att=1
    for (unsigned int i=0; i<vec.size(); i++)
      vec[i]=(vec[i]>=0 && vec[i]<=1)*vec[i]+(vec[i]>1);
  }
  else{            //Use function from Aganj et al, MRM, 2010
    for (unsigned int i=0; i< vec.size(); i++)
      vec[i]=CalculateThreashold(vec[i], delta);
  }
}

template<class TReferenceImagePixelType, class TGradientImagePixelType, class TOdfPixelType, int NOrderL, int NrOdfDirections>
double DiffusionMultiShellQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType,NOrderL, NrOdfDirections>
::CalculateThreashold(const double value, const double delta)
{
  return (value<0)*(0.5*delta) + (value>=0 && value<delta)*(0.5*delta+0.5*(value*value)/delta)
      + (value>=delta && value<1-delta)*value+(value>=1-delta && value<1)*(1-0.5*delta-0.5*((1-value)*(1-value))/delta)
      + (value>=1)*(1-0.5*delta);
}

template<class TReferenceImagePixelType, class TGradientImagePixelType, class TOdfPixelType, int NOrderL, int NrOdfDirections>
void DiffusionMultiShellQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType,NOrderL, NrOdfDirections>
::Projection2( vnl_vector<double> & E1,vnl_vector<double> & E2, vnl_vector<double> & E3, double delta )
{
  // Each direction on its own, without the 4 matrices and 9 vectors per voxel of a vectorized
  // version; the expressions are those of the vectorized one, element by element.
  const double sF = sqrt(5.0);

  for(unsigned int i = 0 ; i < m_MaxDirections; i++)
  {
    // logarithm of the three shells' values
    const double t0 = -log(E1(i));
    const double t1 = -log(E2(i));
    const double t2 = -log(E3(i));

    // sum over the shells, and the first two shells normalized by it
    const double s0 = t0 + t1 + t2;
    const double a0 = t0 / s0;
    const double b0 = t1 / s0;

    const double ta = a0 * 3.0;
    const double tb = b0 * 3.0;
    const double e = tb - (ta * 2.0);
    const double m = (tb * 2.0) + ta;

    bool C[7];
    C[0] = (tb < 1+3*delta && 0.5+1.5*(sF+1)*delta < ta && ta < 1-3* (sF+2) *delta);
    C[1] = (e <= -1+3*(2*sF+5)*delta) && (ta>=1-3*(sF+2)*delta);
    C[2] = (m > 3-3*sF*delta) && (-1+3*(2*sF+5)*delta<e) && (e<-3*sF*delta);
    C[3] = (m >= 3-3*sF*delta && e >= -3 *sF * delta);
    C[4] = (2.5 + 1.5*(5+sF)*delta < m && m < 3-3*sF*delta && e > -3*sF*delta);
    C[5] = (ta <= 0.5+1.5 *(sF+1)*delta && m <= 2.5 + 1.5 *(5+sF) * delta);
    C[6] = !(C[0] || C[1] || C[2] || C[3] || C[4] || C[5]);

    const double A[7] = {C[0] * a0,
                         C[1] * (1.0/3.0-(sF+2)*delta),
                         C[2] * (0.2+0.8*a0-0.4*b0-delta/sF),
                         C[3] * (0.2+delta/sF),
                         C[4] * (0.2*a0+0.4*b0+2*delta/sF),
                         C[5] * (1.0/6.0+0.5*(sF+1)*delta),
                         C[6] * a0};
    const double B[7] = {C[0] * (1.0/3.0+delta),
                         C[1] * (1.0/3.0+delta),
                         C[2] * (0.4-0.4*a0+0.2*b0-2*delta/sF),
                         C[3] * (0.4-3*delta/sF),
                         C[4] * (0.4*a0+0.8*b0-delta/sF),
                         C[5] * (1.0/3.0+delta),
                         C[6] * b0};

    double a = 0;
    double b = 0;
    for(int j = 0 ; j < 7; j++)
    {
      a += A[j];
      b += B[j];
    }

    E1(i) = exp(-(a*s0));
    E2(i) = exp(-(b*s0));
    E3(i) = exp(-((1-a-b)*s0));
  }
}

template<class TReferenceImagePixelType, class TGradientImagePixelType, class TOdfPixelType, int NOrderL, int NrOdfDirections>
void DiffusionMultiShellQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType,NOrderL, NrOdfDirections>
::Projection3( vnl_vector<double> & A, vnl_vector<double> & a, vnl_vector<double> & b, double delta0)
{
  // Each direction on its own, without the matrices of 15 candidates per voxel of a vectorized
  // version; the expressions are those of the vectorized one, element by element. The candidates
  // use delta0, the constraints 0.99 delta0.
  const double s6 = sqrt(6.0);
  const double s15 = s6/2.0;
  const double delta = delta0;
  const double constraintDelta = delta0 * 0.99;

  for (unsigned int i=0; i < A.size(); i++)
  {
    const double Ai = A(i);
    const double ai = a(i);
    const double bi = b(i);

    const double AM[15] = {Ai,
                           Ai,
                           Ai,
                           delta,
                           (((Ai+ai)-bi) - (delta*s6))/3.0,
                           delta,
                           delta,
                           delta,
                           Ai,
                           ((ai*2+Ai) - delta*(2*(s6+1)))*0.2,
                           ((((bi*(-2))+Ai)+2) - delta*(2*(s6+1)))*0.2,
                           delta,
                           delta,
                           delta,
                           0.5-delta*(1+s15)};
    const double aM[15] = {ai,
                           ai,
                           -delta + 1,
                           ai,
                           ((((Ai*2)+(ai*5))+bi)+delta*s6)/6.0,
                           ai,
                           -delta + 1,
                           ((ai+bi)*0.5)+delta*(1+s15),
                           -delta + 1,
                           (((ai*4)+(Ai*2))+delta*(s6+1))*0.2,
                           -delta + 1,
                           delta*(s6+3),
                           -delta + 1,
                           -delta + 1,
                           -delta + 1};
    const double bM[15] = {bi,
                           delta,
                           bi,
                           bi,
                           ((((Ai*(-2))+ai)+(bi*5))-delta*s6)/6.0,
                           delta,
                           bi,
                           ((ai+bi)*0.5)-delta*(1+s15),
                           delta,
                           delta,
                           ((((bi*4)-(Ai*2))+1)-delta*(s6+1))*0.2,
                           delta,
                           delta,
                           ((-delta)*(s6+3)) + 1,
                           delta};

    // The nearest candidate that meets the constraints (R2 < 1e20), however far away: candidate 0
    // is the unconstrained input, for which lambda in the caller can be NaN. The vertex candidate
    // 14 meets the constraints for every input.
    unsigned int index = 0;
    double minvalue = std::numeric_limits<double>::max();
    for(int j = 0 ; j < 15 ; j++)
    {
      double R2 = 1e20;
      if (constraintDelta < AM[j] && 2*(AM[j]+constraintDelta*s15)<aM[j]-bM[j] && bM[j]>constraintDelta && aM[j]<1-constraintDelta)
        R2 = (AM[j]-Ai)*(AM[j]-Ai)+ (aM[j]-ai)*(aM[j]-ai)+(bM[j]-bi)*(bM[j]-bi);
      if(R2 < minvalue){
        minvalue = R2;
        index = j;
      }
    }

    A(i) = AM[index];
    a(i) = aM[index];
    b(i) = bM[index];
  }
}

template<class TReferenceImagePixelType, class TGradientImagePixelType, class TOdfPixelType, int NOrderL, int NrOdfDirections>
void DiffusionMultiShellQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType,NOrderL, NrOdfDirections>
::S_S0Normalization( vnl_vector<double> & vec, double S0 )
{
  for(unsigned int i = 0; i < vec.size(); i++)
  {
    if (S0==0)
      S0 = 0.01;
    vec[i] /= S0;
  }

}

template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::DoubleLogarithm(vnl_vector<double> & vec)
{
  for(unsigned int i = 0; i < vec.size(); i++)
  {
    vec[i] = log(-log(vec[i]));
  }
}

template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::AllocateOutputs()
{
  if (m_ComputeOdf)
  {
    Superclass::AllocateOutputs();
    return;
  }
  // Without ODFs the primary output gets no buffer (NODF values per voxel that nobody reads).
  OdfImageType *output = this->GetOutput();
  output->SetBufferedRegion(OutputImageRegionType());
  output->Allocate();
}

template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::BeforeThreadedGenerateData()
{
  // Only the analytical solution for three equidistant shells (Aganj 2010) is implemented; with any
  // other selection the outputs would stay unwritten.
  if (m_BValueMap.size() != 4)
    itkExceptionMacro(<< "Only three equidistant shells are supported, the selection has " << m_BValueMap.size() - 1 << ".");

  {

    BValueMapIteraotr it = m_BValueMap.begin();
    it++; // skip b0 entry
    const unsigned int bValue_shell1 = it->first;
    const unsigned int size_shell1 = it->second.size();
    IndiciesVector shell1 = it->second;
    it++;
    const unsigned int bValue_shell2 = it->first;
    const unsigned int size_shell2 = it->second.size();
    IndiciesVector shell2 = it->second;
    it++;
    const unsigned int bValue_shell3 = it->first;
    const unsigned int size_shell3 = it->second.size();
    IndiciesVector shell3 = it->second;

    // arithmetic progrssion
    if(!(bValue_shell2 - bValue_shell1 == bValue_shell1 && bValue_shell3 - bValue_shell2 == bValue_shell1))
      itkExceptionMacro(<< "Only three equidistant shells (b, 2b, 3b) are supported, the selection has b = "
                        << bValue_shell1 << ", " << bValue_shell2 << ", " << bValue_shell3 << ".");

    {
      // check if Interpolation is needed
      // if shells with different numbers of directions exist
      m_Interpolation_Flag = false;
      if(size_shell1 != size_shell2 || size_shell2 != size_shell3 || size_shell1 != size_shell3)
      {
        m_Interpolation_Flag = true;
        itkDebugMacro(<< "Shell interpolation: shells with different numbers of directions");
      }
      else
      {
        // else if each shell holds same numbers of directions, but the gradient direction differ more than one 1 degree
        m_Interpolation_Flag = CheckForDifferingShellDirections();
        if(m_Interpolation_Flag) itkDebugMacro(<< "Shell interpolation: gradient directions differ by more than 1 degree");
      }

      if(m_Interpolation_Flag)
      {
        unsigned int interp_SHOrder_shell1 = 12;
        while( ((interp_SHOrder_shell1+1)*(interp_SHOrder_shell1+2)/2) > size_shell1 && interp_SHOrder_shell1 > L )
          interp_SHOrder_shell1 -= 2 ;

        unsigned int interp_SHOrder_shell2 = 12;
        while( ((interp_SHOrder_shell2+1)*(interp_SHOrder_shell2+2)/2) > size_shell2 && interp_SHOrder_shell2 > L )
          interp_SHOrder_shell2 -= 2 ;

        unsigned int interp_SHOrder_shell3 = 12;
        while( ((interp_SHOrder_shell3+1)*(interp_SHOrder_shell3+2)/2) > size_shell3 && interp_SHOrder_shell3 > L )
          interp_SHOrder_shell3 -= 2 ;

        // Create direction container for all directions (no duplicates, different directions from all shells)
        IndiciesVector all_directions_container = mitk::gradients::GetAllUniqueDirections(m_BValueMap, m_GradientDirectionContainer.GetPointer());
        m_MaxDirections = all_directions_container.size();

        // Target SH basis per shell on all directions, and the shell's SH fit (pseudo-inverse of its
        // basis on its own directions): Values = target basis * fit * data of the shell.
        const vnl_matrix<double> Q = mitk::gradients::ComputeSphericalFromCartesian(all_directions_container, m_GradientDirectionContainer);
        const auto shellInterpolation = [&](const IndiciesVector &shell, unsigned int order, vnl_matrix<double> &targetBasis, vnl_matrix<double> &fit) {
          targetBasis = mitk::sh::ComputeSphericalHarmonicsBasis(Q, order);
          const vnl_matrix<double> shellBasis = mitk::sh::ComputeSphericalHarmonicsBasis(
            mitk::gradients::ComputeSphericalFromCartesian(shell, m_GradientDirectionContainer), order);
          fit = vnl_matrix_inverse<double>(shellBasis).inverse();
        };
        shellInterpolation(shell1, interp_SHOrder_shell1, m_TARGET_SH_shell1, m_Interpolation_SHT1_inv);
        shellInterpolation(shell2, interp_SHOrder_shell2, m_TARGET_SH_shell2, m_Interpolation_SHT2_inv);
        shellInterpolation(shell3, interp_SHOrder_shell3, m_TARGET_SH_shell3, m_Interpolation_SHT3_inv);

        ComputeReconstructionMatrix(all_directions_container);

        itkDebugMacro(<< "Shell interpolation, SH orders " << interp_SHOrder_shell1 << ", " << interp_SHOrder_shell2 << ", "
                      << interp_SHOrder_shell3 << " for " << size_shell1 << ", " << size_shell2 << ", " << size_shell3
                      << " directions; " << m_MaxDirections << " directions in all");
        return;
      }else
      {
        // Set here, not in the work units, which would all write it at the same time.
        m_MaxDirections = size_shell1;
        ComputeReconstructionMatrix(shell1);
      }
    }
  }

}


template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
  // BeforeThreadedGenerateData accepts only three equidistant shells.
  AnalyticalThreeShellReconstruction(outputRegionForThread);
}


template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::AnalyticalThreeShellReconstruction(const OutputImageRegionType& outputRegionForThread)
{
  // Input Gradient Image and Output ODF Image
  typedef typename GradientImagesType::PixelType         GradientVectorType;
  typename OdfImageType::Pointer outputImage = static_cast< OdfImageType * >(ProcessObject::GetPrimaryOutput());
  typename GradientImagesType::Pointer gradientImagePointer = static_cast< GradientImagesType * >( ProcessObject::GetInput(0) );

  // Define Image iterators; the ODF output has a buffer only when the ODFs are computed
  // (AllocateOutputs).
  std::optional<ImageRegionIterator< OdfImageType >> odfOutputImageIterator;
  if (m_ComputeOdf)
    odfOutputImageIterator.emplace(outputImage, outputRegionForThread);
  ImageRegionConstIterator< GradientImagesType > gradientInputImageIterator(gradientImagePointer, outputRegionForThread );
  ImageRegionIterator< BZeroImageType > bzeroIterator(m_BZeroImage, outputRegionForThread);
  ImageRegionIterator< CoefficientImageType > coefficientImageIterator(m_CoefficientImage, outputRegionForThread);

  // All iterators seht to Begin of the specific OutputRegion
  coefficientImageIterator.GoToBegin();
  bzeroIterator.GoToBegin();
  gradientInputImageIterator.GoToBegin();

  // Get Shell Indicies for all non-BZero Gradients
  // it MUST be a arithmetic progression eg.: 1000, 2000, 3000
  BValueMapIteraotr it = m_BValueMap.begin();
  it++;
  // it = b-value = 1000
  IndiciesVector Shell1Indiecies = it->second;
  it++;
  // it = b-value = 2000
  IndiciesVector Shell2Indiecies = it->second;
  it++;
  // it = b-value = 3000
  IndiciesVector Shell3Indiecies = it->second;
  const IndiciesVector &BZeroIndicies = m_BValueMap.at(0);

  // m_MaxDirections is set in BeforeThreadedGenerateData, with and without interpolation.

  // Nx3 Signal Matrix with E(0) = Shell 1, E(1) = Shell 2, E(2) = Shell 3
  vnl_vector< double > E1(m_MaxDirections);
  vnl_vector< double > E2(m_MaxDirections);
  vnl_vector< double > E3(m_MaxDirections);

  vnl_vector<double> AlphaValues(m_MaxDirections);
  vnl_vector<double> BetaValues(m_MaxDirections);
  vnl_vector<double> LAValues(m_MaxDirections);
  vnl_vector<double> PValues(m_MaxDirections);

  vnl_vector<double> DataShell1(Shell1Indiecies.size());
  vnl_vector<double> DataShell2(Shell2Indiecies.size());
  vnl_vector<double> DataShell3(Shell3Indiecies.size());

  // Further buffers of this work unit, so that no voxel allocates. The products are those of vnl's
  // matrix * vector operator, in the same order.
  vnl_vector<double> SignalVector(m_MaxDirections);
  vnl_vector<double> coeffs(m_CoeffReconstructionMatrix.rows());
  vnl_vector<double> odfValues(m_ODFSphericalHarmonicBasisMatrix.rows());
  const auto multiply = [](const vnl_matrix<double> &matrix, const vnl_vector<double> &vector, vnl_vector<double> &result) {
    vnl_sse<double>::matrix_x_vector(matrix.data_block(), vector.data_block(), result.data_block(), matrix.rows(), matrix.cols());
  };

  vnl_matrix<double> tempInterpolationMatrixShell1,tempInterpolationMatrixShell2,tempInterpolationMatrixShell3;

  if(m_Interpolation_Flag)
  {
    tempInterpolationMatrixShell1 = m_TARGET_SH_shell1 * m_Interpolation_SHT1_inv;
    tempInterpolationMatrixShell2 = m_TARGET_SH_shell2 * m_Interpolation_SHT2_inv;
    tempInterpolationMatrixShell3 = m_TARGET_SH_shell3 * m_Interpolation_SHT3_inv;
  }

  OdfPixelType odf(0.0);
  typename CoefficientImageType::PixelType coeffPixel(0.0);

  double P2,A,B2,B,P,alpha,beta,lambda, ER1, ER2;



  // iterate overall voxels of the gradient image region
  while( ! gradientInputImageIterator.IsAtEnd() )
  {

    odf.Fill(0.0);
    coeffPixel.Fill(0.0);

    GradientVectorType b = gradientInputImageIterator.Get();

    // The mean of all b0 volumes normalizes every shell (E = S/S0) and is the b0 output. A DWI does
    // not say which b0 volume was acquired with which shell, and their order in the file need not
    // follow the shells' order (ImageReconstruction/test_multi.dwi stores its b0 volumes before the
    // b = 3000, 2000 and 1000 blocks), so no subset of them is assigned to a shell.
    double b0average = 0;
    for (const auto index : BZeroIndicies)
      b0average += b[index];
    b0average /= BZeroIndicies.size();

    bzeroIterator.Set(b0average);
    ++bzeroIterator;

    if( (b0average != 0) && ( b0average >= m_Threshold) )
    {
      // Get the Signal-Value for each Shell at each direction (specified in the ShellIndicies Vector .. this direction corresponse to this shell...)

      for(unsigned int i = 0 ; i < Shell1Indiecies.size(); i++)
        DataShell1[i] = static_cast<double>(b[Shell1Indiecies[i]]);
      for(unsigned int i = 0 ; i < Shell2Indiecies.size(); i++)
        DataShell2[i] = static_cast<double>(b[Shell2Indiecies[i]]);
      for(unsigned int i = 0 ; i < Shell3Indiecies.size(); i++)
        DataShell3[i] = static_cast<double>(b[Shell3Indiecies[i]]);



      // Normalize the Signal: Si/S0
      S_S0Normalization(DataShell1, b0average);
      S_S0Normalization(DataShell2, b0average);
      S_S0Normalization(DataShell3, b0average);


      if(m_Interpolation_Flag)
      {
        multiply(tempInterpolationMatrixShell1, DataShell1, E1);
        multiply(tempInterpolationMatrixShell2, DataShell2, E2);
        multiply(tempInterpolationMatrixShell3, DataShell3, E3);
      }else{
        E1 = (DataShell1);
        E2 = (DataShell2);
        E3 = (DataShell3);
      }

      //Implements Eq. [19] and Fig. 4.
      Projection1(E1);
      Projection1(E2);
      Projection1(E3);
      //inqualities [31]. Taking the lograithm of th first tree inqualities
      //convert the quadratic inqualities to linear ones.
      Projection2(E1,E2,E3);

      for( unsigned int i = 0; i< m_MaxDirections; i++ )
      {
        double e1 = E1.get(i);
        double e2 = E2.get(i);
        double e3 = E3.get(i);

        P2 = e2-e1*e1;
        A = (e3 -e1*e2) / ( 2* P2);
        B2 = A * A -(e1 * e3 - e2 * e2) /P2;
        B = 0;
        if(B2 > 0) B = sqrt(B2);
        P = 0;
        if(P2 > 0) P = sqrt(P2);

        alpha = A + B;
        beta = A - B;

        PValues.put(i, P);
        AlphaValues.put(i, alpha);
        BetaValues.put(i, beta);

      }

      Projection3(PValues, AlphaValues, BetaValues);

      for(unsigned int i = 0 ; i < m_MaxDirections; i++)
      {
        const double fac = (PValues[i] * 2 ) / (AlphaValues[i] - BetaValues[i]);
        lambda = 0.5 + 0.5 * std::sqrt(1 - fac * fac);
        ER1 = std::fabs(lambda * (AlphaValues[i] - BetaValues[i]) + (BetaValues[i] - E1.get(i) ))
            + std::fabs(lambda * (AlphaValues[i] * AlphaValues[i] - BetaValues[i] * BetaValues[i]) + (BetaValues[i] * BetaValues[i] - E2.get(i) ))
            + std::fabs(lambda * (AlphaValues[i] * AlphaValues[i] * AlphaValues[i] - BetaValues[i] * BetaValues[i] * BetaValues[i]) + (BetaValues[i] * BetaValues[i] * BetaValues[i] - E3.get(i) ));
        ER2 = std::fabs((1-lambda) * (AlphaValues[i] - BetaValues[i]) + (BetaValues[i] - E1.get(i) ))
            + std::fabs((1-lambda) * (AlphaValues[i] * AlphaValues[i] - BetaValues[i] * BetaValues[i]) + (BetaValues[i] * BetaValues[i] - E2.get(i) ))
            + std::fabs((1-lambda) * (AlphaValues[i] * AlphaValues[i] * AlphaValues[i] - BetaValues[i] * BetaValues[i] * BetaValues[i]) + (BetaValues[i] * BetaValues[i] * BetaValues[i] - E3.get(i)));
        if(ER1 < ER2)
          LAValues.put(i, lambda);
        else
          LAValues.put(i, 1-lambda);

      }

      DoubleLogarithm(AlphaValues);
      DoubleLogarithm(BetaValues);

      for(unsigned int i = 0 ; i < m_MaxDirections; i++)
        SignalVector[i] = LAValues[i] * (AlphaValues[i] - BetaValues[i]) + BetaValues[i];

      multiply(m_CoeffReconstructionMatrix, SignalVector, coeffs);

      // the first coeff is a fix value
      coeffs[0] = 1.0/(2.0*sqrt(itk::Math::pi));
      for(unsigned int i = 0 ; i < coeffs.size(); i++)
        coeffPixel[i] = static_cast<TO>(coeffs[i]);

      // Cast the Signal-Type from double to float for the ODF-Image
      if (m_ComputeOdf)
      {
        multiply(m_ODFSphericalHarmonicBasisMatrix, coeffs, odfValues);
        for(unsigned int i = 0 ; i < odfValues.size(); i++)
          odf[i] = static_cast<TO>(odfValues[i]);
        odf *= ((itk::Math::pi*4)/NODF);
      }
    }

    // set ODF to ODF-Image
    coefficientImageIterator.Set(coeffPixel);
    if (odfOutputImageIterator)
    {
      odfOutputImageIterator->Set( odf );
      ++(*odfOutputImageIterator);
    }
    ++coefficientImageIterator;
    ++gradientInputImageIterator;
  }

}



template< class T, class TG, class TO, int L, int NOdfDirections>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NOdfDirections>
::ComputeReconstructionMatrix(IndiciesVector const & refVector)
{

  int numberOfGradientDirections = refVector.size();

  if( numberOfGradientDirections <= (((L+1)*(L+2))/2) || numberOfGradientDirections < 6  )
  {
    itkExceptionMacro( << "At least (L+1)(L+2)/2 gradient directions for each shell are required; current : " << numberOfGradientDirections );
  }

  CheckDuplicateDiffusionGradients();

  const int NumberOfCoeffs = (L*L + L + 2)/2 + L;
  const vnl_matrix<double> SHBasisMatrix = mitk::sh::ComputeSphericalHarmonicsBasis(
    mitk::gradients::ComputeSphericalFromCartesian(refVector, m_GradientDirectionContainer), L);

  // Per coefficient of order k and phase m: the Laplace-Beltrami penalty, weighted with the norm of
  // MITK's basis functions so that the regularization is rotation invariant; the eigenvalue
  // -k(k+1) of the Laplace-Beltrami operator; the Funk-Radon transform 2 pi P_k(0).
  vnl_matrix<double> LaplacianBaltrami(NumberOfCoeffs,NumberOfCoeffs, 0.0);
  vnl_matrix<double> FRTMatrix(NumberOfCoeffs,NumberOfCoeffs, 0.0);
  vnl_matrix<double> SHEigenvalues(NumberOfCoeffs,NumberOfCoeffs, 0.0);
  for (int k = 0; k <= L; k += 2)
  {
    for (int m = -k; m <= k; ++m)
    {
      const int j = (k * k + k + 2) / 2 + m - 1;
      LaplacianBaltrami(j, j) = mitk::LaplaceBeltramiPenalty(k, m);
      SHEigenvalues(j, j) = -k * (k + 1);
      FRTMatrix(j, j) = 2.0 * itk::Math::pi * mitk::sh::legendre0(k);
    }
  }

  const vnl_matrix<double> inverse = vnl_matrix_inverse<double>(SHBasisMatrix.transpose() * SHBasisMatrix + m_Lambda * LaplacianBaltrami).inverse();

  const double factor = (1.0/(16.0*itk::Math::pi*itk::Math::pi));
  const vnl_matrix<double> SignalReonstructionMatrix = inverse * SHBasisMatrix.transpose();
  m_CoeffReconstructionMatrix = factor * (FRTMatrix * (SHEigenvalues * SignalReonstructionMatrix));

  // SH Basis for ODF-reconstruction
  // Spherical coordinates of the ODF directions
  vnl_matrix<double> U = mitk::pointshell::Directions(NOdfDirections);
  for(int i=0; i<NOdfDirections; i++)
  {
    double cart[3];
    mitk::gradients::Cart2Sph(U(0,i),U(1,i),U(2,i),cart);
    U(0,i) = cart[0];
    U(1,i) = cart[1];
    U(2,i) = cart[2];
  }
  m_ODFSphericalHarmonicBasisMatrix = mitk::sh::ComputeSphericalHarmonicsBasis(U, L);
}

template< class T, class TG, class TO, int L, int NODF>
bool DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::CheckDuplicateDiffusionGradients()
{
  bool value = false;

  BValueMapIteraotr mapIterator = m_BValueMap.begin();
  mapIterator++;
  while(mapIterator != m_BValueMap.end())
  {
    for(auto it1 = mapIterator->second.cbegin(); it1 != mapIterator->second.cend(); ++it1)
    {
      // Every pair of the shell: it2 starts anew for each it1.
      for(auto it2 = mapIterator->second.cbegin(); it2 != mapIterator->second.cend(); ++it2)
      {
        if(m_GradientDirectionContainer->ElementAt(*it1) == m_GradientDirectionContainer->ElementAt(*it2) && it1 != it2)
        {
          itkWarningMacro( << "Some of the Diffusion Gradients equal each other. Corresponding image data should be averaged before calling this filter." );
          value = true;
        }
      }
    }
    ++mapIterator;
  }
  return value;
}

// corresponding directions between shells (e.g. dir1_shell1 vs dir1_shell2) differ more than 1 degree.
template< class T, class TG, class TO, int L, int NODF>
bool DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::CheckForDifferingShellDirections()
{
  bool interp_flag = false;

  BValueMapIteraotr mapIterator = m_BValueMap.begin();
  mapIterator++;
  IndiciesVector shell1 = mapIterator->second;
  mapIterator++;
  IndiciesVector shell2 = mapIterator->second;
  mapIterator++;
  IndiciesVector shell3 = mapIterator->second;

  for (unsigned int i=0; i< shell1.size(); i++)
    if (AbsoluteCosine(m_GradientDirectionContainer->ElementAt(shell1[i]), m_GradientDirectionContainer->ElementAt(shell2[i])) <= 0.9998) {interp_flag=true; break;}
  for (unsigned int i=0; i< shell1.size(); i++)
    if (AbsoluteCosine(m_GradientDirectionContainer->ElementAt(shell1[i]), m_GradientDirectionContainer->ElementAt(shell3[i])) <= 0.9998) {interp_flag=true; break;}
  for (unsigned int i=0; i< shell1.size(); i++)
    if (AbsoluteCosine(m_GradientDirectionContainer->ElementAt(shell2[i]), m_GradientDirectionContainer->ElementAt(shell3[i])) <= 0.9998) {interp_flag=true; break;}
  return interp_flag;
}


template< class T, class TG, class TO, int L, int NODF>
void DiffusionMultiShellQballReconstructionImageFilter<T,TG,TO,L,NODF>
::PrintSelf(std::ostream& os, Indent indent) const
{
  std::locale C("C");
  std::locale originalLocale = os.getloc();
  os.imbue(C);

  Superclass::PrintSelf(os,indent);

  //os << indent << "OdfReconstructionMatrix: " << m_ReconstructionMatrix << std::endl;
  if ( m_GradientDirectionContainer )
  {
    os << indent << "GradientDirectionContainer: "
       << m_GradientDirectionContainer << std::endl;
  }
  else
  {
    os << indent <<
          "GradientDirectionContainer: (Gradient directions not set)" << std::endl;
  }
  os << indent << "NumberOfGradientDirections: " <<
        m_NumberOfGradientDirections << std::endl;
  os << indent << "NumberOfBaselineImages: " <<
        m_NumberOfBaselineImages << std::endl;
  os << indent << "Threshold for reference B0 image: " << m_Threshold << std::endl;
  os << indent << "BValue: " << m_BValue << std::endl;

  os.imbue( originalLocale );
}




}

#endif // itkDiffusionMultiShellQballReconstructionImageFilter_cpp
