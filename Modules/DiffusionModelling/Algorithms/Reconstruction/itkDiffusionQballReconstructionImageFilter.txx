/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkDiffusionQballReconstructionImageFilter_txx
#define __itkDiffusionQballReconstructionImageFilter_txx

#include "itkDiffusionQballReconstructionImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"
#include "itkPointShell.h"

#include <memory>
#include "mitkSingleShellGradients.h"

namespace itk {

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::DiffusionQballReconstructionImageFilter() :
    m_GradientDirectionContainer(nullptr),
    m_NumberOfGradientDirections(0),
    m_NumberOfEquatorSamplingPoints(0),
    m_NumberOfBaselineImages(1),
    m_Threshold(NumericTraits< ReferencePixelType >::NonpositiveMin()),
    m_BValue(1.0),
    m_DirectionsDuplicated(false),
    m_NormalizationMethod(QBR_STANDARD)
  {
    // At least 1 inputs is necessary for a vector image.
    // For images added one at a time we need at least six
    this->SetNumberOfRequiredInputs( 1 );


  }

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    void DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::BeforeThreadedGenerateData()
  {
    // Compute reconstruction matrix that is multiplied to the data-vector
    // each voxel in order to reconstruct the ODFs
    this->ComputeReconstructionMatrix();

    // Allocate the b-zero image
    m_BZeroImage = BZeroImageType::New();
    typename GradientImagesType::Pointer img = static_cast< GradientImagesType * >(
      this->ProcessObject::GetInput(0) );
    m_BZeroImage->SetSpacing( img->GetSpacing() );   // Set the image spacing
    m_BZeroImage->SetOrigin( img->GetOrigin() );     // Set the image origin
    m_BZeroImage->SetDirection( img->GetDirection() );  // Set the image direction
    m_BZeroImage->SetLargestPossibleRegion( img->GetLargestPossibleRegion());
    m_BZeroImage->SetBufferedRegion( img->GetLargestPossibleRegion() );

    m_BZeroImage->Allocate();

  }

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    typename itk::DiffusionQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType, NrOdfDirections, NrBasisFunctionCenters>::OdfPixelType
    itk::DiffusionQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType, NrOdfDirections, NrBasisFunctionCenters>::Normalize( OdfPixelType odf, double b0 )
  {
    switch( m_NormalizationMethod )
    {

      // divide by the Euclidean norm (itk::Vector::Normalize), not by the sum
    case QBR_STANDARD:
      {
        odf.Normalize();
        return odf;
        break;
      }
      // ADC style
    case QBR_B_ZERO_B_VALUE:
      {
        for(int i=0; i<NrOdfDirections; i++)
        {
          odf[i] = ((TOdfPixelType)log((TOdfPixelType)b0)-odf[i])/m_BValue;
        }
        return odf;
        break;
      }
      // divide by b-zero value of voxel
    case QBR_B_ZERO:
      {
        odf *= 1.0/b0;
        return odf;
        break;
      }
      // no normalization of ODF
    case QBR_NONE:
      {
        return odf;
        break;
      }
    }

    return odf;
  }


  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters >
    vnl_vector<TOdfPixelType> itk::DiffusionQballReconstructionImageFilter<TReferenceImagePixelType, TGradientImagePixelType, TOdfPixelType, NrOdfDirections, NrBasisFunctionCenters>::PreNormalize( vnl_vector<TOdfPixelType> vec )
  {
    switch( m_NormalizationMethod )
    {
      // standard: no normalization before reconstruction
    case QBR_STANDARD:
      {
        return vec;
        break;
      }
      // log of signal
    case QBR_B_ZERO_B_VALUE:
      {
        int n = vec.size();
        for(int i=0; i<n; i++)
        {
          vec[i] = log(vec[i]);
        }
        return vec;
        break;
      }
      // no normalization before reconstruction here
    case QBR_B_ZERO:
      {
        return vec;
        break;
      }
      // no normalization before reconstruction here
    case QBR_NONE:
      {
        return vec;
        break;
      }
    }

    return vec;
  }

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    void DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
  {
    // init output and b-zero iterators
    typename OutputImageType::Pointer outputImage =
      static_cast< OutputImageType * >(this->ProcessObject::GetPrimaryOutput());
    ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
    oit.GoToBegin();
    ImageRegionIterator< BZeroImageType > oit2(m_BZeroImage, outputRegionForThread);
    oit2.GoToBegin();
    vnl_vector<TOdfPixelType> B(m_NumberOfGradientDirections);

    // init input iterator
    typedef ImageRegionConstIterator< GradientImagesType > GradientIteratorType;
    typedef typename GradientImagesType::PixelType         GradientVectorType;
    typename GradientImagesType::Pointer gradientImagePointer = nullptr;
    // dynamic_cast would be nice, static because of SGI
    gradientImagePointer = static_cast< GradientImagesType * >(
      this->ProcessObject::GetInput(0) );
    GradientIteratorType git(gradientImagePointer, outputRegionForThread );
    git.GoToBegin();

    // set of indicies each for the baseline images and gradient images, duplicated with the shell
    const std::vector<unsigned int> &baselineind = m_BaselineIndices;
    const std::vector<unsigned int> &gradientind = m_GradientIndices;

    // Following loop does the actual reconstruction work in each voxel
    // (Tuch, Q-Ball Reconstruction [1])
    while( !git.IsAtEnd() )
    {
      // current vector of diffusion measurements
      GradientVectorType b = git.Get();

      // average of current b-zero reference values, in double so that the mean of several b0
      // volumes is not truncated (it normalizes the signal and is the b0 output)
      double b0 = 0;
      for(unsigned int i = 0; i < baselineind.size(); ++i)
      {
        b0 += b[baselineind[i]];
      }
      b0 /= this->m_NumberOfBaselineImages;

      // init resulting ODF
      OdfPixelType odf(0.0);

      // threshold on reference value to suppress noisy regions
      if( (b0 != 0) && (b0 >= m_Threshold) )
      {
        for( unsigned int i = 0; i< m_NumberOfGradientDirections; i++ )
        {
          B[i] = static_cast<TOdfPixelType>(b[gradientind[i]]);
        }

        // pre-normalization according to m_NormalizationMethod
        B = PreNormalize(B);

        // actual reconstruction
        odf = ( m_ReconstructionMatrix * B ).data_block();

        // post-normalization according to m_NormalizationMethod
        odf = Normalize(odf, b0);

      }

      for (unsigned int i=0; i<odf.Size(); i++)
          if (odf.GetElement(i)!=odf.GetElement(i))
              odf.Fill(0.0);

      // set and increment output iterators
      oit.Set( odf );
      ++oit;
      oit2.Set( static_cast<BZeroPixelType>(b0) );
      ++oit2;
      ++git; // Gradient  image iterator
    }

    std::cout << "One Thread finished reconstruction" << std::endl;
  }

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    void DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::ComputeReconstructionMatrix()
  {
    // Derived per execution: the b-value may be set after the gradients, and the half-sphere
    // duplication below doubles the count, which a second Update() must not double again. b0
    // volumes are those whose b-value rounds to 0, as in the b-value map; the model needs one shell.
    const mitk::SingleShellGradients shell = mitk::GetSingleShellGradients(m_GradientDirectionContainer, m_BValue);
    m_BaselineIndices = shell.baselineIndices;
    m_GradientIndices = shell.diffusionWeightedIndices;
    m_NumberOfBaselineImages = m_BaselineIndices.size();
    m_NumberOfGradientDirections = m_GradientIndices.size();
    m_DirectionsDuplicated = false;

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
        if (warning) break;
      }

      // handle acquisition schemes where only half of the spherical
      // shell is sampled by the gradient directions. In this case,
      // each gradient direction is duplicated in negative direction.
      vnl_vector<double> centerMass(3);
      centerMass.fill(0.0);
      int count = 0;
      for (const auto &direction : shell.directions)
      {
        centerMass += direction.as_vector();
        count ++;
      }
      centerMass /= count;
      if(centerMass.two_norm() > 0.1)
      {
        m_DirectionsDuplicated = true;
        m_NumberOfGradientDirections *= 2;
        m_GradientIndices.insert(m_GradientIndices.end(), shell.diffusionWeightedIndices.begin(), shell.diffusionWeightedIndices.end());
      }
    }

    // number of equator sampling points, derived per execution from this execution's gradient
    // count (there is no setter; a value kept from an earlier input would not fit a new one)
    this->m_NumberOfEquatorSamplingPoints
      = (int) ceil((double)sqrt(8*itk::Math::pi*this->m_NumberOfGradientDirections));

    // The matrices of this function are values (they leaked about 27 MB per execution for 64
    // gradients); smart pointers would slow down the loops below in Debug builds.
    vnl_matrix<double> Q(3, m_NumberOfGradientDirections);

    {
      // Fill matrix Q with gradient directions
      int i = 0;
      for (const auto &direction : shell.directions)
      {
        Q(0,i) = direction[0];
        Q(1,i) = direction[1];
        Q(2,i++) = direction[2];
      }
      if(m_DirectionsDuplicated)
      {
        for (const auto &direction : shell.directions)
        {
          Q(0,i) = -direction[0];
          Q(1,i) = -direction[1];
          Q(2,i++) = -direction[2];
        }
      }

    }

    // DistributePointShell returns a new matrix that the caller owns.
    std::unique_ptr<vnl_matrix_fixed<double, 3, NOdfDirections>> UOwner(
      itk::PointShell<NOdfDirections, vnl_matrix_fixed<double, 3, NOdfDirections> >::DistributePointShell());
    const vnl_matrix_fixed<double, 3, NOdfDirections> &U = *UOwner;

    std::unique_ptr<vnl_matrix_fixed<double, 3, NBasisFunctionCenters>> VOwner(
      itk::PointShell<NBasisFunctionCenters, vnl_matrix_fixed<double, 3, NBasisFunctionCenters> >::DistributePointShell());
    const vnl_matrix_fixed<double, 3, NBasisFunctionCenters> &V = *VOwner;

    // calculate sampling points on the equator perpendicular to z-axis
    vnl_matrix<double> C(3, m_NumberOfEquatorSamplingPoints);
    for(unsigned int i=0; i<m_NumberOfEquatorSamplingPoints; i++)
    {
      double theta = i * (2*itk::Math::pi / m_NumberOfEquatorSamplingPoints);
      C(0,i) = cos(theta);
      C(1,i) = sin(theta);
      C(2,i) = NumericTraits<double>::Zero;
    }

    // rotate the sampling points to each directions of interest
    vnl_matrix<double> S(3,m_NumberOfEquatorSamplingPoints*NOdfDirections);

    {
      vnl_vector_fixed<double,3> z(NumericTraits<double>::Zero);
      z.put(2,NumericTraits<double>::One);
      vnl_matrix_fixed<double,3,3> eye(NumericTraits<double>::Zero);
      eye.fill_diagonal(NumericTraits<double>::One);
      for(int i=0; i<NOdfDirections; i++)
      {
        vnl_vector_fixed<double,3> ui   = U.get_column(i);
        vnl_matrix<double> RC(3,m_NumberOfEquatorSamplingPoints);

        if( (z(0)*ui(0)+z(1)*ui(1)+z(2)*ui(2)+1) != 0 )
        {
          vnl_matrix_fixed<double,3,3> R;
          R.set_column(0, (z+ui)*(z(0)+ui(0)));
          R.set_column(1, (z+ui)*(z(1)+ui(1)));
          R.set_column(2, (z+ui)*(z(2)+ui(2)));
          R /= (z(0)*ui(0)+z(1)*ui(1)+z(2)*ui(2)+1);
          R -= eye;
          RC = R*C;
        }
        else
        {
          RC = C;
        }
        S.set_columns(i*m_NumberOfEquatorSamplingPoints, RC);
      }
    }

    // determine interpolation kernel width first
    // use to determine diffusion measurement contribution to each of the kernels
    vnl_matrix<double> H_plus(NBasisFunctionCenters,m_NumberOfGradientDirections);

    double maxSigma = itk::Math::pi/6;
    double bestSigma = maxSigma;

    {
      double stepsize = 0.01;
      double start = 0.01;
      double minCondition = NumericTraits<double>::max();

      vnl_matrix<double> H(m_NumberOfGradientDirections,NBasisFunctionCenters,(double)0);

      {

        int increasing = 0;
        for( double sigma=start; sigma<maxSigma; sigma+=stepsize)
        {
          vnl_matrix<double> tmpH(m_NumberOfGradientDirections,NBasisFunctionCenters);

          for(unsigned int r=0; r<m_NumberOfGradientDirections; r++)
          {
            for(int c=0; c<NBasisFunctionCenters; c++)
            {
              double qtv =  Q(0,r)*V(0,c)
                + Q(1,r)*V(1,c)
                + Q(2,r)*V(2,c);
              qtv = (qtv<-1.0) ? -1.0 : ( (qtv>1.0) ? 1.0 : qtv);
              double x = acos(qtv);
              tmpH(r,c) = (1.0/(sigma*sqrt(2.0*itk::Math::pi)))
                *exp((-x*x)/(2*sigma*sigma));
            }
          }

          std::unique_ptr<vnl_svd<double>> solver;
          if(m_NumberOfGradientDirections>NBasisFunctionCenters)
          {
            solver = std::make_unique<vnl_svd<double>>(tmpH);
          }
          else
          {
            solver = std::make_unique<vnl_svd<double>>(tmpH.transpose());
          }
          double condition = solver->sigma_max() / solver->sigma_min();

          std::cout << sigma << ": " << condition << std::endl;

          if( condition < minCondition )
          {
            minCondition = condition;
            bestSigma = sigma;
            H.update(tmpH);
          }
          else
          {
            // optimum assumed to be hit after condition increased 3 times
            if (++increasing>3) break;
          }
        }
      }


      auto pseudoInverse = std::make_unique<vnl_matrix_inverse<double>>(H);
      H_plus = pseudoInverse->pinverse();

      std::cout << "choosing sigma = " << bestSigma << std::endl;

    }

    // this is the contribution of each kernel to each sampling point on the
    // equator
    vnl_matrix<double> G(m_NumberOfEquatorSamplingPoints*NOdfDirections,NBasisFunctionCenters);

    {
      for(unsigned int r=0; r<m_NumberOfEquatorSamplingPoints*NOdfDirections; r++)
      {
        for(int c=0; c<NBasisFunctionCenters; c++)
        {
          double stv =  S(0,r)*V(0,c)
            + S(1,r)*V(1,c)
            + S(2,r)*V(2,c);
          stv = (stv<-1.0) ? -1.0 : ( (stv>1.0) ? 1.0 : stv);
          double x = acos(stv);
          G(r,c) = (1.0/(bestSigma*sqrt(2.0*itk::Math::pi)))
            *exp((-x*x)/(2*bestSigma*bestSigma));
        }
      }
    }

    vnl_matrix<double> GH_plus(m_NumberOfEquatorSamplingPoints*NOdfDirections,m_NumberOfGradientDirections);

    // simple matrix multiplication, manual cause of stack overflow using operator
    for (unsigned i = 0; i < m_NumberOfEquatorSamplingPoints*NOdfDirections; ++i)
    {
      for (unsigned j = 0; j < m_NumberOfGradientDirections; ++j)
      {
        double accum = G(i,0) * H_plus(0,j);
        for (unsigned k = 1; k < NOdfDirections; ++k)
          accum += G(i,k) * H_plus(k,j);
        GH_plus(i,j) = accum;
      }
    }

    typename vnl_matrix<double>::iterator it3;
    for( it3 = GH_plus.begin(); it3 != GH_plus.end(); it3++)
    {
      if(*it3<0.0)
        *it3 = 0;
    }

    // this is an addition to the original tuch algorithm
    for(unsigned int i=0; i<NOdfDirections*m_NumberOfEquatorSamplingPoints; i++)
    {
      vnl_vector< double > r = GH_plus.get_row(i);
      r /= r.sum();
      GH_plus.set_row(i,r);
    }

    m_ReconstructionMatrix = vnl_matrix<TOdfPixelType>(NOdfDirections,m_NumberOfGradientDirections,0.0);
    for(int i=0; i<NOdfDirections; i++)
    {
      for(unsigned int j=0; j<m_NumberOfGradientDirections; j++)
      {
        for(unsigned int k=0; k<m_NumberOfEquatorSamplingPoints; k++)
        {
          m_ReconstructionMatrix(i,j) += (TOdfPixelType)GH_plus(m_NumberOfEquatorSamplingPoints*i+k,j);
        }
      }
    }

    // this is also an addition to the original tuch algorithm
    for(int i=0; i<NOdfDirections; i++)
    {
      vnl_vector< TOdfPixelType > r = m_ReconstructionMatrix.get_row(i);
      r /= r.sum();
      m_ReconstructionMatrix.set_row(i,r);
    }
    std::cout << "Reconstruction Matrix computed." << std::endl;
  }

  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    void DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::SetGradientImage(const GradientDirectionContainerType *gradientDirection,
    const GradientImagesType *gradientImage )
  {
    // The gradients are classified and normalized per execution (ComputeReconstructionMatrix).
    this->m_GradientDirectionContainer = GradientDirectionContainerType::New();
    for(GradientDirectionContainerType::ConstIterator it = gradientDirection->Begin();
      it != gradientDirection->End(); it++)
    {
      this->m_GradientDirectionContainer->push_back(it.Value());
    }

    // ensure that the gradient image we received has as many components as
    // the number of gradient directions
    if( gradientImage->GetVectorLength() != gradientDirection->Size() )
    {
      itkExceptionMacro( << gradientDirection->Size() << " directions specified but image has "
        << gradientImage->GetVectorLength() << " components.");
    }

    this->ProcessObject::SetNthInput( 0,
      const_cast< GradientImagesType* >(gradientImage) );
  }


  template< class TReferenceImagePixelType,
  class TGradientImagePixelType,
  class TOdfPixelType,
    int NrOdfDirections,
    int NrBasisFunctionCenters>
    void DiffusionQballReconstructionImageFilter< TReferenceImagePixelType,
    TGradientImagePixelType, TOdfPixelType, NrOdfDirections,
    NrBasisFunctionCenters>
    ::PrintSelf(std::ostream& os, Indent indent) const
  {
    std::locale C("C");
    std::locale originalLocale = os.getloc();
    os.imbue(C);

    Superclass::PrintSelf(os,indent);

    os << indent << "OdfReconstructionMatrix: " << m_ReconstructionMatrix.rows() << " x " << m_ReconstructionMatrix.cols() << std::endl;
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

#endif // __itkDiffusionQballReconstructionImageFilter_txx
