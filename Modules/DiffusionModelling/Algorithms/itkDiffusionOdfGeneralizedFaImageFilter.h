/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

/*===================================================================

This file is based heavily on a corresponding ITK filter.

===================================================================*/
#ifndef __itkDiffusionOdfGeneralizedFaImageFilter_h_
#define __itkDiffusionOdfGeneralizedFaImageFilter_h_

//#include "MitkDiffusionImagingMBIExports.h"
#include "itkImageToImageFilter.h"
//#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

namespace itk{
/** \class DiffusionOdfGeneralizedFaImageFilter
 */

template< class TOdfPixelType,
          class TGfaPixelType,
          unsigned int NrOdfDirections>
class DiffusionOdfGeneralizedFaImageFilter :
  public ImageToImageFilter< Image< Vector< TOdfPixelType, NrOdfDirections >, 3 >,
                              Image< TGfaPixelType, 3 > >
{

public:

  enum GfaComputationMethod
  {
    GFA_STANDARD,
    GFA_QUANTILES_HIGH_LOW,
    GFA_QUANTILE_HIGH,
    GFA_MAX_ODF_VALUE,
    GFA_DECONVOLUTION_COEFFS,
    GFA_MIN_MAX_NORMALIZED_STANDARD,
    GFA_NORMALIZED_ENTROPY,
    GFA_NEMATIC_ORDER_PARAMETER,
    GFA_QUANTILE_LOW,
    GFA_MIN_ODF_VALUE,
    GFA_QUANTILES_LOW_HIGH,
    GFA_STD_BY_MAX,
    GFA_PRINCIPLE_CURVATURE,
    GFA_GENERALIZED_GFA
  };

  typedef DiffusionOdfGeneralizedFaImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< Image< Vector< TOdfPixelType, NrOdfDirections >, 3 >,
    Image< TGfaPixelType, 3 > >
                          Superclass;

   /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(DiffusionOdfGeneralizedFaImageFilter,
                                                   ImageToImageFilter);

  typedef TOdfPixelType                 OdfComponentType;

  typedef TGfaPixelType                  DirectionPixelType;

  typedef typename Superclass::InputImageType      InputImageType;

  typedef typename Superclass::OutputImageType      OutputImageType;

  typedef typename Superclass::OutputImageRegionType
                                                   OutputImageRegionType;

  itkStaticConstMacro(NOdfDirections,int,NrOdfDirections);

  void SetOdfImage( const InputImageType *image );

  itkGetMacro(ComputationMethod, GfaComputationMethod);
  itkSetMacro(ComputationMethod, GfaComputationMethod);

  itkGetMacro(Param1, double);
  itkSetMacro(Param1, double);

  itkGetMacro(Param2, double);
  itkSetMacro(Param2, double);

protected:
  DiffusionOdfGeneralizedFaImageFilter();
  ~DiffusionOdfGeneralizedFaImageFilter() override {};
  void PrintSelf(std::ostream& os, Indent indent) const override;

  void BeforeThreadedGenerateData() override;
  void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread) override;

  GfaComputationMethod m_ComputationMethod;

  double m_Param1;
  double m_Param2;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionOdfGeneralizedFaImageFilter.txx"
#endif

#endif //__itkDiffusionOdfGeneralizedFaImageFilter_h_

