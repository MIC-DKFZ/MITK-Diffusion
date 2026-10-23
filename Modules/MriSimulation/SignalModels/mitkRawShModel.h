/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_RawShModel_H
#define _MITK_RawShModel_H

#include <mitkDiffusionSignalModel.h>
#include <mitkImage.h>
#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>

namespace mitk {

/**
  * \brief The spherical harmonic representation of a prototype diffusion weighted MR signal is used to obtain the direction dependent signal.
  *
  */

template< class ScalarType = double >
class RawShModel : public DiffusionSignalModel< ScalarType >
{
public:

  RawShModel();
  ~RawShModel() override;

  DiffusionSignalModel< ScalarType >* Clone() const override { return new RawShModel(*this); }

  typedef itk::Image< double, 3 >                                         ItkDoubleImageType;
  typedef itk::Image< unsigned char, 3 >                                  ItkUcharImageType;
  typedef itk::Image< itk::DiffusionTensor3D< double >, 3 >               TensorImageType;
  typedef itk::AnalyticalDiffusionQballReconstructionImageFilter<short,short,float,2,ODF_SAMPLING_SIZE> QballFilterType;
  typedef typename DiffusionSignalModel< ScalarType >::PixelType          PixelType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientType       GradientType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientListType   GradientListType;
  typedef itk::Matrix< double, 3, 3 >                                     MatrixType;

  /** Actual signal generation **/
  PixelType SimulateMeasurement(GradientType& fiberDirection) override;
  ScalarType SimulateMeasurement(unsigned int dir, GradientType& fiberDirection) override;

  bool SetShCoefficients(vnl_vector< double > shCoefficients, double b0);
  vnl_matrix<double> SetFiberDirection(GradientType& fiberDirection);
  void SetFaRange(double min, double max){ m_FaRange.first = min; m_FaRange.second = max; }
  void SetAdcRange(double min, double max){ m_AdcRange.first = min; m_AdcRange.second = max; }
  void SetMaxNumKernels(unsigned int max){ m_MaxNumKernels = max; }
  unsigned int GetNumberOfKernels();
  std::pair< double, double > GetFaRange(){ return m_FaRange; }
  std::pair< double, double > GetAdcRange(){ return m_AdcRange; }
  unsigned int GetMaxNumKernels(){ return m_MaxNumKernels; }
  void Clear();

  std::vector< vnl_vector< double > > GetShCoefficients(){ return m_ShCoefficients; }
  std::vector< double > GetB0Signals(){ return m_B0Signal; }
  unsigned int GetShOrder(){ return m_ShOrder; }
  int GetModelIndex(){ return m_ModelIndex; }

  double GetBaselineSignal(int index){ return m_B0Signal.at(index); }
  vnl_vector< double > GetCoefficients(int listIndex){ return m_ShCoefficients.at(listIndex); }

  bool SampleKernels(Image::Pointer diffImg, ItkUcharImageType::Pointer maskImage, TensorImageType::Pointer tensorImage=nullptr, QballFilterType::CoefficientImageType::Pointer coeffImage=nullptr, ItkDoubleImageType::Pointer adcImage=nullptr);

protected:

  vnl_matrix<double> Cart2Sph( GradientListType& gradients );
  void RandomModel();

  std::vector< vnl_vector< double > > m_ShCoefficients;
  std::vector< double >               m_B0Signal;
  std::vector< GradientType >         m_PrototypeMaxDirection;
  std::pair< double, double >         m_AdcRange;
  std::pair< double, double >         m_FaRange;
  unsigned int                        m_ShOrder;
  int                                 m_ModelIndex;
  unsigned int                        m_MaxNumKernels;
};

}

#include "mitkRawShModel.cpp"

#endif

