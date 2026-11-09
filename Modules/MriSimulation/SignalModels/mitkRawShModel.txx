/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include <mitkRawShModel.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkAdcImageFilter.h>
#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkImageCast.h>
#include <mitkProperties.h>
#include <vnl/vnl_quaternion.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

template< class ScalarType >
mitk::RawShModel< ScalarType >::RawShModel()
  : m_ShOrder(0)
  , m_MaxNumKernels(1000)
{
  this->m_RandGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  this->m_RandGen->SetSeed();
  m_AdcRange.first = 0;
  m_AdcRange.second = 0.004;
  m_FaRange.first = 0;
  m_FaRange.second = 1;
}

template< class ScalarType >
mitk::RawShModel< ScalarType >::~RawShModel()
{

}

template< class ScalarType >
void mitk::RawShModel< ScalarType >::Clear()
{
  m_ShCoefficients.clear();
  m_PrototypeMaxDirection.clear();
  m_B0Signal.clear();
}

template< class ScalarType >
unsigned int mitk::RawShModel< ScalarType >::GetNumberOfKernels() const
{
  return m_B0Signal.size();
}

template< class ScalarType >
unsigned int mitk::RawShModel< ScalarType >::DrawKernel()
{
  if (m_B0Signal.empty())
    mitkThrow() << "Prototype signal model without kernels.";
  if (m_B0Signal.size()==1)
    return 0;
  return this->m_RandGen->GetIntegerVariate(m_B0Signal.size()-1);
}

template< class ScalarType >
bool mitk::RawShModel< ScalarType >::SampleKernels(Image::Pointer diffImg, ItkUcharImageType::Pointer maskImage, TensorImageType::Pointer tensorImage, QballFilterType::CoefficientImageType::Pointer itkFeatureImage, ItkDoubleImageType::Pointer adcImage)
{
  if (diffImg.IsNull())
    return false;

  typedef itk::VectorImage<short, 3> ITKDiffusionImageType;
  ITKDiffusionImageType::Pointer itkVectorImagePointer = ITKDiffusionImageType::New();
  mitk::CastToItkImage(diffImg, itkVectorImagePointer);

  const int shOrder = 2;

  if (tensorImage.IsNull())
  {
    typedef itk::DiffusionTensor3DReconstructionImageFilter< short, short, double > TensorReconstructionImageFilterType;

    // ITK takes only exact zero vectors as b0 volumes and fits every other volume at the reference
    // b-value, so it gets a copy in which the b0 volumes of MITK's rule are zero vectors, as the
    // Q-ball and ADC filters below and the b0 normalization take them.
    const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg);
    const double referenceBValue = mitk::DiffusionPropertyHelper::GetReferenceBValue(diffImg);
    auto gradientsForItk = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
    for (const auto& g : *gradients)
      gradientsForItk->push_back(mitk::DiffusionPropertyHelper::IsB0(g, referenceBValue) ? mitk::DiffusionPropertyHelper::GradientDirectionType(0.0) : g);

    TensorReconstructionImageFilterType::Pointer filter = TensorReconstructionImageFilterType::New();
    filter->SetBValue(referenceBValue);
    filter->SetGradientImage(gradientsForItk, itkVectorImagePointer );
    filter->Update();
    tensorImage = filter->GetOutput();
  }

  const int NumCoeffs = (shOrder*shOrder + shOrder + 2)/2 + shOrder;
  if (itkFeatureImage.IsNull())
  {
    QballFilterType::Pointer qballfilter = QballFilterType::New();
    qballfilter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(diffImg));
    qballfilter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg), itkVectorImagePointer );
    qballfilter->SetLambda(0.006);
    qballfilter->SetNormalizationMethod(QballFilterType::QBAR_RAW_SIGNAL);
    qballfilter->Update();
    itkFeatureImage = qballfilter->GetCoefficientImage();
  }

  if (adcImage.IsNull())
  {
    itk::AdcImageFilter< short, double >::Pointer adcFilter = itk::AdcImageFilter< short, double >::New();
    adcFilter->SetInput(itkVectorImagePointer);
    adcFilter->SetGradientDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg));
    adcFilter->SetB_value(mitk::DiffusionPropertyHelper::GetReferenceBValue(diffImg));
    adcFilter->Update();
    adcImage = adcFilter->GetOutput();
  }

  // The first b0 volume (mitk::DiffusionPropertyHelper::IsB0) normalizes the kernels.
  int b0Index = 0;
  const auto b0Indices = mitk::DiffusionPropertyHelper::GetB0Indices(mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg), mitk::DiffusionPropertyHelper::GetReferenceBValue(diffImg));
  if (!b0Indices.empty())
    b0Index = b0Indices.front();

  double max = 0;
  {
    itk::ImageRegionIterator< itk::VectorImage< short, 3 > >  it(itkVectorImagePointer, itkVectorImagePointer->GetLargestPossibleRegion());
    while(!it.IsAtEnd())
    {
      if (maskImage.IsNotNull() && maskImage->GetPixel(it.GetIndex())<=0)
      {
        ++it;
        continue;
      }
      if (it.Get()[b0Index]>max)
        max = it.Get()[b0Index];
      ++it;
    }
  }

  MITK_INFO << "Sampling signal kernels.";
  unsigned int count = 0;
  itk::ImageRegionIterator< itk::Image< itk::DiffusionTensor3D< double >, 3 > >  it(tensorImage, tensorImage->GetLargestPossibleRegion());
  while(!it.IsAtEnd())
  {
    bool skipPixel = false;
    if (maskImage.IsNotNull() && maskImage->GetPixel(it.GetIndex())<=0)
    {
      ++it;
      continue;
    }
    for (unsigned int i=0; i<mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg)->Size(); i++)
    {
      if (itkVectorImagePointer->GetPixel(it.GetIndex())[i]!=itkVectorImagePointer->GetPixel(it.GetIndex())[i] || itkVectorImagePointer->GetPixel(it.GetIndex())[i]<=0 || itkVectorImagePointer->GetPixel(it.GetIndex())[i]>itkVectorImagePointer->GetPixel(it.GetIndex())[b0Index])
      {
        skipPixel = true;
        break;
      }
    }
    if (skipPixel)
    {
      ++it;
      continue;
    }

    typedef itk::DiffusionTensor3D<double>    TensorType;
    TensorType::EigenValuesArrayType eigenvalues;
    TensorType::EigenVectorsMatrixType eigenvectors;
    TensorType tensor = it.Get();
    double FA = tensor.GetFractionalAnisotropy();
    double ADC = adcImage->GetPixel(it.GetIndex());
    QballFilterType::CoefficientImageType::PixelType itkv = itkFeatureImage->GetPixel(it.GetIndex());
    vnl_vector_fixed< double, NumCoeffs > coeffs;
    for (unsigned int c=0; c<itkv.Size(); c++)
      coeffs[c] = itkv[c]/max;

    if ( this->GetMaxNumKernels()>this->GetNumberOfKernels() && FA>m_FaRange.first && FA<m_FaRange.second && ADC>m_AdcRange.first && ADC<m_AdcRange.second)
    {
      tensor.ComputeEigenAnalysis(eigenvalues, eigenvectors);
      GradientType dir;
      dir[0] = eigenvectors(2, 0);
      dir[1] = eigenvectors(2, 1);
      dir[2] = eigenvectors(2, 2);
      this->AddKernel(coeffs.as_vector(), (double)itkVectorImagePointer->GetPixel(it.GetIndex())[b0Index]/max, dir);
      count++;
      MITK_INFO << "KERNEL: " << it.GetIndex() << " (" << count << ")";
    }
    ++it;
  }

  if (m_ShCoefficients.size()>0)
    return true;
  return false;
}

template< class ScalarType >
double mitk::RawShModel< ScalarType >::EvaluateKernel(const vnl_vector< double >& coefficients, unsigned int order, const GradientType& direction)
{
  // the analytical Q-ball filter fits the kernels with this basis and these spherical coordinates
  double spherical[3];
  mitk::gradients::Cart2Sph(direction[0], direction[1], direction[2], spherical);
  vnl_matrix<double> sphericalMatrix(3, 1);
  for (unsigned int r=0; r<3; ++r)
    sphericalMatrix(r, 0) = spherical[r];
  const vnl_matrix<double> basis = mitk::sh::ComputeSphericalHarmonicsBasis(sphericalMatrix, order);
  double value = 0;
  for (unsigned int j=0; j<basis.cols(); ++j)
    value += coefficients[j]*basis(0, j);
  return value;
}

template< class ScalarType >
void mitk::RawShModel< ScalarType >::AddKernel(const vnl_vector< double >& shCoefficients, double b0, const GradientType& direction)
{
  const unsigned int order = mitk::sh::ShOrder(shCoefficients.size());
  if (!m_ShCoefficients.empty() && order!=m_ShOrder)
    mitkThrow() << "Prototype kernels of SH order " << order << " and " << m_ShOrder << " cannot be mixed.";
  if (direction.GetNorm()<0.0001)
    mitkThrow() << "Prototype kernel without direction.";
  m_ShOrder = order;
  m_B0Signal.push_back(b0);
  m_ShCoefficients.push_back(shCoefficients);
  GradientType normalized = direction;
  normalized.Normalize();
  m_PrototypeMaxDirection.push_back(normalized);
}

template< class ScalarType >
void mitk::RawShModel< ScalarType >::AddKernel(const vnl_vector< double >& shCoefficients, double b0)
{
  // The order-2 part of the expansion is a quadratic form u^T Q u on the unit sphere; six values
  // give Q. The principal diffusion direction is where the signal is lowest: the eigenvector of the
  // smallest eigenvalue.
  const unsigned int order = std::min(2u, mitk::sh::ShOrder(shCoefficients.size()));
  auto f = [&](double x, double y, double z)
  {
    GradientType u; u[0] = x; u[1] = y; u[2] = z;
    u.Normalize();
    return EvaluateKernel(shCoefficients, order, u);
  };
  vnl_matrix<double> Q(3, 3);
  Q(0,0) = f(1,0,0); Q(1,1) = f(0,1,0); Q(2,2) = f(0,0,1);
  Q(0,1) = Q(1,0) = f(1,1,0) - (Q(0,0)+Q(1,1))/2;
  Q(0,2) = Q(2,0) = f(1,0,1) - (Q(0,0)+Q(2,2))/2;
  Q(1,2) = Q(2,1) = f(0,1,1) - (Q(1,1)+Q(2,2))/2;
  const vnl_symmetric_eigensystem<double> eigen(Q);
  GradientType direction;
  direction[0] = eigen.V(0, 0); direction[1] = eigen.V(1, 0); direction[2] = eigen.V(2, 0);
  AddKernel(shCoefficients, b0, direction);
}

template< class ScalarType >
double mitk::RawShModel< ScalarType >::KernelSignal(unsigned int kernel, const GradientType& gradient, const GradientType& fiberDirection) const
{
  if (gradient.GetNorm()<=0.001)
    return m_B0Signal.at(kernel);

  // The rotation R maps the fiber direction onto the prototype direction p, so the signal of the
  // turned kernel for gradient g is the kernel's signal for R g. The signal is antipodally
  // symmetric: the fiber is flipped towards p, which keeps the angle at most 90 degrees.
  const GradientType& p = m_PrototypeMaxDirection.at(kernel);
  GradientType f = fiberDirection;
  f.Normalize();
  if (f*p < 0)
    f *= -1;
  GradientType g = gradient;
  g.Normalize();

  GradientType axis = itk::CrossProduct(f, p);
  const double sinAngle = axis.GetNorm();
  if (sinAngle>1e-12)
  {
    axis /= sinAngle;
    const double angle = std::atan2(sinAngle, std::min(1.0, f*p));
    vnl_quaternion<double> rotation(axis.GetVnlVector(), angle);
    const vnl_vector_fixed<double, 3> rotated = rotation.rotate(g.GetVnlVector());
    g[0] = rotated[0]; g[1] = rotated[1]; g[2] = rotated[2];
  }
  return EvaluateKernel(m_ShCoefficients.at(kernel), m_ShOrder, g);
}

template< class ScalarType >
ScalarType mitk::RawShModel< ScalarType >::SimulateMeasurement(unsigned int dir, GradientType& fiberDirection)
{
  if (dir>=this->m_GradientList.size())
    return 0;
  return KernelSignal(DrawKernel(), this->m_GradientList[dir], fiberDirection);
}

template< class ScalarType >
typename mitk::RawShModel< ScalarType >::PixelType mitk::RawShModel< ScalarType >::SimulateMeasurement(GradientType& fiberDirection)
{
  const unsigned int kernel = DrawKernel();
  PixelType signal;
  signal.SetSize(this->m_GradientList.size());
  for (unsigned int p=0; p<this->m_GradientList.size(); p++)
    signal[p] = KernelSignal(kernel, this->m_GradientList[p], fiberDirection);
  return signal;
}
