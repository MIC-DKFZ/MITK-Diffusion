/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionModellingHelperFunctions.h"

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <itkShToOdfImageFilter.h>
#include <itkTensorImageToOdfImageFilter.h>
#include <mitkITKImageImport.h>
#include <mitkImageToItk.h>

namespace
{
  /** Basis matrix for directions in spherical coordinates (phi in row 0, theta in row 1). */
  template <typename T>
  vnl_matrix<T> BasisForSphericalDirections(const vnl_matrix<double> &spherical, unsigned int order, bool mrtrix)
  {
    vnl_matrix<T> basis(spherical.cols(), (order + 1) * (order + 2) / 2, 0);
    for (unsigned int i = 0; i < spherical.cols(); ++i)
    {
      for (int k = 0; k <= static_cast<int>(order); k += 2)
      {
        for (int m = -k; m <= k; ++m)
        {
          const int j = (k * k + k + 2) / 2 + m - 1;
          basis(i, j) = mitk::sh::Yj(m, k, spherical(1, i), spherical(0, i), mrtrix);
        }
      }
    }
    return basis;
  }

  template <int ShOrder>
  mitk::OdfImage::ItkOdfImageType::Pointer ShToOdf(const mitk::ShImage *shImage)
  {
    using ShConverterType = itk::ShToOdfImageFilter<float, ShOrder>;
    typename ShConverterType::InputImageType::Pointer coefficients = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(shImage, coefficients);
    auto converter = ShConverterType::New();
    converter->SetInput(coefficients);
    converter->SetToolkit(shImage->GetShConvention());
    converter->Update();
    return converter->GetOutput();
  }
}

//------------------------- SH-function ------------------------------------

double mitk::sh::legendre0(int l)
{
  if( l%2 != 0 )
  {
    return 0;
  }
  else
  {
    double prod1 = 1.0;
    for(int i=1;i<l;i+=2) prod1 *= i;
    double prod2 = 1.0;
    for(int i=2;i<=l;i+=2) prod2 *= i;
    return pow(-1.0,l/2.0)*(prod1/prod2);
  }
}

double mitk::sh::Yj(int m, int k, float theta, float phi, bool mrtrix)
{
  if (!mrtrix)
  {
    if (m<0)
      return sqrt(2.0)*static_cast<double>(::boost::math::spherical_harmonic_r(static_cast<unsigned int>(k), -m, theta, phi));
    else if (m==0)
      return static_cast<double>(::boost::math::spherical_harmonic_r(static_cast<unsigned int>(k), m, theta, phi));
    else
      return pow(-1.0,m)*sqrt(2.0)*static_cast<double>(::boost::math::spherical_harmonic_i(static_cast<unsigned int>(k), m, theta, phi));
  }

  double plm = static_cast<double>(::boost::math::legendre_p<float>(k,abs(m),-cos(theta)));
  double mag = sqrt((2.0*k+1.0)/(4.0*itk::Math::pi)*::boost::math::factorial<double>(k-abs(m))/::boost::math::factorial<double>(k+abs(m)))*plm;
  if (m>0)
    return mag*static_cast<double>(cos(m*phi));
  else if (m==0)
    return mag;
  else
    return mag*static_cast<double>(sin(-m*phi));
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromShImage(const mitk::Image *mitkImage)
{
  const auto *shImage = dynamic_cast<const mitk::ShImage*>(mitkImage);
  if (shImage == nullptr)
    mitkThrow() << "Input image is not a SH image!";
  switch (shImage->ShOrder())
  {
  case 2:
    return ShToOdf<2>(shImage);
  case 4:
    return ShToOdf<4>(shImage);
  case 6:
    return ShToOdf<6>(shImage);
  case 8:
    return ShToOdf<8>(shImage);
  case 10:
    return ShToOdf<10>(shImage);
  case 12:
    return ShToOdf<12>(shImage);
  default:
    mitkThrow() << "SH orders higher than 12 are not supported!";
  }
}

mitk::OdfImage::Pointer mitk::convert::GetOdfFromShImage(const mitk::Image *mitkImage)
{
  auto image = mitk::OdfImage::New();
  mitk::GrabItkImageMemory(GetItkOdfFromShImage(mitkImage).GetPointer(), image.GetPointer());
  return image;
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromTensorImage(const mitk::Image *mitkImage)
{
  typedef itk::TensorImageToOdfImageFilter< float, float > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( GetItkTensorFromTensorImage(mitkImage) );
  filter->Update();
  return filter->GetOutput();
}

mitk::TensorImage::ItkTensorImageType::Pointer mitk::convert::GetItkTensorFromTensorImage(const mitk::Image *mitkImage)
{
  typedef mitk::ImageToItk< mitk::TensorImage::ItkTensorImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(mitkImage);
  caster->Update();
  return caster->GetOutput();
}

mitk::OdfImage::Pointer mitk::convert::GetOdfFromTensorImage(const mitk::Image *mitkImage)
{
  auto image = mitk::OdfImage::New();
  mitk::GrabItkImageMemory(GetItkOdfFromTensorImage(mitkImage).GetPointer(), image.GetPointer());
  return image;
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromOdfImage(const mitk::Image *mitkImage)
{
  typedef mitk::ImageToItk< mitk::OdfImage::ItkOdfImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(mitkImage);
  caster->Update();
  return caster->GetOutput();
}

vnl_matrix<float> mitk::sh::CalcShBasisForDirections(unsigned int sh_order, const vnl_matrix<double> &U, bool mrtrix)
{
  vnl_matrix<double> spherical(3, U.cols());
  for (unsigned int i = 0; i < U.cols(); ++i)
  {
    double coordinates[3];
    mitk::gradients::Cart2Sph(U(0, i), U(1, i), U(2, i), coordinates);
    for (unsigned int r = 0; r < 3; ++r)
      spherical(r, i) = coordinates[r];
  }
  return BasisForSphericalDirections<float>(spherical, sh_order, mrtrix);
}

unsigned int mitk::sh::ShOrder(int num_coeffs)
{
  // An even order L has (L+1)(L+2)/2 coefficients.
  for (int order = 0; (order+1)*(order+2)/2 <= num_coeffs; order += 2)
  {
    if ((order+1)*(order+2)/2 == num_coeffs)
      return static_cast<unsigned int>(order);
  }
  mitkThrow() << num_coeffs << " is not the number of coefficients of an even SH order.";
}

float mitk::sh::GetValue(const vnl_vector<float> &coefficients, const int &sh_order, const double theta, const double phi, const bool mrtrix)
{
  if (sh_order < 0 || coefficients.size() < static_cast<unsigned int>((sh_order+1)*(sh_order+2)/2))
    mitkThrow() << coefficients.size() << " coefficients are too few for SH order " << sh_order << ".";

  float val = 0;
  for(int k=0; k<=sh_order; k+=2)
  {
    for(int m=-k; m<=k; m++)
    {
      unsigned int j = static_cast<unsigned int>((k*k + k + 2)/2 + m - 1);
      val += coefficients[j] * mitk::sh::Yj(m, k, theta, phi, mrtrix);
    }
  }

  return val;
}

float mitk::sh::GetValue(const vnl_vector<float> &coefficients, const int &sh_order, const vnl_vector_fixed<double, 3> &dir, const bool mrtrix)
{
  double spherical[3];
  mitk::gradients::Cart2Sph(dir[0], dir[1], dir[2], spherical);
  return GetValue(coefficients, sh_order, spherical[1], spherical[0], mrtrix);
}

vnl_matrix<double> mitk::sh::ComputeSphericalHarmonicsBasis(const vnl_matrix<double> & QBallReference, const unsigned int & LOrder)
{
  return BasisForSphericalDirections<double>(QBallReference, LOrder, true);
}
