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
#include "mitkNumericTypes.h"

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/version.hpp>
#include <boost/algorithm/string.hpp>
#include <itkPointShell.h>
#include "itkVectorContainer.h"
#include "vnl/vnl_vector.h"
#include <vtkBox.h>
#include <vtkMath.h>
#include <itksys/SystemTools.hxx>
#include <itkShToOdfImageFilter.h>
#include <itkTensorImageToOdfImageFilter.h>

//------------------------- SH-function ------------------------------------

double mitk::sh::factorial(int number) {
  if(number <= 1) return 1;
  double result = 1.0;
  for(int i=1; i<=number; i++)
    result *= i;
  return result;
}

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
  else
  {
    double plm = static_cast<double>(::boost::math::legendre_p<float>(k,abs(m),-cos(theta)));
    double mag = sqrt((2.0*k+1.0)/(4.0*itk::Math::pi)*::boost::math::factorial<double>(k-abs(m))/::boost::math::factorial<double>(k+abs(m)))*plm;
    if (m>0)
      return mag*static_cast<double>(cos(m*phi));
    else if (m==0)
      return mag;
    else
      return mag*static_cast<double>(sin(-m*phi));
  }

  return 0;
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromShImage(mitk::Image::Pointer mitkImage)
{
  mitk::ShImage::Pointer mitkShImage = dynamic_cast<mitk::ShImage*>(mitkImage.GetPointer());
  if (mitkShImage.IsNull())
    mitkThrow() << "Input image is not a SH image!";
  mitk::OdfImage::ItkOdfImageType::Pointer output;
  switch (mitkShImage->ShOrder())
  {
  case 2:
  {
    typedef itk::ShToOdfImageFilter< float, 2 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  case 4:
  {
    typedef itk::ShToOdfImageFilter< float, 4 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  case 6:
  {
    typedef itk::ShToOdfImageFilter< float, 6 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  case 8:
  {
    typedef itk::ShToOdfImageFilter< float, 8 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  case 10:
  {
    typedef itk::ShToOdfImageFilter< float, 10 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  case 12:
  {
    typedef itk::ShToOdfImageFilter< float, 12 > ShConverterType;
    typename ShConverterType::InputImageType::Pointer itkvol = ShConverterType::InputImageType::New();
    mitk::CastToItkImage(mitkImage, itkvol);
    typename ShConverterType::Pointer converter = ShConverterType::New();
    converter->SetInput(itkvol);
    converter->SetToolkit(mitkShImage->GetShConvention());
    converter->Update();
    output = converter->GetOutput();
    break;
  }
  default:
    mitkThrow() << "SH orders higher than 12 are not supported!";
  }

  return output;
}

mitk::OdfImage::Pointer mitk::convert::GetOdfFromShImage(mitk::Image::Pointer mitkImage)
{
  mitk::OdfImage::Pointer image = mitk::OdfImage::New();
  auto img = GetItkOdfFromShImage(mitkImage);
  image->InitializeByItk( img.GetPointer() );
  image->SetVolume( img->GetBufferPointer() );
  return image;
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromTensorImage(mitk::Image::Pointer mitkImage)
{
  typedef itk::TensorImageToOdfImageFilter< float, float > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( GetItkTensorFromTensorImage(mitkImage) );
  filter->Update();
  return filter->GetOutput();
}

mitk::TensorImage::ItkTensorImageType::Pointer mitk::convert::GetItkTensorFromTensorImage(mitk::Image::Pointer mitkImage)
{
  typedef mitk::ImageToItk< mitk::TensorImage::ItkTensorImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(mitkImage);
  caster->Update();
  return caster->GetOutput();
}

mitk::OdfImage::Pointer mitk::convert::GetOdfFromTensorImage(mitk::Image::Pointer mitkImage)
{
  mitk::OdfImage::Pointer image = mitk::OdfImage::New();
  auto img = GetItkOdfFromTensorImage(mitkImage);
  image->InitializeByItk( img.GetPointer() );
  image->SetVolume( img->GetBufferPointer() );
  return image;
}

mitk::OdfImage::ItkOdfImageType::Pointer mitk::convert::GetItkOdfFromOdfImage(mitk::Image::Pointer mitkImage)
{
  typedef mitk::ImageToItk< mitk::OdfImage::ItkOdfImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(mitkImage);
  caster->Update();
  return caster->GetOutput();
}

vnl_matrix<float> mitk::sh::CalcShBasisForDirections(unsigned int sh_order, vnl_matrix<double> U, bool mrtrix)
{
  vnl_matrix<float> sh_basis  = vnl_matrix<float>(U.cols(), (sh_order*sh_order + sh_order + 2)/2 + sh_order );
  for(unsigned int i=0; i<U.cols(); i++)
  {
    double x = U(0,i);
    double y = U(1,i);
    double z = U(2,i);
    double spherical[3];
    mitk::gradients::Cart2Sph(x,y,z,spherical);
    U(0,i) = spherical[0];
    U(1,i) = spherical[1];
    U(2,i) = spherical[2];
  }

  for(unsigned int i=0; i<U.cols(); i++)
  {
    for(int k=0; k<=static_cast<int>(sh_order); k+=2)
    {
      for(int m=-k; m<=k; m++)
      {
        int j = (k*k + k + 2)/2 + m - 1;
        double phi = U(0,i);
        double th = U(1,i);
        sh_basis(i,j) = mitk::sh::Yj(m,k,th,phi, mrtrix);
      }
    }
  }

  return sh_basis;
}

unsigned int mitk::sh::ShOrder(int num_coeffs)
{
  int c=3, d=2-2*num_coeffs;
  int D = c*c-4*d;
  if (D>0)
  {
    int s = (-c+static_cast<int>(sqrt(D)))/2;
    if (s<0)
      s = (-c-static_cast<int>(sqrt(D)))/2;
    return static_cast<unsigned int>(s);
  }
  else if (D==0)
    return static_cast<unsigned int>(-c/2);
  return 0;
}

float mitk::sh::GetValue(const vnl_vector<float> &coefficients, const int &sh_order, const double theta, const double phi, const bool mrtrix)
{
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

  float val = 0;
  for(int k=0; k<=sh_order; k+=2)
  {
    for(int m=-k; m<=k; m++)
    {
      int j = (k*k + k + 2)/2 + m - 1;
      val += coefficients[j] * mitk::sh::Yj(m, k, spherical[1], spherical[0], mrtrix);
    }
  }

  return val;
}



vnl_matrix<double> mitk::sh::ComputeSphericalHarmonicsBasis(const vnl_matrix<double> & QBallReference, const unsigned int & LOrder)
{
    vnl_matrix<double> SHBasisOutput(QBallReference.cols(), (LOrder+1)*(LOrder+2)*0.5);
    SHBasisOutput.fill(0.0);
    for(unsigned int i=0; i<SHBasisOutput.rows(); i++)
        for(int k = 0; k <= (int)LOrder; k += 2)
            for(int m =- k; m <= k; m++)
            {
                int j = ( k * k + k + 2 ) / 2.0 + m - 1;
                double phi = QBallReference(0,i);
                double th = QBallReference(1,i);
                double val = mitk::sh::Yj(m,k,th,phi);
                SHBasisOutput(i,j) = val;
            }
    return SHBasisOutput;
}
