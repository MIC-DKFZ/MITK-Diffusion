/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionFunctionCollection_h_
#define __mitkDiffusionFunctionCollection_h_


#include <MitkDiffusionCoreExports.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkVectorContainer.h>
#include <itkImage.h>
#include <itkLinearInterpolateImageFunction.h>
#include <mitkImage.h>
#include <mitkShImage.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkOdfImage.h>
#include <mitkTensorImage.h>
#include <itkOrientationDistributionFunction.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkDiffusionImageHelperFunctions.h>

namespace mitk{

class MITKDIFFUSIONCORE_EXPORT convert
{
public:
  static mitk::OdfImage::ItkOdfImageType::Pointer GetItkOdfFromTensorImage(mitk::Image::Pointer mitkImage);
  static mitk::OdfImage::Pointer GetOdfFromTensorImage(mitk::Image::Pointer mitkImage);
  static mitk::TensorImage::ItkTensorImageType::Pointer GetItkTensorFromTensorImage(mitk::Image::Pointer mitkImage);

  template<unsigned int NUM_COEFFS>
  static typename itk::Image< itk::Vector<float, NUM_COEFFS>, 3>::Pointer GetItkShFromShImage(mitk::Image::Pointer mitkImage)
  {
    mitk::ShImage::Pointer mitkShImage = dynamic_cast<mitk::ShImage*>(mitkImage.GetPointer());
    if (mitkShImage.IsNull())
      mitkThrow() << "Input image is not a SH image!";

    typename itk::Image< itk::Vector<float, NUM_COEFFS>, 3>::Pointer itkvol;
    mitk::CastToItkImage(mitkImage, itkvol);

    return itkvol;
  }

  static mitk::OdfImage::ItkOdfImageType::Pointer GetItkOdfFromShImage(mitk::Image::Pointer mitkImage);
  static mitk::OdfImage::Pointer GetOdfFromShImage(mitk::Image::Pointer mitkImage);
  static mitk::OdfImage::ItkOdfImageType::Pointer GetItkOdfFromOdfImage(mitk::Image::Pointer mitkImage);
};

class MITKDIFFUSIONCORE_EXPORT sh
{
public:
  static double factorial(int number);
  static double legendre0(int l);
  static double spherical_harmonic(int m,int l,double theta,double phi, bool complexPart);
  static double Yj(int m, int k, float theta, float phi, bool mrtrix=true);
  static vnl_matrix<float> CalcShBasisForDirections(unsigned int sh_order, vnl_matrix<double> U, bool mrtrix=true);
  static float GetValue(const vnl_vector<float>& coefficients, const int& sh_order, const vnl_vector_fixed<double, 3>& dir, const bool mrtrix);
  static float GetValue(const vnl_vector<float> &coefficients, const int &sh_order, const double theta, const double phi, const bool mrtrix);
  static unsigned int ShOrder(int num_coeffs);
  static vnl_matrix<double> ComputeSphericalHarmonicsBasis(const vnl_matrix<double> & QBallReference, const unsigned int & LOrder);

  template <unsigned int N, typename TComponent=float>
  static void SampleOdf(const vnl_vector<float>& coefficients, itk::OrientationDistributionFunction<TComponent, N>& odf)
  {
    auto dirs = odf.GetDirections();
    auto basis = CalcShBasisForDirections(ShOrder(coefficients.size()), dirs->as_matrix());
    auto odf_vals = basis * coefficients;
    for (unsigned int i=0; i<odf_vals.size(); ++i)
      odf[i] = odf_vals[i];
  }
};

}

#endif //__mitkDiffusionFunctionCollection_h_

