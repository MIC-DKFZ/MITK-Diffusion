/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkDiffusionModellingHelperFunctions_h
#define mitkDiffusionModellingHelperFunctions_h

#include <MitkDiffusionModellingExports.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkImageCast.h>
#include <mitkOdfImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>

namespace mitk{

class MITKDIFFUSIONMODELLING_EXPORT convert
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

  /** ODFs sampled from the SH coefficients in the image's SH convention; orders 2 to 12. */
  static mitk::OdfImage::ItkOdfImageType::Pointer GetItkOdfFromShImage(mitk::Image::Pointer mitkImage);
  static mitk::OdfImage::Pointer GetOdfFromShImage(mitk::Image::Pointer mitkImage);
  static mitk::OdfImage::ItkOdfImageType::Pointer GetItkOdfFromOdfImage(mitk::Image::Pointer mitkImage);
};

/**
 * Real spherical harmonics of even order L, with (L+1)(L+2)/2 coefficients. The convention is
 * always explicit: mrtrix = true is MITK's own basis, the one of MRtrix 0.2 (orthogonal, norm 1/2
 * for m != 0); mrtrix = false is the orthonormal FSL/Dipy basis. Values are evaluated in float
 * precision, which the reconstruction and tractography references depend on.
 */
class MITKDIFFUSIONMODELLING_EXPORT sh
{
public:
  static double legendre0(int l);
  /** Basis function of order k and phase m at polar angle theta and azimuth phi. */
  static double Yj(int m, int k, float theta, float phi, bool mrtrix);
  /** Basis matrix, one row per direction; the directions are the columns of U (Cartesian). */
  static vnl_matrix<float> CalcShBasisForDirections(unsigned int sh_order, const vnl_matrix<double> &U, bool mrtrix);
  static float GetValue(const vnl_vector<float>& coefficients, const int& sh_order, const vnl_vector_fixed<double, 3>& dir, const bool mrtrix);
  static float GetValue(const vnl_vector<float> &coefficients, const int &sh_order, const double theta, const double phi, const bool mrtrix);
  /** Order of the given number of coefficients; throws if it belongs to no even order. */
  static unsigned int ShOrder(int num_coeffs);
  /**
   * Basis matrix in MITK's convention as CalcShBasisForDirections, in double, for directions
   * given in spherical coordinates (phi in row 0, theta in row 1).
   */
  static vnl_matrix<double> ComputeSphericalHarmonicsBasis(const vnl_matrix<double> & QBallReference, const unsigned int & LOrder);
};

}

#endif
