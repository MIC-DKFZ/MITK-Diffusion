/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkModelDirectionRotation.h"

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkImageWriteAccessor.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkPointShell.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <itkDiffusionTensor3D.h>

#include <vnl/algo/vnl_svd.h>

namespace
{
  using MatrixType = mitk::ModelDirectionRotation::MatrixType;

  /** SH order of the fit through which ODF images are turned; see the header. */
  constexpr unsigned int OdfShOrder = 12;

  /** Directions the SH coefficients of an SH image are fitted on; ample for order 12 (91 coefficients). */
  constexpr unsigned int ShFitDirections = 1002;

  /** Number of float values of the image: all components of all voxels and time steps. */
  std::size_t NumberOfValues(const mitk::Image *image)
  {
    if (image->GetPixelType().GetComponentType() != itk::IOComponentEnum::FLOAT)
      mitkThrow() << "Only images with float components are supported.";
    std::size_t count = image->GetPixelType().GetNumberOfComponents();
    for (unsigned int d = 0; d < image->GetDimension(); ++d)
      count *= image->GetDimension(d);
    return count;
  }

  vnl_matrix<double> ToDouble(const vnl_matrix<float> &m)
  {
    vnl_matrix<double> result(m.rows(), m.cols());
    for (unsigned int r = 0; r < m.rows(); ++r)
      for (unsigned int c = 0; c < m.cols(); ++c)
        result(r, c) = m(r, c);
    return result;
  }

  /**
   * Maps the SH coefficients c of f to those of f'(u) = f(R^T u): the least-squares fit, on the
   * given directions u, of the basis sampled at R^T u.
   */
  vnl_matrix<double> ShRotation(unsigned int order, const MatrixType &rotation, bool mrtrix, const vnl_matrix<double> &directions)
  {
    const vnl_matrix<double> turnedBack = rotation.transpose().as_matrix() * directions;
    const vnl_matrix<double> basis = ToDouble(mitk::sh::CalcShBasisForDirections(order, directions, mrtrix));
    const vnl_matrix<double> turnedBasis = ToDouble(mitk::sh::CalcShBasisForDirections(order, turnedBack, mrtrix));
    return vnl_svd<double>(basis).pinverse() * turnedBasis;
  }

  /** Applies the matrix to each pixel's vector of components (interleaved per voxel). */
  void TransformPixelVectors(mitk::Image *image, const vnl_matrix<double> &transform)
  {
    const std::size_t numberOfValues = NumberOfValues(image);
    const unsigned int components = transform.cols();
    mitk::ImageWriteAccessor accessor(image);
    auto *data = static_cast<float *>(accessor.GetData());

    vnl_vector<double> in(components);
    for (std::size_t offset = 0; offset + components <= numberOfValues; offset += components)
    {
      for (unsigned int c = 0; c < components; ++c)
        in[c] = data[offset + c];
      const vnl_vector<double> out = transform * in;
      for (unsigned int c = 0; c < transform.rows(); ++c)
        data[offset + c] = static_cast<float>(out[c]);
    }
  }

  void RotateTensors(mitk::Image *image, const MatrixType &rotation)
  {
    itk::Matrix<double, 3, 3> matrix;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        matrix[r][c] = rotation(r, c);

    const std::size_t numberOfValues = NumberOfValues(image);
    mitk::ImageWriteAccessor accessor(image);
    auto *data = static_cast<float *>(accessor.GetData());
    for (std::size_t offset = 0; offset + 6 <= numberOfValues; offset += 6)
    {
      itk::DiffusionTensor3D<float> tensor;
      for (unsigned int c = 0; c < 6; ++c)
        tensor[c] = data[offset + c];
      const auto turned = tensor.Rotate(matrix);
      for (unsigned int c = 0; c < 6; ++c)
        data[offset + c] = turned[c];
    }
  }

  void RotatePeaks(mitk::Image *image, const MatrixType &rotation)
  {
    // 4D image: the fourth index is 3 * peak + axis, so the three axes of a peak lie one volume apart.
    if (image->GetDimension() != 4 || image->GetDimension(3) % 3 != 0)
      mitkThrow() << "A peak image needs a fourth dimension of three values per peak.";
    const std::size_t volume = static_cast<std::size_t>(image->GetDimension(0)) * image->GetDimension(1) * image->GetDimension(2);
    const unsigned int numberOfPeaks = image->GetDimension(3) / 3;

    NumberOfValues(image); // throws unless the components are float, as the pointer below assumes
    mitk::ImageWriteAccessor accessor(image);
    auto *data = static_cast<float *>(accessor.GetData());
    for (unsigned int peak = 0; peak < numberOfPeaks; ++peak)
    {
      float *x = data + 3 * peak * volume;
      float *y = x + volume;
      float *z = y + volume;
      for (std::size_t voxel = 0; voxel < volume; ++voxel)
      {
        const vnl_vector_fixed<double, 3> turned = rotation * vnl_vector_fixed<double, 3>(x[voxel], y[voxel], z[voxel]);
        x[voxel] = static_cast<float>(turned[0]);
        y[voxel] = static_cast<float>(turned[1]);
        z[voxel] = static_cast<float>(turned[2]);
      }
    }
  }

  void RotateSh(mitk::ShImage *image, const MatrixType &rotation)
  {
    const bool mrtrix = image->GetShConvention() == mitk::ShImage::MRTRIX;
    TransformPixelVectors(image, ShRotation(image->ShOrder(), rotation, mrtrix, mitk::pointshell::Directions(ShFitDirections)));
  }

  void RotateOdfs(mitk::Image *image, const MatrixType &rotation)
  {
    // The SH fit of the samples f(u), pinv(B(u)) f, evaluated at the turned-back directions R^T u.
    const vnl_matrix<double> &directions = mitk::pointshell::Directions(ODF_SAMPLING_SIZE);
    const vnl_matrix<double> turnedBack = rotation.transpose().as_matrix() * directions;
    const vnl_matrix<double> basis = ToDouble(mitk::sh::CalcShBasisForDirections(OdfShOrder, directions, true));
    const vnl_matrix<double> turnedBasis = ToDouble(mitk::sh::CalcShBasisForDirections(OdfShOrder, turnedBack, true));
    TransformPixelVectors(image, turnedBasis * vnl_svd<double>(basis).pinverse());
  }
}

bool mitk::ModelDirectionRotation::HasModelDirections(const Image *image)
{
  return dynamic_cast<const TensorImage *>(image) != nullptr || dynamic_cast<const ShImage *>(image) != nullptr ||
         dynamic_cast<const OdfImage *>(image) != nullptr || dynamic_cast<const PeakImage *>(image) != nullptr;
}

void mitk::ModelDirectionRotation::Rotate(Image *image, const MatrixType &rotation)
{
  if (!HasModelDirections(image))
    mitkThrow() << "Not a tensor, SH, ODF or peak image.";

  MatrixType identity;
  identity.set_identity();
  if ((rotation * rotation.transpose() - identity).absolute_value_max() > 1e-6)
    mitkThrow() << "The matrix is not orthogonal.";
  // Without a rotation the ODF fit below would still smooth the ODFs; the image stays as it is.
  if ((rotation - identity).absolute_value_max() < 1e-9)
    return;

  if (dynamic_cast<TensorImage *>(image) != nullptr)
    RotateTensors(image, rotation);
  else if (auto *shImage = dynamic_cast<ShImage *>(image))
    RotateSh(shImage, rotation);
  else if (dynamic_cast<OdfImage *>(image) != nullptr)
    RotateOdfs(image, rotation);
  else
    RotatePeaks(image, rotation);

  // The colour maps (RgbCachedImage) and peak polydata are rebuilt from the modification time.
  image->Modified();
}
