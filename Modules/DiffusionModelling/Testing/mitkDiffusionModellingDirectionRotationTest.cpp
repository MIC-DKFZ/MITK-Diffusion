/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkITKImageImport.h>
#include <mitkImageReadAccessor.h>
#include <mitkModelDirectionRotation.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkPointShell.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkOrientationDistributionFunction.h>

#include <vnl/vnl_rotation_matrix.h>

#include <cmath>
#include <random>
#include <vector>

class mitkDiffusionModellingDirectionRotationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingDirectionRotationTestSuite);
  MITK_TEST(Tensor_Rotation90AboutZ_TurnsXToY);
  MITK_TEST(Tensor_ObliqueRotation_IsRDRt);
  MITK_TEST(Peaks_Rotation90AboutZ_TurnEachPeak);
  MITK_TEST(Sh_ObliqueRotation_ValuesAtTurnedDirections);
  MITK_TEST(Sh_FslConvention_ValuesAtTurnedDirections);
  MITK_TEST(Odf_BandLimited_ValuesAtTurnedDirections);
  MITK_TEST(Odf_TensorOdf_PrincipalDirectionTurned);
  MITK_TEST(Odf_Identity_Unchanged);
  MITK_TEST(Rotate_NotOrthogonal_Throws);
  MITK_TEST(Rotate_PlainImage_Throws);
  CPPUNIT_TEST_SUITE_END();

  using MatrixType = mitk::ModelDirectionRotation::MatrixType;
  using Vector3 = vnl_vector_fixed<double, 3>;

  static MatrixType Rotation(Vector3 axis, double degrees)
  {
    return vnl_rotation_matrix(axis.normalize() * (degrees * M_PI / 180.0));
  }

  /** An oblique rotation, so that no direction table maps onto itself. */
  static MatrixType Oblique() { return Rotation(Vector3(1, 2, 3), 37); }

  template <class TImage, class TPixel>
  static typename TImage::Pointer OneVoxelImage(const TPixel &pixel)
  {
    using ItkImageType = itk::Image<TPixel, 3>;
    auto itkImage = ItkImageType::New();
    typename ItkImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    itkImage->FillBuffer(pixel);
    auto image = TImage::New();
    image->InitializeByItk(itkImage.GetPointer());
    image->SetVolume(itkImage->GetBufferPointer());
    return image;
  }

  static std::vector<float> Values(const mitk::Image *image, unsigned int count)
  {
    mitk::ImageReadAccessor accessor(image);
    const auto *data = static_cast<const float *>(accessor.GetData());
    return std::vector<float>(data, data + count);
  }

  static MatrixType TensorMatrix(const std::vector<float> &v)
  {
    MatrixType m;
    m(0, 0) = v[0]; m(0, 1) = v[1]; m(0, 2) = v[2];
    m(1, 0) = v[1]; m(1, 1) = v[3]; m(1, 2) = v[4];
    m(2, 0) = v[2]; m(2, 1) = v[4]; m(2, 2) = v[5];
    return m;
  }

  static mitk::TensorImage::Pointer TensorImage(const MatrixType &d)
  {
    mitk::TensorImage::PixelType tensor;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = r; c < 3; ++c)
        tensor(r, c) = d(r, c);
    return OneVoxelImage<mitk::TensorImage>(tensor);
  }

  static std::vector<float> RandomCoefficients(unsigned int count)
  {
    std::mt19937 engine(42);
    std::uniform_real_distribution<float> distribution(-1, 1);
    std::vector<float> coefficients(count);
    for (auto &c : coefficients)
      c = distribution(engine);
    return coefficients;
  }

  /** Largest difference between the turned SH function and the original at R^T u, on 1002 directions. */
  static double ShRotationError(mitk::ShImage::SH_CONVENTION convention)
  {
    constexpr unsigned int order = 8;
    constexpr unsigned int count = 45;
    const auto coefficients = RandomCoefficients(count);
    itk::Vector<float, count> pixel;
    for (unsigned int i = 0; i < count; ++i)
      pixel[i] = coefficients[i];
    auto image = OneVoxelImage<mitk::ShImage>(pixel);
    image->SetShConvention(convention);
    const bool mrtrix = convention == mitk::ShImage::MRTRIX;

    const MatrixType rotation = Oblique();
    mitk::ModelDirectionRotation::Rotate(image, rotation);
    const auto turned = Values(image, count);

    const vnl_vector<float> before(coefficients.data(), count);
    const vnl_vector<float> after(turned.data(), count);
    const auto &directions = mitk::pointshell::Directions(1002);
    double maxError = 0;
    for (unsigned int i = 0; i < directions.cols(); ++i)
    {
      const Vector3 u = directions.get_column(i);
      const double expected = mitk::sh::GetValue(before, order, rotation.transpose() * u, mrtrix);
      maxError = std::max(maxError, std::abs(mitk::sh::GetValue(after, order, u, mrtrix) - expected));
    }
    return maxError;
  }

public:
  // A tensor along x turned by 90 degrees about z lies along y.
  void Tensor_Rotation90AboutZ_TurnsXToY()
  {
    MatrixType d(0.0);
    d(0, 0) = 3;
    d(1, 1) = 1;
    d(2, 2) = 1;
    auto image = TensorImage(d);
    mitk::ModelDirectionRotation::Rotate(image, Rotation(Vector3(0, 0, 1), 90));

    MatrixType expected(0.0);
    expected(0, 0) = 1;
    expected(1, 1) = 3;
    expected(2, 2) = 1;
    CPPUNIT_ASSERT((TensorMatrix(Values(image, 6)) - expected).absolute_value_max() < 1e-6);
  }

  void Tensor_ObliqueRotation_IsRDRt()
  {
    MatrixType d;
    d(0, 0) = 1.7; d(0, 1) = 0.2; d(0, 2) = -0.3;
    d(1, 0) = 0.2; d(1, 1) = 0.6; d(1, 2) = 0.1;
    d(2, 0) = -0.3; d(2, 1) = 0.1; d(2, 2) = 0.4;
    auto image = TensorImage(d);
    const MatrixType rotation = Oblique();
    mitk::ModelDirectionRotation::Rotate(image, rotation);
    CPPUNIT_ASSERT((TensorMatrix(Values(image, 6)) - rotation * d * rotation.transpose()).absolute_value_max() < 1e-6);
  }

  // Peaks are stored as x, y, z volumes per peak; each peak turns on its own.
  void Peaks_Rotation90AboutZ_TurnEachPeak()
  {
    mitk::PeakImage::ItkPeakImageType::Pointer itkImage = mitk::PeakImage::ItkPeakImageType::New();
    mitk::PeakImage::ItkPeakImageType::RegionType region;
    region.SetSize({{2, 1, 1, 6}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    itkImage->FillBuffer(0);
    // Voxel 0: peak 0 along x (length 2), peak 1 along z; voxel 1: peak 0 along y.
    itkImage->SetPixel({{0, 0, 0, 0}}, 2);
    itkImage->SetPixel({{0, 0, 0, 5}}, 1);
    itkImage->SetPixel({{1, 0, 0, 1}}, 1);
    // As the peak reader and filters build peak images: SetVolume would fill only time step 0.
    auto image = mitk::PeakImage::New();
    mitk::GrabItkImageMemoryChannel(itkImage.GetPointer(), nullptr, image.GetPointer(), false);

    mitk::ModelDirectionRotation::Rotate(image, Rotation(Vector3(0, 0, 1), 90));

    // Layout: index = voxel + 2 * (3 * peak + axis).
    const auto v = Values(image, 12);
    const std::vector<float> expected = {0, -1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0};
    for (unsigned int i = 0; i < 12; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], v[i], 1e-6);
  }

  void Sh_ObliqueRotation_ValuesAtTurnedDirections()
  {
    const double error = ShRotationError(mitk::ShImage::MRTRIX);
    MITK_INFO << "SH rotation, MITK basis, max error " << error;
    CPPUNIT_ASSERT(error < 1e-4);
  }

  void Sh_FslConvention_ValuesAtTurnedDirections()
  {
    const double error = ShRotationError(mitk::ShImage::FSL);
    MITK_INFO << "SH rotation, FSL basis, max error " << error;
    CPPUNIT_ASSERT(error < 1e-4);
  }

  // An ODF sampled from SH of order 4 (as the Q-ball reconstructions give) turns exactly.
  void Odf_BandLimited_ValuesAtTurnedDirections()
  {
    const auto coefficients = RandomCoefficients(15);
    const vnl_vector<float> sh(coefficients.data(), 15);
    const auto &directions = mitk::pointshell::Directions(ODF_SAMPLING_SIZE);
    mitk::OdfImage::PixelType odf;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = mitk::sh::GetValue(sh, 4, Vector3(directions.get_column(i)), true);
    auto image = OneVoxelImage<mitk::OdfImage>(odf);

    const MatrixType rotation = Oblique();
    mitk::ModelDirectionRotation::Rotate(image, rotation);
    const auto turned = Values(image, ODF_SAMPLING_SIZE);

    double maxError = 0;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      const double expected = mitk::sh::GetValue(sh, 4, rotation.transpose() * Vector3(directions.get_column(i)), true);
      maxError = std::max(maxError, std::abs(turned[i] - expected));
    }
    MITK_INFO << "Band-limited ODF rotation, max error " << maxError;
    CPPUNIT_ASSERT(maxError < 1e-4);
  }

  // A tensor ODF is not band-limited; turned, its maximum lies at the turned principal direction.
  void Odf_TensorOdf_PrincipalDirectionTurned()
  {
    MatrixType d(0.0);
    d(0, 0) = 1.7e-3;
    d(1, 1) = 0.3e-3;
    d(2, 2) = 0.3e-3;
    itk::DiffusionTensor3D<float> tensor;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = r; c < 3; ++c)
        tensor(r, c) = d(r, c);
    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> odf;
    odf.InitFromTensor(tensor);
    mitk::OdfImage::PixelType pixel;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      pixel[i] = odf[i];
    auto image = OneVoxelImage<mitk::OdfImage>(pixel);

    const MatrixType rotation = Oblique();
    mitk::ModelDirectionRotation::Rotate(image, rotation);
    const auto turned = Values(image, ODF_SAMPLING_SIZE);

    // Expected: the ODF of the turned tensor, sampled at the same directions.
    itk::DiffusionTensor3D<float> turnedTensor;
    const MatrixType rdrt = rotation * d * rotation.transpose();
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = r; c < 3; ++c)
        turnedTensor(r, c) = rdrt(r, c);
    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> expected;
    expected.InitFromTensor(turnedTensor);

    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> result;
    double maxError = 0;
    double maxValue = 0;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      result[i] = turned[i];
      maxError = std::max(maxError, static_cast<double>(std::abs(turned[i] - expected[i])));
      maxValue = std::max(maxValue, static_cast<double>(expected[i]));
    }
    MITK_INFO << "Tensor ODF rotation, max error " << maxError << " of max value " << maxValue;
    CPPUNIT_ASSERT_EQUAL(expected.GetPrincipalDiffusionDirectionIndex(), result.GetPrincipalDiffusionDirectionIndex());
    CPPUNIT_ASSERT(maxError < 0.05 * maxValue);
  }

  void Odf_Identity_Unchanged()
  {
    mitk::OdfImage::PixelType pixel;
    const auto values = RandomCoefficients(ODF_SAMPLING_SIZE);
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      pixel[i] = values[i];
    auto image = OneVoxelImage<mitk::OdfImage>(pixel);
    MatrixType identity;
    identity.set_identity();
    mitk::ModelDirectionRotation::Rotate(image, identity);
    CPPUNIT_ASSERT(Values(image, ODF_SAMPLING_SIZE) == values);
  }

  void Rotate_NotOrthogonal_Throws()
  {
    auto image = TensorImage(MatrixType(0.0));
    MatrixType scaling;
    scaling.set_identity();
    scaling(0, 0) = 2;
    CPPUNIT_ASSERT_THROW(mitk::ModelDirectionRotation::Rotate(image, scaling), mitk::Exception);
  }

  void Rotate_PlainImage_Throws()
  {
    auto image = OneVoxelImage<mitk::Image>(1.0f);
    CPPUNIT_ASSERT(!mitk::ModelDirectionRotation::HasModelDirections(image));
    CPPUNIT_ASSERT_THROW(mitk::ModelDirectionRotation::Rotate(image, Oblique()), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingDirectionRotation)
