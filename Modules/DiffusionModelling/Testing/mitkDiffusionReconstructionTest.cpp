/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionIndices.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionPeaks.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkException.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkPointShell.h>

#include <itkAdcImageFilter.h>
#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <itkDiffusionMultiShellQballReconstructionImageFilter.h>
#include <itkDiffusionTensorPrincipalDirectionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkOdfMaximaExtractionFilter.h>
#include <itkTensorDerivedMeasurementsFilter.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>

#include <vnl/algo/vnl_svd.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/** Tests of the reconstruction functions that the reconstruction views and command line apps use. */
class mitkDiffusionReconstructionTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionTestSuite);
  MITK_TEST(Qball_EachMethod_MatchesReference);
  MITK_TEST(Qball_SampledOdf_OnlyOnRequest);
  MITK_TEST(Qball_ShOrder_GivesCoefficientCount);
  MITK_TEST(Qball_UnsupportedShOrder_Throws);
  MITK_TEST(Qball_NoDwi_Throws);
  MITK_TEST(Qball_SingleShellMethods_RejectMultiShell);
  MITK_TEST(Qball_MultiShell_SameAsFilter);
  MITK_TEST(Qball_MultiShell_ShellSelection);
  MITK_TEST(Tensors_LinearLeastSquares_MatchesReference);
  MITK_TEST(Tensors_EigenvalueCorrection_SameAsFilter);
  MITK_TEST(Tensors_ZeroNegativeEigenvalues);
  MITK_TEST(Tensors_NoDwi_Throws);
  MITK_TEST(PointShellGradients_TableAndB0);
  MITK_TEST(TensorsToDwi_WorldGradients);
  MITK_TEST(Peaks_FslConvention_SameAsMrtrix);
  MITK_TEST(Peaks_EachAxisOnce);
  MITK_TEST(Peaks_SameAsFilter_AllVolumes);
  MITK_TEST(Peaks_MaskOfOtherSize_Throws);
  MITK_TEST(TensorPeaks_SameAsFilter);
  MITK_TEST(TensorIndices_MatchReferences);
  MITK_TEST(Gfa_MatchesReference);
  MITK_TEST(Adc_SameAsFilter);
  MITK_TEST(Indices_Scale);
  MITK_TEST(Indices_WrongInput_Throws);
  CPPUNIT_TEST_SUITE_END();

  using GradientContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using DwiType = mitk::DiffusionPropertyHelper::ImageType;
  using Direction = vnl_vector_fixed<double, 3>;
  using TensorImageType = mitk::TensorImage::ItkTensorImageType;

  mitk::Image::Pointer m_Dwi;

  /** n unit directions spread over the sphere (Fibonacci lattice). */
  static std::vector<Direction> Directions(unsigned int n)
  {
    std::vector<Direction> directions;
    const double golden = itk::Math::pi * (3.0 - std::sqrt(5.0));
    for (unsigned int i = 0; i < n; ++i)
    {
      const double z = 1.0 - 2.0 * (i + 0.5) / n;
      const double r = std::sqrt(1.0 - z * z);
      directions.emplace_back(r * std::cos(golden * i), r * std::sin(golden * i), z);
    }
    return directions;
  }

  /**
   * A 3x1x1 DWI with two b0 volumes and 30 directions per b-value (b = referenceBValue * |g|^2),
   * as an MITK image with world gradients; the signal follows cylindrical tensors.
   */
  static mitk::Image::Pointer SimulateDwi(const std::vector<double> &bValues, double referenceBValue)
  {
    auto gradients = GradientContainerType::New();
    gradients->push_back(Direction(0.0));
    gradients->push_back(Direction(0.0));
    for (double b : bValues)
      for (const auto &d : Directions(30))
        gradients->push_back(d * std::sqrt(b / referenceBValue));

    auto dwi = DwiType::New();
    DwiType::RegionType region;
    region.SetSize({{3, 1, 1}});
    dwi->SetRegions(region);
    dwi->SetVectorLength(gradients->Size());
    dwi->Allocate();
    for (unsigned int x = 0; x < 3; ++x)
    {
      const Direction e(std::cos(0.5 * x), std::sin(0.5 * x), 0.2);
      const Direction u = e / e.two_norm();
      DwiType::PixelType pixel(gradients->Size());
      for (unsigned int i = 0; i < gradients->Size(); ++i)
      {
        const auto &g = gradients->ElementAt(i);
        const double gu = dot_product(g, u);
        const double gDg = 1.4e-3 * gu * gu + 0.3e-3 * dot_product(g, g);
        pixel[i] = static_cast<short>(std::round(1000 * std::exp(-referenceBValue * gDg)));
      }
      dwi->SetPixel({{x, 0, 0}}, pixel);
    }

    auto image = mitk::GrabItkImageMemory(dwi.GetPointer());
    mitk::DiffusionPropertyHelper::MeasurementFrameType identity;
    identity.set_identity();
    mitk::DiffusionPropertyHelper::InitializeImage(image, gradients.GetPointer(), referenceBValue, identity, false);
    return image;
  }

  static mitk::QballParameters Parameters(mitk::QballMethod method)
  {
    mitk::QballParameters parameters;
    parameters.method = method;
    parameters.sampledOdf = true;
    return parameters;
  }

  template <class TImage>
  static bool SameVoxels(const TImage *expected, const TImage *actual)
  {
    itk::ImageRegionConstIterator<TImage> itExpected(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TImage> itActual(actual, actual->GetLargestPossibleRegion());
    for (; !itExpected.IsAtEnd(); ++itExpected, ++itActual)
    {
      if (itExpected.Get() != itActual.Get())
        return false;
    }
    return expected->GetLargestPossibleRegion() == actual->GetLargestPossibleRegion();
  }

  /** Largest difference of two ODF images, each relative to the largest absolute value of the reference ODF in its voxel. */
  static double LargestRelativeOdfDifference(mitk::OdfImage *actual, mitk::OdfImage *reference)
  {
    auto a = mitk::convert::GetItkOdfFromOdfImage(actual);
    auto b = mitk::convert::GetItkOdfFromOdfImage(reference);
    double largest = 0;
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> itA(a, a->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> itB(b, b->GetLargestPossibleRegion());
    for (; !itB.IsAtEnd(); ++itA, ++itB)
    {
      double scale = 0;
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
        scale = std::max(scale, std::abs(static_cast<double>(itB.Get()[i])));
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE && scale > 0; ++i)
        largest = std::max(largest, std::abs(static_cast<double>(itA.Get()[i]) - itB.Get()[i]) / scale);
    }
    return largest;
  }

public:
  void setUp() override
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"), &functor);
  }

  void tearDown() override { m_Dwi = nullptr; }

  // The view's methods give the references of mitkImageReconstructionTest (lambda 0.006, order 4,
  // same tolerances).
  void Qball_EachMethod_MatchesReference()
  {
    const std::vector<std::tuple<mitk::QballMethod, std::string, double>> cases = {
      {mitk::QballMethod::Numerical, "test_QN0.qbi", 1e-4},
      {mitk::QballMethod::Standard, "test_QA0_v2.qbi", 1e-4},
      {mitk::QballMethod::AdcProfile, "test_QA4_v2.qbi", 1e-4},
      {mitk::QballMethod::RawSignal, "test_QA5_v2.qbi", 0.1},
      {mitk::QballMethod::SolidAngle, "test_QA6_v2.qbi", 1e-4}};
    for (const auto &[method, file, tolerance] : cases)
    {
      const auto result = mitk::ReconstructQball(m_Dwi, Parameters(method));
      auto reference = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/" + file));
      CPPUNIT_ASSERT_MESSAGE(file, result.odf.IsNotNull());
      CPPUNIT_ASSERT_MESSAGE(file, mitk::Equal(*result.odf, *reference, tolerance, true));
      // The absolute tolerances are percents of an ODF value (about 1/252); the ODFs agree to float rounding.
      CPPUNIT_ASSERT_MESSAGE(file, LargestRelativeOdfDifference(result.odf, reference) < 1e-5);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(file, method != mitk::QballMethod::Numerical, result.sh.IsNotNull());
    }
  }

  void Qball_SampledOdf_OnlyOnRequest()
  {
    auto parameters = Parameters(mitk::QballMethod::Standard);
    parameters.sampledOdf = false;
    const auto result = mitk::ReconstructQball(m_Dwi, parameters);
    CPPUNIT_ASSERT(result.sh.IsNotNull());
    CPPUNIT_ASSERT(result.odf.IsNull());

    // The SH image equals the filter's coefficient image.
    using FilterType = itk::AnalyticalDiffusionQballReconstructionImageFilter<short, short, float, 4, ODF_SAMPLING_SIZE>;
    DwiType::Pointer itkDwi;
    mitk::CastToItkImage(m_Dwi, itkDwi);
    auto filter = FilterType::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(m_Dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(m_Dwi), itkDwi);
    filter->SetLambda(0.006);
    filter->Update();
    FilterType::CoefficientImageType::Pointer coefficients;
    mitk::CastToItkImage(result.sh.GetPointer(), coefficients);
    CPPUNIT_ASSERT(SameVoxels(filter->GetCoefficientImage().GetPointer(), coefficients.GetPointer()));
    CPPUNIT_ASSERT_EQUAL(4u, result.sh->ShOrder());
  }

  void Qball_ShOrder_GivesCoefficientCount()
  {
    for (const auto &[order, count] : std::vector<std::pair<unsigned int, std::size_t>>{{2, 6}, {4, 15}, {6, 28}})
    {
      auto parameters = Parameters(mitk::QballMethod::SolidAngle);
      parameters.shOrder = order;
      const auto result = mitk::ReconstructQball(m_Dwi, parameters);
      CPPUNIT_ASSERT_EQUAL(count, result.sh->GetPixelType().GetNumberOfComponents());
    }
  }

  void Qball_UnsupportedShOrder_Throws()
  {
    for (unsigned int order : {0u, 3u, 14u})
    {
      auto parameters = Parameters(mitk::QballMethod::Standard);
      parameters.shOrder = order;
      CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(m_Dwi, parameters), mitk::Exception);
    }
  }

  void Qball_NoDwi_Throws()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(tensors, Parameters(mitk::QballMethod::Standard)), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(nullptr, Parameters(mitk::QballMethod::Standard)), mitk::Exception);
  }

  // The single-shell models would fit every volume with the reference b-value.
  void Qball_SingleShellMethods_RejectMultiShell()
  {
    auto dwi = SimulateDwi({1000, 2000, 3000}, 3000);
    for (auto method : {mitk::QballMethod::Numerical, mitk::QballMethod::Standard, mitk::QballMethod::SolidAngle,
                        mitk::QballMethod::AdcProfile, mitk::QballMethod::RawSignal})
      CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(dwi, Parameters(method)), itk::ExceptionObject);
  }

  void Qball_MultiShell_SameAsFilter()
  {
    auto dwi = SimulateDwi({1000, 2000, 3000}, 3000);
    const auto result = mitk::ReconstructQball(dwi, Parameters(mitk::QballMethod::MultiShell));

    using FilterType = itk::DiffusionMultiShellQballReconstructionImageFilter<short, short, float, 4, ODF_SAMPLING_SIZE>;
    DwiType::Pointer itkDwi;
    mitk::CastToItkImage(dwi, itkDwi);
    auto filter = FilterType::New();
    filter->SetBValueMap(mitk::DiffusionPropertyHelper::GetBValueMap(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi, 3000);
    filter->SetLambda(0.006);
    filter->Update();

    FilterType::CoefficientImageType::Pointer coefficients;
    mitk::CastToItkImage(result.sh.GetPointer(), coefficients);
    CPPUNIT_ASSERT(SameVoxels(filter->GetCoefficientImage().GetPointer(), coefficients.GetPointer()));
    mitk::OdfImage::ItkOdfImageType::Pointer odfs;
    mitk::CastToItkImage(result.odf.GetPointer(), odfs);
    CPPUNIT_ASSERT(SameVoxels(filter->GetOutput(), odfs.GetPointer()));
  }

  void Qball_MultiShell_ShellSelection()
  {
    auto dwi = SimulateDwi({1000, 2000, 3000, 4000}, 4000);
    auto parameters = Parameters(mitk::QballMethod::MultiShell);

    // All four shells, or three that are not b, 2b, 3b: no analytical solution.
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(dwi, parameters), itk::ExceptionObject);
    parameters.shells = {1000, 2000, 4000};
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(dwi, parameters), itk::ExceptionObject);
    parameters.shells = {2000, 3000, 4000};
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(dwi, parameters), itk::ExceptionObject);

    // A shell the DWI does not have.
    parameters.shells = {1000, 2000, 5000};
    CPPUNIT_ASSERT_THROW(mitk::ReconstructQball(dwi, parameters), mitk::Exception);

    parameters.shells = {1000, 2000, 3000};
    const auto result = mitk::ReconstructQball(dwi, parameters);
    CPPUNIT_ASSERT(result.sh.IsNotNull());
    CPPUNIT_ASSERT(result.odf.IsNotNull());
  }

  // The view's linear least-squares tensors of a single-shell DWI give the reference of
  // mitkImageReconstructionTest.
  void Tensors_LinearLeastSquares_MatchesReference()
  {
    auto tensors = mitk::ReconstructTensors(m_Dwi, mitk::TensorFit::LinearLeastSquares, 0, false);
    auto reference = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT(mitk::Equal(*tensors, *reference, 1e-4, true));
  }

  static mitk::Image::Pointer LoadDwi(const std::string &path)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(path), &functor);
  }

  static TensorImageType::Pointer ItkTensors(const mitk::TensorImage *image)
  {
    TensorImageType::Pointer tensors;
    mitk::CastToItkImage(image, tensors);
    return tensors;
  }

  static bool HasNegativeEigenvalue(const TensorImageType::PixelType &tensor)
  {
    TensorImageType::PixelType::EigenValuesArrayType eigenvalues;
    tensor.ComputeEigenValues(eigenvalues);
    return eigenvalues[0] < 0 || eigenvalues[1] < 0 || eigenvalues[2] < 0;
  }

  // diffusionImage.hdwi has voxels whose eigenvalues the correction cannot make positive (see
  // mitkDiffusionModellingAlgorithmsTest); without zeroing, the function gives the filter's tensors.
  void Tensors_EigenvalueCorrection_SameAsFilter()
  {
    auto dwi = LoadDwi("diffusionImage.hdwi");
    DwiType::Pointer itkDwi;
    mitk::CastToItkImage(dwi, itkDwi);
    auto filter = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi);
    filter->SetB0Threshold(0);
    filter->Update();

    auto tensors = mitk::ReconstructTensors(dwi, mitk::TensorFit::EigenvalueCorrection, 0, false);
    CPPUNIT_ASSERT(SameVoxels(filter->GetOutput(), ItkTensors(tensors).GetPointer()));
  }

  // Zeroing sets exactly the tensors with a negative eigenvalue to zero, for both fits.
  void Tensors_ZeroNegativeEigenvalues()
  {
    auto dwi = LoadDwi("diffusionImage.hdwi");
    for (auto fit : {mitk::TensorFit::LinearLeastSquares, mitk::TensorFit::EigenvalueCorrection})
    {
      auto raw = ItkTensors(mitk::ReconstructTensors(dwi, fit, 0, false));
      auto zeroed = ItkTensors(mitk::ReconstructTensors(dwi, fit, 0, true));
      unsigned int negatives = 0;
      itk::ImageRegionConstIterator<TensorImageType> itRaw(raw, raw->GetLargestPossibleRegion());
      itk::ImageRegionConstIterator<TensorImageType> itZeroed(zeroed, zeroed->GetLargestPossibleRegion());
      for (; !itRaw.IsAtEnd(); ++itRaw, ++itZeroed)
      {
        if (HasNegativeEigenvalue(itRaw.Get()))
        {
          ++negatives;
          CPPUNIT_ASSERT(itZeroed.Get() == TensorImageType::PixelType(0.0));
        }
        else
        {
          CPPUNIT_ASSERT(itZeroed.Get() == itRaw.Get());
        }
      }
      // The linear fit of this DWI has tensors with negative eigenvalues, so the zeroing is tested.
      if (fit == mitk::TensorFit::LinearLeastSquares)
        CPPUNIT_ASSERT(negatives > 0);
    }
  }

  void Tensors_NoDwi_Throws()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::ReconstructTensors(tensors, mitk::TensorFit::LinearLeastSquares, 0, false), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::ReconstructTensors(nullptr, mitk::TensorFit::EigenvalueCorrection, 0, false), mitk::Exception);
  }

  void PointShellGradients_TableAndB0()
  {
    for (unsigned int n : {12u, 92u, 252u})
    {
      const auto gradients = mitk::PointShellGradients(n);
      const auto &table = mitk::pointshell::Directions(n);
      CPPUNIT_ASSERT_EQUAL(n + 1, static_cast<unsigned int>(gradients->Size()));
      for (unsigned int i = 0; i < n; ++i)
        CPPUNIT_ASSERT(gradients->ElementAt(i) == Direction(table(0, i), table(1, i), table(2, i)));
      CPPUNIT_ASSERT(gradients->ElementAt(n) == Direction(0.0));
    }
    CPPUNIT_ASSERT_THROW(mitk::PointShellGradients(100), mitk::Exception);
  }

  // The simulated DWI carries world gradients: identity frame, image rotation off, the modified
  // gradients equal to the table.
  void TensorsToDwi_WorldGradients()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    auto dwi = mitk::TensorsToDwi(tensors, 92, 1000);

    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    CPPUNIT_ASSERT_EQUAL(93u, static_cast<unsigned int>(dwi->GetPixelType().GetNumberOfComponents()));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi), 0.0);
    // The flag behind DiffusionPropertyHelper::SetApplyMatrixToGradients.
    bool applyRotation = true;
    CPPUNIT_ASSERT(dwi->GetPropertyList()->GetBoolProperty("DWMRI.ApplyMatrixToGradients", applyRotation));
    CPPUNIT_ASSERT(!applyRotation);
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::GetMeasurementFrame(dwi).is_identity());

    const auto expected = mitk::PointShellGradients(92);
    const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(dwi);
    for (unsigned int i = 0; i < expected->Size(); ++i)
      CPPUNIT_ASSERT((gradients->ElementAt(i) - expected->ElementAt(i)).two_norm() < 1e-12);

    CPPUNIT_ASSERT_THROW(mitk::TensorsToDwi(nullptr, 92, 1000), mitk::Exception);
  }

  using ShItkImageType = itk::Image<itk::Vector<float, 15>, 3>;
  using PeakItkImageType = mitk::PeakImage::ItkPeakImageType;

  static mitk::ShImage::Pointer LoadShImage()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &functor);
    shImage->SetShConvention(mitk::ShImage::MRTRIX);
    return shImage;
  }

  static vnl_matrix<double> Basis(bool mrtrix)
  {
    const vnl_matrix<float> basis = mitk::sh::CalcShBasisForDirections(4, mitk::pointshell::Directions(1002), mrtrix);
    vnl_matrix<double> result(basis.rows(), basis.cols());
    for (unsigned int r = 0; r < basis.rows(); ++r)
      for (unsigned int c = 0; c < basis.cols(); ++c)
        result(r, c) = basis(r, c);
    return result;
  }

  /**
   * The same ODFs in the FSL convention: the coefficients refitted by least squares on the 1002
   * directions, which determine order-4 functions.
   */
  static mitk::ShImage::Pointer ToFsl(const mitk::ShImage *mrtrixImage)
  {
    const vnl_matrix<double> conversion = vnl_svd<double>(Basis(false)).pinverse() * Basis(true);
    auto fslImage = mrtrixImage->Clone();
    auto coefficients = mitk::ImageToItkImage<itk::Vector<float, 15>, 3>(fslImage.GetPointer());
    for (itk::ImageRegionIterator<ShItkImageType> it(coefficients, coefficients->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      vnl_vector<double> c(15);
      for (unsigned int i = 0; i < 15; ++i)
        c[i] = it.Get()[i];
      const vnl_vector<double> converted = conversion * c;
      ShItkImageType::PixelType pixel;
      for (unsigned int i = 0; i < 15; ++i)
        pixel[i] = static_cast<float>(converted[i]);
      it.Set(pixel);
    }
    fslImage->SetShConvention(mitk::ShImage::FSL);
    return fslImage;
  }

  /** The non-zero peaks of a voxel. */
  static std::vector<Direction> Peaks(const PeakItkImageType *image, const itk::Index<3> &voxel)
  {
    std::vector<Direction> peaks;
    const unsigned int count = image->GetLargestPossibleRegion().GetSize(3) / 3;
    for (unsigned int p = 0; p < count; ++p)
    {
      Direction peak;
      for (unsigned int c = 0; c < 3; ++c)
        peak[c] = image->GetPixel({{voxel[0], voxel[1], voxel[2], 3 * p + c}});
      if (peak.two_norm() > 0)
        peaks.push_back(peak);
    }
    return peaks;
  }

  /**
   * Largest distance between the peaks of two peak images, independent of the order of the peaks
   * in a voxel (nearly equal peaks may swap) and of their signs; infinite where the counts differ.
   */
  static double PeakDifference(mitk::PeakImage *a, mitk::PeakImage *b)
  {
    auto peaksA = a->GetItkPeakImage();
    auto peaksB = b->GetItkPeakImage();
    CPPUNIT_ASSERT(peaksA->GetLargestPossibleRegion() == peaksB->GetLargestPossibleRegion());
    const auto size = peaksA->GetLargestPossibleRegion().GetSize();
    double largest = 0;
    for (itk::IndexValueType z = 0; z < static_cast<itk::IndexValueType>(size[2]); ++z)
      for (itk::IndexValueType y = 0; y < static_cast<itk::IndexValueType>(size[1]); ++y)
        for (itk::IndexValueType x = 0; x < static_cast<itk::IndexValueType>(size[0]); ++x)
        {
          const auto voxelA = Peaks(peaksA, {{x, y, z}});
          const auto voxelB = Peaks(peaksB, {{x, y, z}});
          if (voxelA.size() != voxelB.size())
            return std::numeric_limits<double>::infinity();
          for (const auto &peak : voxelA)
          {
            double nearest = std::numeric_limits<double>::infinity();
            for (const auto &other : voxelB)
              nearest = std::min({nearest, (peak - other).two_norm(), (peak + other).two_norm()});
            largest = std::max(largest, nearest);
          }
        }
    return largest;
  }

  /**
   * Mask of the voxels whose ODF varies by more than a tenth of its maximum: the isotropic
   * background of sh_image_test.nii.gz has no peaks, and the broad maxima of weakly anisotropic
   * ODFs move with the last bits of the coefficients.
   */
  static mitk::Image::Pointer AnisotropyMask(const mitk::ShImage *shImage, unsigned int &count)
  {
    auto odfs = mitk::convert::GetItkOdfFromShImage(shImage);
    auto mask = itk::Image<unsigned char, 3>::New();
    mask->CopyInformation(odfs);
    mask->SetRegions(odfs->GetLargestPossibleRegion());
    mask->Allocate();
    count = 0;
    itk::ImageRegionIterator<itk::Image<unsigned char, 3>> itMask(mask, mask->GetLargestPossibleRegion());
    for (itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> it(odfs, odfs->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it, ++itMask)
    {
      const auto odf = it.Get();
      const auto [min, max] = std::minmax_element(odf.Begin(), odf.End());
      const bool anisotropic = *max - *min > 0.1 * *max;
      itMask.Set(anisotropic ? 1 : 0);
      count += anisotropic ? 1 : 0;
    }
    return mitk::GrabItkImageMemory(mask.GetPointer());
  }

  // The peaks follow the SH convention stored in the image (PL-M5): the same ODFs in the FSL
  // convention give the same peaks, the FSL coefficients read as MRtrix ones do not (0.55 apart).
  // The gradient descent that refines each peak (f tolerance 1e-6) stops up to 0.05 (3 degrees)
  // apart on the broad lobes of this image when the coefficients differ in the last bits.
  void Peaks_FslConvention_SameAsMrtrix()
  {
    auto mrtrixImage = LoadShImage();
    auto fslImage = ToFsl(mrtrixImage);
    unsigned int anisotropic = 0;
    auto mask = AnisotropyMask(mrtrixImage, anisotropic);
    CPPUNIT_ASSERT(anisotropic > 0);
    mitk::PeakParameters parameters;

    const auto expected = mitk::ExtractPeaks(mrtrixImage, parameters, mask);
    const auto actual = mitk::ExtractPeaks(fslImage, parameters, mask);
    CPPUNIT_ASSERT(PeakDifference(expected.peaks, actual.peaks) < 0.1);

    fslImage->SetShConvention(mitk::ShImage::MRTRIX);
    const auto misread = mitk::ExtractPeaks(fslImage, parameters, mask);
    CPPUNIT_ASSERT(PeakDifference(expected.peaks, misread.peaks) > 0.1);
  }

  // Of two peak candidates closer than the angular threshold only the larger is kept, also when
  // both have exactly the same length, as the two opposite (antipodal) candidates of one axis in
  // voxel (6, 4, 1) in MITK's convention. Isotropic voxels (the background, whose coefficients
  // other than the first are 1e-8 of it) have no peak.
  void Peaks_EachAxisOnce()
  {
    auto mrtrixImage = LoadShImage();
    auto fslImage = ToFsl(mrtrixImage);
    const mitk::PeakParameters parameters;
    const double closest = std::cos(parameters.angularThreshold * itk::Math::pi / 180);
    const itk::Index<3> antipodalPair = {{6, 4, 1}};

    for (mitk::ShImage *shImage : {mrtrixImage.GetPointer(), fslImage.GetPointer()})
    {
      const auto result = mitk::ExtractPeaks(shImage, parameters);
      auto peaks = result.peaks->GetItkPeakImage();
      itk::Image<unsigned char, 3>::Pointer numberOfPeaks;
      mitk::CastToItkImage(result.numberOfPeaks, numberOfPeaks);
      auto coefficients = mitk::ImageToItkImage<itk::Vector<float, 15>, 3>(mrtrixImage.GetPointer());

      for (itk::ImageRegionConstIteratorWithIndex<ShItkImageType> it(coefficients, coefficients->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
      {
        const auto voxel = it.GetIndex();
        const unsigned int count = numberOfPeaks->GetPixel(voxel);
        std::vector<Direction> axes;
        for (unsigned int p = 0; p < count; ++p)
        {
          Direction peak;
          for (unsigned int c = 0; c < 3; ++c)
            peak[c] = peaks->GetPixel({{voxel[0], voxel[1], voxel[2], 3 * p + c}});
          for (const auto &other : axes)
            CPPUNIT_ASSERT(std::abs(dot_product(peak, other)) <= closest * peak.two_norm() * other.two_norm());
          axes.push_back(peak);
        }

        double anisotropic = 0;
        for (unsigned int j = 1; j < 15; ++j)
          anisotropic = std::max(anisotropic, std::abs(static_cast<double>(it.Get()[j])));
        if (anisotropic < 1e-6 * it.Get()[0])
          CPPUNIT_ASSERT_EQUAL(0u, count);
        if (voxel == antipodalPair)
          CPPUNIT_ASSERT_EQUAL(1u, count);
      }
    }
  }

  // The peak image is 4D; mitk::GrabItkImageMemory would take over only its first volume and leave
  // the others uninitialized, so ExtractPeaks takes the whole channel. Every volume, the last
  // included, equals the filter's.
  void Peaks_SameAsFilter_AllVolumes()
  {
    auto shImage = LoadShImage();
    mitk::PeakParameters parameters;
    parameters.maxPeaks = 2;
    const auto result = mitk::ExtractPeaks(shImage, parameters);

    using FilterType = itk::OdfMaximaExtractionFilter<float, 4, ODF_SAMPLING_SIZE>;
    auto filter = FilterType::New();
    filter->SetInput(mitk::ImageToItkImage<ShItkImageType::PixelType, 3>(shImage.GetPointer()));
    filter->SetToolkit(FilterType::MRTRIX);
    filter->SetAngularThreshold(std::cos(parameters.angularThreshold * itk::Math::pi / 180));
    filter->SetMaxNumPeaks(parameters.maxPeaks);
    filter->SetRelativePeakThreshold(parameters.relativeThreshold);
    filter->SetAbsolutePeakThreshold(parameters.absoluteThreshold);
    filter->SetNormalizationMethod(FilterType::MAX_VEC_NORM);
    filter->Update();

    auto peaks = result.peaks->GetItkPeakImage();
    CPPUNIT_ASSERT(SameVoxels(filter->GetPeakImage().GetPointer(), peaks.GetPointer()));
    const auto size = peaks->GetLargestPossibleRegion().GetSize();
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(6), size[3]);
    bool lastVolumeSet = false;
    for (itk::IndexValueType z = 0; z < static_cast<itk::IndexValueType>(size[2]); ++z)
      for (itk::IndexValueType y = 0; y < static_cast<itk::IndexValueType>(size[1]); ++y)
        for (itk::IndexValueType x = 0; x < static_cast<itk::IndexValueType>(size[0]); ++x)
          lastVolumeSet = lastVolumeSet || peaks->GetPixel({{x, y, z, 5}}) != 0;
    CPPUNIT_ASSERT(lastVolumeSet);
  }

  static mitk::Image::Pointer MaskLargerThan(const mitk::Image *image)
  {
    auto mask = itk::Image<unsigned char, 3>::New();
    itk::ImageRegion<3> region;
    for (unsigned int d = 0; d < 3; ++d)
      region.SetSize(d, image->GetDimension(d) + (d == 0 ? 1 : 0));
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(1);
    return mitk::GrabItkImageMemory(mask.GetPointer());
  }

  void Peaks_MaskOfOtherSize_Throws()
  {
    auto shImage = LoadShImage();
    CPPUNIT_ASSERT_THROW(mitk::ExtractPeaks(shImage, mitk::PeakParameters(), MaskLargerThan(shImage)), itk::ExceptionObject);
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::ExtractTensorPeaks(tensors, 0.2, true, MaskLargerThan(tensors)), itk::ExceptionObject);
  }

  void TensorPeaks_SameAsFilter()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    for (bool normalize : {false, true})
    {
      auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
      filter->SetInput(ItkTensors(tensors));
      filter->SetUsePolarCoordinates(false);
      filter->SetNormalizeVectors(normalize);
      filter->SetFaThreshold(0.1);
      filter->Update();

      const auto result = mitk::ExtractTensorPeaks(tensors, 0.1, normalize);
      CPPUNIT_ASSERT(SameVoxels(filter->GetPeakImage().GetPointer(), result.peaks->GetItkPeakImage().GetPointer()));
      itk::Image<unsigned char, 3>::Pointer numberOfPeaks;
      mitk::CastToItkImage(result.numberOfPeaks, numberOfPeaks);
      CPPUNIT_ASSERT(SameVoxels(filter->GetOutput(), numberOfPeaks.GetPointer()));
    }
  }

  static mitk::Image::Pointer LoadImage(const std::string &path)
  {
    return mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(path));
  }

  // mitk-diffusion-test-data holds the indices of tensorImage.dti; the function reproduces FA, AD and RD.
  // tensorImage_CA.nrrd and tensorImage_MD.nrrd are older than this repository's history and differ
  // from today's filter (CA: 106 voxels, up to 0.031; MD: 12096 voxels, up to 1.5e-3 at values of
  // about 1e-3, so not the mean eigenvalue), and tensorImage_RA.nrrd holds ITK's RA formula, which
  // the filter replaced by Basser and Pierpaoli's (mitkDiffusionModellingAlgorithmsTest checks it
  // against the eigenvalues), so CA, MD and RA are compared with the filter.
  void TensorIndices_MatchReferences()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti"));
    const std::vector<std::pair<mitk::TensorIndex, std::string>> cases = {
      {mitk::TensorIndex::FA, "FA"}, {mitk::TensorIndex::AD, "AD"}, {mitk::TensorIndex::RD, "RD"}};
    for (const auto &[index, name] : cases)
    {
      auto reference = LoadImage("tensorImage_" + name + ".nrrd");
      auto result = mitk::TensorIndexImage(tensors, index, 1.0);
      CPPUNIT_ASSERT_MESSAGE(name, mitk::Equal(*result, *reference, 1e-5, true));
    }

    using FilterType = itk::TensorDerivedMeasurementsFilter<float>;
    for (const auto &[index, measure] : {std::pair{mitk::TensorIndex::CA, FilterType::CA},
                                         std::pair{mitk::TensorIndex::MD, FilterType::MD},
                                         std::pair{mitk::TensorIndex::RA, FilterType::RA}})
    {
      auto filter = FilterType::New();
      filter->SetInput(ItkTensors(tensors));
      filter->SetMeasure(measure);
      filter->Update();
      itk::Image<float, 3>::Pointer result;
      mitk::CastToItkImage(mitk::TensorIndexImage(tensors, index, 1.0), result);
      CPPUNIT_ASSERT(SameVoxels(filter->GetOutput(), result.GetPointer()));
    }
  }

  void Gfa_MatchesReference()
  {
    auto odfs = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi"));
    auto reference = LoadImage("qBallImage_GFA.nrrd");
    CPPUNIT_ASSERT(mitk::Equal(*mitk::GfaImage(odfs, 1.0), *reference, 1e-5, true));
  }

  void Adc_SameAsFilter()
  {
    DwiType::Pointer itkDwi;
    mitk::CastToItkImage(m_Dwi, itkDwi);
    for (bool fit : {false, true})
    {
      auto filter = itk::AdcImageFilter<short, double>::New();
      filter->SetInput(itkDwi);
      filter->SetGradientDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(m_Dwi));
      filter->SetB_value(mitk::DiffusionPropertyHelper::GetReferenceBValue(m_Dwi));
      filter->SetFitSignal(fit);
      filter->Update();

      itk::Image<double, 3>::Pointer adc;
      mitk::CastToItkImage(mitk::AdcImage(m_Dwi, fit, 1.0), adc);
      CPPUNIT_ASSERT(SameVoxels(filter->GetOutput(), adc.GetPointer()));
    }
  }

  template <class TPixel>
  static void CheckScaled(const mitk::Image *unscaled, const mitk::Image *scaled, double factor)
  {
    typename itk::Image<TPixel, 3>::Pointer a;
    typename itk::Image<TPixel, 3>::Pointer b;
    mitk::CastToItkImage(unscaled, a);
    mitk::CastToItkImage(scaled, b);
    itk::ImageRegionConstIterator<itk::Image<TPixel, 3>> itA(a, a->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<itk::Image<TPixel, 3>> itB(b, b->GetLargestPossibleRegion());
    unsigned int nonZero = 0;
    for (; !itA.IsAtEnd(); ++itA, ++itB)
    {
      CPPUNIT_ASSERT_EQUAL(static_cast<TPixel>(itA.Get() * factor), itB.Get());
      nonZero += itA.Get() != 0 ? 1 : 0;
    }
    CPPUNIT_ASSERT(nonZero > 0);
  }

  // The view's scale factor multiplies the index (e.g. diffusivities in 10^-3 mm^2/s).
  void Indices_Scale()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CheckScaled<float>(mitk::TensorIndexImage(tensors, mitk::TensorIndex::MD, 1), mitk::TensorIndexImage(tensors, mitk::TensorIndex::MD, 1000), 1000);
    auto odfs = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi"));
    CheckScaled<float>(mitk::GfaImage(odfs, 1), mitk::GfaImage(odfs, 2.5), 2.5);
    CheckScaled<double>(mitk::AdcImage(m_Dwi, false, 1), mitk::AdcImage(m_Dwi, false, 1000), 1000);
  }

  void Indices_WrongInput_Throws()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT_THROW(mitk::TensorIndexImage(nullptr, mitk::TensorIndex::FA, 1), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::GfaImage(tensors, 1), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::GfaImage(nullptr, 1), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::AdcImage(tensors, true, 1), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstruction)
