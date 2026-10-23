/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkException.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkPointShell.h>

#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <itkDiffusionMultiShellQballReconstructionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>

#include <cmath>
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

public:
  void setUp() override
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"), &functor);
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
      auto reference = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/" + file));
      CPPUNIT_ASSERT_MESSAGE(file, result.odf.IsNotNull());
      CPPUNIT_ASSERT_MESSAGE(file, mitk::Equal(*result.odf, *reference, tolerance, true));
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
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
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
    auto reference = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    CPPUNIT_ASSERT(mitk::Equal(*tensors, *reference, 1e-4, true));
  }

  static mitk::Image::Pointer LoadDwi(const std::string &path)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath(path), &functor);
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
    auto dwi = LoadDwi("DiffusionImaging/diffusionImage.hdwi");
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
    auto dwi = LoadDwi("DiffusionImaging/diffusionImage.hdwi");
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
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
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
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstruction)
