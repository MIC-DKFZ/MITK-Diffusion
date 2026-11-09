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
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>
#include <itkImageRegionConstIterator.h>

#include <cmath>
#include <vector>

/** Tests of the tensor reconstruction function used by the reconstruction view. */
class mitkDiffusionReconstructionTensorTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionTensorTestSuite);
  MITK_TEST(LinearLeastSquares_MultiShell_RecoversTensor);
  MITK_TEST(LinearLeastSquares_ShellBelowReference_RecoversTensor);
  MITK_TEST(LinearLeastSquares_ShortB0Vector_IsB0);
  MITK_TEST(LinearLeastSquares_SingleShell_SameAsItk);
  MITK_TEST(LinearLeastSquares_InvalidGradients_Throw);
  MITK_TEST(EigenvalueCorrection_LowBVolumes_AreDiffusionWeighted);
  MITK_TEST(EigenvalueCorrection_B0AtThreshold_IsFitted);
  CPPUNIT_TEST_SUITE_END();

  using GradientContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using DwiType = mitk::DiffusionPropertyHelper::ImageType;
  using TensorImageType = mitk::TensorImage::ItkTensorImageType;
  using Direction = vnl_vector_fixed<double, 3>;

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

  /** Two b0 volumes followed by the given directions scaled to b = bValue for each b-value. */
  static GradientContainerType::Pointer Gradients(const std::vector<double> &bValues, double referenceBValue)
  {
    auto gradients = GradientContainerType::New();
    gradients->push_back(Direction(0.0));
    gradients->push_back(Direction(0.0));
    for (double b : bValues)
      for (const auto &d : Directions(30))
        gradients->push_back(d * std::sqrt(b / referenceBValue));
    return gradients;
  }

  /** The known tensor of voxel x of the simulated DWI: a cylinder along a direction rotating with x. */
  static itk::DiffusionTensor3D<double> Tensor(unsigned int x)
  {
    const double angle = 0.4 * x;
    const Direction e(std::cos(angle), std::sin(angle), 0.3);
    const Direction u = e / e.two_norm();
    itk::DiffusionTensor3D<double> tensor;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = i; j < 3; ++j)
        tensor(i, j) = 1.4e-3 * u[i] * u[j] + (i == j ? 0.3e-3 : 0.0);
    return tensor;
  }

  /** A 4x1x1 DWI of the tensors of Tensor(), S0 = 10000; b = referenceBValue * |g|^2. */
  static DwiType::Pointer Simulate(const GradientContainerType *gradients, double referenceBValue)
  {
    auto dwi = DwiType::New();
    DwiType::RegionType region;
    region.SetSize({{4, 1, 1}});
    dwi->SetRegions(region);
    dwi->SetVectorLength(gradients->Size());
    dwi->Allocate();
    for (unsigned int x = 0; x < 4; ++x)
    {
      const auto tensor = Tensor(x);
      DwiType::PixelType pixel(gradients->Size());
      for (unsigned int i = 0; i < gradients->Size(); ++i)
      {
        const auto &g = gradients->ElementAt(i);
        double gDg = 0;
        for (unsigned int r = 0; r < 3; ++r)
          for (unsigned int c = 0; c < 3; ++c)
            gDg += g[r] * tensor(r, c) * g[c];
        pixel[i] = static_cast<short>(std::round(10000 * std::exp(-referenceBValue * gDg)));
      }
      dwi->SetPixel({{x, 0, 0}}, pixel);
    }
    return dwi;
  }

  static TensorImageType::Pointer CorrectEigenvalues(const GradientContainerType *gradients, const DwiType *dwi, double referenceBValue, double b0Threshold)
  {
    auto filter = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>::New();
    filter->SetBValue(referenceBValue);
    filter->SetGradientImage(gradients, dwi);
    filter->SetB0Threshold(b0Threshold);
    filter->Update();
    return filter->GetOutput();
  }

  static void CheckRecovered(const TensorImageType *tensors, double tolerance)
  {
    for (unsigned int x = 0; x < 4; ++x)
    {
      const auto expected = Tensor(x);
      const auto actual = tensors->GetPixel({{x, 0, 0}});
      for (unsigned int i = 0; i < 6; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], actual[i], tolerance);
    }
  }

public:
  // ITK's filter fitted the b = 1000 and 2000 volumes of a DWI with reference b = 3000 as if they
  // had b = 3000: the tensors were far too small.
  void LinearLeastSquares_MultiShell_RecoversTensor()
  {
    auto gradients = Gradients({1000, 2000, 3000}, 3000);
    auto dwi = Simulate(gradients, 3000);
    CheckRecovered(mitk::ReconstructTensorsLinearLeastSquares(dwi, gradients, 3000, 0), 2e-5);
  }

  // One shell at b = 1000 in a DWI whose reference b-value is 3000, e.g. after shell extraction.
  void LinearLeastSquares_ShellBelowReference_RecoversTensor()
  {
    auto gradients = Gradients({1000}, 3000);
    auto dwi = Simulate(gradients, 3000);
    CheckRecovered(mitk::ReconstructTensorsLinearLeastSquares(dwi, gradients, 3000, 0), 2e-5);
  }

  // A b0 volume stored with a short non-zero vector (b = 0.1 here) is a b0 volume, as in the b-value
  // map; ITK's filter took it for a volume at the reference b-value.
  void LinearLeastSquares_ShortB0Vector_IsB0()
  {
    auto gradients = Gradients({1000}, 1000);
    gradients->ElementAt(1) = Direction(0.01, 0.0, 0.0);
    auto dwi = Simulate(gradients, 1000);
    CheckRecovered(mitk::ReconstructTensorsLinearLeastSquares(dwi, gradients, 1000, 0), 2e-5);
  }

  // DWIs whose diffusion-weighted volumes all have the reference b-value keep ITK's results.
  void LinearLeastSquares_SingleShell_SameAsItk()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"), &functor);
    auto dwi = DwiType::New();
    mitk::CastToItkImage(image, dwi);
    const double bValue = mitk::DiffusionPropertyHelper::GetReferenceBValue(image);
    auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(image);

    auto copy = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      copy->push_back(gradients->ElementAt(i));
    auto itkFilter = itk::DiffusionTensor3DReconstructionImageFilter<short, short, float>::New();
    itkFilter->SetBValue(bValue);
    itkFilter->SetGradientImage(copy, dwi);
    itkFilter->SetThreshold(50);
    itkFilter->Update();

    auto tensors = mitk::ReconstructTensorsLinearLeastSquares(dwi, gradients, bValue, 50);
    itk::ImageRegionConstIterator<TensorImageType> itExpected(itkFilter->GetOutput(), dwi->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> itActual(tensors, dwi->GetLargestPossibleRegion());
    for (; !itExpected.IsAtEnd(); ++itExpected, ++itActual)
      CPPUNIT_ASSERT(itExpected.Get() == itActual.Get());
  }

  // ITK's filter divided by a zero b0 count (SIGFPE) for a DWI without b0.
  void LinearLeastSquares_InvalidGradients_Throw()
  {
    auto gradients = Gradients({1000, 2000}, 2000);
    auto dwi = Simulate(gradients, 2000);

    auto noB0 = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      noB0->push_back(i < 2 ? Direction(0.0, 0.0, 1.0) : gradients->ElementAt(i));
    CPPUNIT_ASSERT_THROW(mitk::ReconstructTensorsLinearLeastSquares(dwi, noB0, 2000, 0), mitk::Exception);

    auto tooFew = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      tooFew->push_back(i < 7 ? gradients->ElementAt(i) : Direction(0.0));
    CPPUNIT_ASSERT_THROW(mitk::ReconstructTensorsLinearLeastSquares(dwi, tooFew, 2000, 0), mitk::Exception);

    auto shorter = GradientContainerType::New();
    for (unsigned int i = 1; i < gradients->Size(); ++i)
      shorter->push_back(gradients->ElementAt(i));
    CPPUNIT_ASSERT_THROW(mitk::ReconstructTensorsLinearLeastSquares(dwi, shorter, 2000, 0), mitk::Exception);
  }

  // The eigenvalue-corrected filter took every volume with b < 100 for a b0 volume, so volumes at
  // b = 90 lowered the S0 it divides by and all diffusivities came out too high. A volume is a b0
  // volume when its b-value rounds to 0, as in the b-value map and the linear least-squares fit.
  void EigenvalueCorrection_LowBVolumes_AreDiffusionWeighted()
  {
    auto gradients = Gradients({90, 1000}, 1000);
    auto dwi = Simulate(gradients, 1000);
    CheckRecovered(CorrectEigenvalues(gradients, dwi, 1000, 0), 2e-5);
  }

  // The filter fitted only voxels whose b0 mean exceeded the threshold, while the view's preview
  // keeps those that reach it.
  void EigenvalueCorrection_B0AtThreshold_IsFitted()
  {
    auto gradients = Gradients({1000}, 1000);
    auto dwi = Simulate(gradients, 1000);
    CheckRecovered(CorrectEigenvalues(gradients, dwi, 1000, 10000), 2e-5);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionTensor)
