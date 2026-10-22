/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkTensorImage.h>

#include <itkAdcImageFilter.h>
#include <itkDiffusionTensorPrincipalDirectionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>
#include <itkTensorToL2NormImageFilter.h>

#ifdef _OPENMP
#include <omp.h>
#endif

class mitkDiffusionModellingAlgorithmsTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingAlgorithmsTestSuite);
  MITK_TEST(EigenvalueCorrection_ThreadCount_SameResult);
  MITK_TEST(PrincipalDirection_PeakImageGeometry);
  MITK_TEST(EigenvalueCorrection_InvalidGradients_Throws);
  MITK_TEST(TensorToL2Norm_FrobeniusNorm);
  MITK_TEST(Adc_NearZeroBValue_CountsAsBaseline);
  CPPUNIT_TEST_SUITE_END();

  using TensorImageType = itk::Image<itk::DiffusionTensor3D<float>, 3>;

  mitk::Image::Pointer LoadDwi(const std::string &path)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath(path), &functor);
  }

  TensorImageType::Pointer CorrectEigenvalues(mitk::Image::Pointer dwi)
  {
    auto itkDwi = itk::VectorImage<short, 3>::New();
    mitk::CastToItkImage(dwi, itkDwi);

    using FilterType = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>;
    auto filter = FilterType::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi);
    filter->SetB0Threshold(0);
    filter->Update();
    return filter->GetOutput();
  }

public:
  // The filter smooths the DWI in place, voxel after voxel, so its result must not depend on the
  // number of OpenMP threads. diffusionImage.hdwi has about 2000 voxels with negative eigenvalues,
  // so the correction loop runs several iterations.
  void EigenvalueCorrection_ThreadCount_SameResult()
  {
    auto dwi = LoadDwi("DiffusionImaging/diffusionImage.hdwi");

#ifdef _OPENMP
    const int maxThreads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    auto serial = CorrectEigenvalues(dwi);
#ifdef _OPENMP
    omp_set_num_threads(std::max(maxThreads, 4));
#endif
    auto parallel1 = CorrectEigenvalues(dwi);
    auto parallel2 = CorrectEigenvalues(dwi);
#ifdef _OPENMP
    omp_set_num_threads(maxThreads);
#endif

    unsigned int differing = 0;
    unsigned int nonZero = 0;
    itk::ImageRegionConstIterator<TensorImageType> it0(serial, serial->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> it1(parallel1, parallel1->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> it2(parallel2, parallel2->GetLargestPossibleRegion());
    for (; !it0.IsAtEnd(); ++it0, ++it1, ++it2)
    {
      for (unsigned int i = 0; i < 6; ++i)
      {
        if (it0.Get()[i] != it1.Get()[i] || it0.Get()[i] != it2.Get()[i])
          ++differing;
        if (it0.Get()[i] != 0)
          ++nonZero;
      }
    }
    CPPUNIT_ASSERT_MESSAGE("Corrected tensors are not all zero", nonZero > 0);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Tensor components that differ between 1 and N threads", 0u, differing);
  }

  /** Runs the eigenvalue correction on a constant 2x2x2 DWI with the given gradients (b=1000). */
  static void CorrectEigenvalues(const std::vector<vnl_vector_fixed<double, 3>> &directions)
  {
    using FilterType = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>;
    auto dwi = itk::VectorImage<short, 3>::New();
    itk::VectorImage<short, 3>::RegionType region;
    region.SetSize({{2, 2, 2}});
    dwi->SetRegions(region);
    dwi->SetVectorLength(directions.size());
    dwi->Allocate();
    itk::VariableLengthVector<short> pixel(directions.size());
    pixel.Fill(500);
    dwi->FillBuffer(pixel);

    auto gradients = FilterType::GradientDirectionContainerType::New();
    for (const auto &direction : directions)
      gradients->push_back(direction);

    auto filter = FilterType::New();
    filter->SetBValue(1000);
    filter->SetGradientImage(gradients.GetPointer(), dwi);
    filter->Update();
  }

  // Without a b0 volume every voxel was divided by zero and the output was silently all zero;
  // fewer than six diffusion-weighted volumes gave a rank-deficient fit.
  void EigenvalueCorrection_InvalidGradients_Throws()
  {
    std::vector<vnl_vector_fixed<double, 3>> directions;
    const double unit[7][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.7071, 0.7071, 0}, {0.7071, 0, 0.7071}, {0, 0.7071, 0.7071}, {0.5774, 0.5774, 0.5774}};
    for (const auto &u : unit)
      directions.emplace_back(u[0], u[1], u[2]);
    CPPUNIT_ASSERT_THROW(CorrectEigenvalues(directions), itk::ExceptionObject);

    directions.resize(5);
    directions.emplace_back(0, 0, 0);
    CPPUNIT_ASSERT_THROW(CorrectEigenvalues(directions), itk::ExceptionObject);
  }

  // The norm used for the synthesized baseline of TensorImageToDiffusionImageFilter mixed
  // products of different off-diagonal elements instead of their squares.
  void TensorToL2Norm_FrobeniusNorm()
  {
    auto image = TensorImageType::New();
    TensorImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    image->SetRegions(region);
    image->Allocate();
    const float elements[6] = {1, 2, 3, 4, 5, 6}; // xx, xy, xz, yy, yz, zz
    image->FillBuffer(itk::DiffusionTensor3D<float>(elements));

    auto filter = itk::TensorToL2NormImageFilter<TensorImageType, itk::Image<float, 3>>::New();
    filter->SetInput(image);
    filter->Update();

    const double expected = std::sqrt(1.0 + 16 + 36 + 2 * (4 + 9 + 25));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, filter->GetOutput()->GetPixel({{0, 0, 0}}), 1e-5);
  }

  static double Adc(bool fitSignal)
  {
    // b_ref = 1000. Volume 1 has |g| = 0.02, i.e. b = 0.4 s/mm^2: unweighted, but above the
    // |g| > 0.001 limit that one part of the filter used. Six volumes at b = 1000 with D = 0.001.
    const double D = 0.001;
    const double g[8][3] = {{0, 0, 0}, {0.02, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.7071068, 0.7071068, 0}, {0.7071068, 0, 0.7071068}, {0, 0.7071068, 0.7071068}};

    auto dwi = itk::VectorImage<short, 3>::New();
    itk::VectorImage<short, 3>::RegionType region;
    region.SetSize({{1, 1, 1}});
    dwi->SetRegions(region);
    dwi->SetVectorLength(8);
    dwi->Allocate();
    itk::VariableLengthVector<short> pixel(8);
    auto gradients = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < 8; ++i)
    {
      mitk::DiffusionPropertyHelper::GradientDirectionType direction;
      direction[0] = g[i][0];
      direction[1] = g[i][1];
      direction[2] = g[i][2];
      gradients->push_back(direction);
      const double b = 1000 * direction.squared_magnitude();
      pixel[i] = static_cast<short>(std::round(10000 * std::exp(-(b > 1 ? b : 0) * D)));
    }
    dwi->FillBuffer(pixel);

    auto filter = itk::AdcImageFilter<short, double>::New();
    filter->SetInput(dwi);
    filter->SetGradientDirections(gradients.GetPointer());
    filter->SetB_value(1000);
    filter->SetFitSignal(fitSignal);
    filter->Update();
    return filter->GetOutput()->GetPixel({{0, 0, 0}});
  }

  void Adc_NearZeroBValue_CountsAsBaseline()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.001, Adc(false), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.001, Adc(true), 1e-6);
  }

  // The 4D peak image takes the tensor image's geometry; its fourth origin component was left
  // uninitialized while a write went past the 3D origin.
  void PrincipalDirection_PeakImageGeometry()
  {
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/FWE/ReferenceTensorImage.dti"));
    auto itkTensorImage = mitk::convert::GetItkTensorFromTensorImage(tensorImage.GetPointer());

    auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    filter->SetInput(itkTensorImage);
    filter->Update();
    auto peaks = filter->GetPeakImage();

    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(itkTensorImage->GetOrigin()[i], peaks->GetOrigin()[i], 1e-4);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(itkTensorImage->GetSpacing()[i], peaks->GetSpacing()[i], 1e-9);
    }
    CPPUNIT_ASSERT_EQUAL(0.0, peaks->GetOrigin()[3]);
    CPPUNIT_ASSERT_EQUAL(1.0, peaks->GetSpacing()[3]);
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(3), peaks->GetLargestPossibleRegion().GetSize(3));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingAlgorithms)
