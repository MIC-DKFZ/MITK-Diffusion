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
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <itkAdcImageFilter.h>
#include <itkDiffusionOdfGeneralizedFaImageFilter.h>
#include <itkDiffusionTensorPrincipalDirectionImageFilter.h>
#include <itkFlipPeaksFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkOdfMaximaExtractionFilter.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkOdfToRgbImageFilter.h>
#include <itkShCoefficientImageImporter.h>
#include <itkShToRgbImageFilter.h>
#include <itkTensorDerivedMeasurementsFilter.h>
#include <itkTensorImageToDiffusionImageFilter.h>
#include <itkTensorImageToOdfImageFilter.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>
#include <itkTensorToL2NormImageFilter.h>
#include <itkTensorToRgbImageFilter.h>

#include <algorithm>
#include <cstdint>
#include <cstring>

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
  MITK_TEST(PrincipalDirection_MaskOfOtherSize_Throws);
  MITK_TEST(PrincipalDirection_SecondInput_SameAsFreshFilter);
  MITK_TEST(MaximaExtraction_MaskOfOtherSize_Throws);
  MITK_TEST(ShToRgb_MatchesOdfOfBothConventions);
  MITK_TEST(TensorToDwi_RoundTrip);
  MITK_TEST(TensorToDwi_SecondInput_SameAsFreshFilter);
  MITK_TEST(TensorToDwi_Mask_AppliesToBaseline);
  MITK_TEST(OdfToRgb_ZeroOdf_Black);
  MITK_TEST(Sh_InvalidCoefficientCounts_Throw);
  MITK_TEST(PrincipalDirection_FaThreshold_NoDirection);
  MITK_TEST(TensorMeasures_MatchReferences);
  MITK_TEST(TensorMeasures_RelativeAnisotropy_UnitIndependent);
  MITK_TEST(Gfa_MatchesReference);
  MITK_TEST(ShToOdf_MatchesReference);
  MITK_TEST(ShBasis_ConventionsAndEvaluationAgree);
  MITK_TEST(Convert_WrongImageType_Throws);
  MITK_TEST(RgbFilters_RecordedOutput);
  MITK_TEST(TensorToOdf_RecordedOutput);
  MITK_TEST(PrincipalDirection_RecordedOutput);
  MITK_TEST(FlipPeaks_NegatesChosenAxes);
  MITK_TEST(MaximaExtraction_RecordedOutput);
  MITK_TEST(Filters_OutputTakenBeforeUpdate_GetsResult);
  MITK_TEST(EigenvalueCorrection_RecordedOutput);
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

  // Corrected tensors of diffusionImage.hdwi as computed before the filter was restructured: a hash
  // of the bits of all tensor components, so any numerical change shows.
  void EigenvalueCorrection_RecordedOutput()
  {
    auto tensors = CorrectEigenvalues(LoadDwi("DiffusionImaging/diffusionImage.hdwi"));
    std::uint64_t hash = 14695981039346656037ull;
    itk::ImageRegionConstIterator<TensorImageType> it(tensors, tensors->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      for (unsigned int i = 0; i < 6; ++i)
      {
        float value = it.Get()[i];
        std::uint32_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        hash = (hash ^ bits) * 1099511628211ull;
      }
    CPPUNIT_ASSERT_EQUAL(std::uint64_t(5156811477493567636ull), hash);
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

  using MaskImageType = itk::Image<unsigned char, 3>;

  template <class TImage>
  static MaskImageType::Pointer MaskLargerThan(const TImage *image)
  {
    auto mask = MaskImageType::New();
    auto region = image->GetLargestPossibleRegion();
    region.SetSize(0, region.GetSize(0) + 1);
    mask->SetRegions(region);
    mask->SetSpacing(image->GetSpacing());
    mask->SetOrigin(image->GetOrigin());
    mask->SetDirection(image->GetDirection());
    mask->Allocate();
    mask->FillBuffer(1);
    return mask;
  }

  mitk::TensorImage::ItkTensorImageType::Pointer LoadItkTensors(const std::string &path)
  {
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath(path));
    return mitk::convert::GetItkTensorFromTensorImage(tensorImage.GetPointer());
  }

  // Masks are read with the input's voxel indices; another size meant out-of-bounds access.
  void PrincipalDirection_MaskOfOtherSize_Throws()
  {
    auto tensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    filter->SetInput(tensors);
    filter->SetMaskImage(MaskLargerThan(tensors.GetPointer()));
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  // The all-ones mask generated without a user mask was kept and reused for the next input.
  void PrincipalDirection_SecondInput_SameAsFreshFilter()
  {
    auto smallTensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    auto tensors = LoadItkTensors("DiffusionImaging/FWE/ReferenceTensorImage.dti");

    auto reused = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    reused->SetInput(smallTensors);
    reused->Update();
    reused->SetInput(tensors);
    reused->Update();

    auto fresh = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    fresh->SetInput(tensors);
    fresh->Update();

    using PeakImageType = itk::Image<float, 4>;
    itk::ImageRegionConstIterator<PeakImageType> itReused(reused->GetPeakImage(), reused->GetPeakImage()->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<PeakImageType> itFresh(fresh->GetPeakImage(), fresh->GetPeakImage()->GetLargestPossibleRegion());
    CPPUNIT_ASSERT(reused->GetPeakImage()->GetLargestPossibleRegion() == fresh->GetPeakImage()->GetLargestPossibleRegion());
    for (; !itFresh.IsAtEnd(); ++itFresh, ++itReused)
      CPPUNIT_ASSERT_EQUAL(itFresh.Get(), itReused.Get());
  }

  void MaximaExtraction_MaskOfOtherSize_Throws()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    auto coefficients = mitk::convert::GetItkShFromShImage<15>(shImage.GetPointer());

    auto filter = itk::OdfMaximaExtractionFilter<float, 4, 20242>::New();
    filter->SetInput(coefficients);
    filter->SetMaskImage(MaskLargerThan(coefficients.GetPointer()));
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  using RgbaImageType = itk::Image<itk::RGBAPixel<unsigned char>, 3>;

  static RgbaImageType::Pointer ShToRgb(itk::Image<itk::Vector<float, 15>, 3> *coefficients, mitk::ShImage::SH_CONVENTION convention, unsigned int workUnits)
  {
    auto filter = itk::ShToRgbImageFilter<itk::Image<itk::Vector<float, 15>, 3>, 4>::New();
    filter->SetInput(coefficients);
    filter->SetShConvention(convention);
    filter->SetNumberOfWorkUnits(workUnits);
    filter->Update();
    return filter->GetOutput();
  }

  // Each work unit used to reallocate the shared output and recompute the SH basis; the colour of
  // each voxel is its ODF's principal direction times GFA, in the image's SH convention.
  void ShToRgb_MatchesOdfOfBothConventions()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    auto coefficients = mitk::convert::GetItkShFromShImage<15>(shImage.GetPointer());

    for (auto convention : {mitk::ShImage::MRTRIX, mitk::ShImage::FSL})
    {
      auto rgbSingle = ShToRgb(coefficients, convention, 1);
      auto rgbMulti = ShToRgb(coefficients, convention, 8);

      shImage->SetShConvention(convention);
      auto odfs = mitk::convert::GetItkOdfFromShImage(shImage.GetPointer());

      CPPUNIT_ASSERT(rgbMulti->GetLargestPossibleRegion() == coefficients->GetLargestPossibleRegion());
      itk::ImageRegionConstIterator<RgbaImageType> itSingle(rgbSingle, rgbSingle->GetLargestPossibleRegion());
      itk::ImageRegionConstIterator<RgbaImageType> itMulti(rgbMulti, rgbMulti->GetLargestPossibleRegion());
      itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> itOdf(odfs, odfs->GetLargestPossibleRegion());
      unsigned int coloured = 0;
      for (; !itSingle.IsAtEnd(); ++itSingle, ++itMulti, ++itOdf)
      {
        CPPUNIT_ASSERT(itSingle.Get() == itMulti.Get());

        itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> odf;
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = itOdf.Get()[i];
        const int pd = odf.GetPrincipalDiffusionDirectionIndex();
        const auto dir = pd == -1 ? vnl_vector_fixed<double, 3>(0.0) : odf.GetDirection(pd);
        const double gfa = std::clamp<double>(odf.GetGeneralizedFractionalAnisotropy(), 0, 1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(std::clamp(std::fabs(dir[0]) * gfa, 0.0, 1.0) * 255, itSingle.Get().GetRed(), 1.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(std::clamp(std::fabs(dir[2]) * gfa, 0.0, 1.0) * 255, itSingle.Get().GetBlue(), 1.0);
        if (itSingle.Get().GetAlpha() > 0)
          ++coloured;
      }
      CPPUNIT_ASSERT(coloured > 0);
    }
  }

  using DwiToTensorFilterType = itk::TensorImageToDiffusionImageFilter<float, short>;

  static DwiToTensorFilterType::GradientListType::Pointer TensorToDwiGradients()
  {
    auto gradients = DwiToTensorFilterType::GradientListType::New();
    const double g[13][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.7071068, 0.7071068, 0}, {0.7071068, 0, 0.7071068},
                             {0, 0.7071068, 0.7071068}, {0.7071068, -0.7071068, 0}, {0.7071068, 0, -0.7071068}, {0, 0.7071068, -0.7071068},
                             {0.5773503, 0.5773503, 0.5773503}, {-0.5773503, 0.5773503, 0.5773503}, {0.5773503, -0.5773503, 0.5773503}};
    for (const auto &direction : g)
      gradients->push_back(DwiToTensorFilterType::GradientType(direction[0], direction[1], direction[2]));
    return gradients;
  }

  static itk::VectorImage<short, 3>::Pointer TensorToDwi(TensorImageType *tensors, itk::Image<short, 3> *mask = nullptr)
  {
    auto filter = DwiToTensorFilterType::New();
    filter->SetInput(tensors);
    filter->SetGradientList(TensorToDwiGradients().GetPointer());
    filter->SetBValue(1000);
    if (mask != nullptr)
      filter->SetMaskImage(mask);
    filter->Update();
    return filter->GetOutput();
  }

  // Simulated signals fitted with ITK's tensor reconstruction give the input tensors back.
  void TensorToDwi_RoundTrip()
  {
    auto tensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    auto dwi = TensorToDwi(tensors);
    CPPUNIT_ASSERT_EQUAL(13u, dwi->GetNumberOfComponentsPerPixel());

    auto reconstruction = itk::DiffusionTensor3DReconstructionImageFilter<short, short, float>::New();
    reconstruction->SetGradientImage(TensorToDwiGradients().GetPointer(), dwi);
    reconstruction->SetBValue(1000);
    reconstruction->Update();

    unsigned int compared = 0;
    itk::ImageRegionConstIterator<TensorImageType> itIn(tensors, tensors->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> itOut(reconstruction->GetOutput(), tensors->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<itk::VectorImage<short, 3>> itDwi(dwi, dwi->GetLargestPossibleRegion());
    for (; !itIn.IsAtEnd(); ++itIn, ++itOut, ++itDwi)
    {
      // Weak baselines lose precision in the short signal.
      if (itDwi.Get()[0] < 2000)
        continue;
      for (unsigned int i = 0; i < 6; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(itIn.Get()[i], itOut.Get()[i], 1e-5);
      ++compared;
    }
    CPPUNIT_ASSERT(compared > 10);
  }

  // The baseline generated for the first input was kept and used for the next one.
  void TensorToDwi_SecondInput_SameAsFreshFilter()
  {
    auto smallTensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    auto tensors = LoadItkTensors("DiffusionImaging/FWE/ReferenceTensorImage.dti");

    auto reused = DwiToTensorFilterType::New();
    reused->SetGradientList(TensorToDwiGradients().GetPointer());
    reused->SetBValue(1000);
    reused->SetInput(smallTensors);
    reused->Update();
    reused->SetInput(tensors);
    reused->UpdateLargestPossibleRegion();

    auto fresh = TensorToDwi(tensors);
    CPPUNIT_ASSERT(reused->GetOutput()->GetLargestPossibleRegion() == fresh->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<itk::VectorImage<short, 3>> itReused(reused->GetOutput(), fresh->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<itk::VectorImage<short, 3>> itFresh(fresh, fresh->GetLargestPossibleRegion());
    for (; !itFresh.IsAtEnd(); ++itFresh, ++itReused)
      CPPUNIT_ASSERT(itFresh.Get() == itReused.Get());
  }

  // The mask zeroed only the diffusion-weighted volumes.
  void TensorToDwi_Mask_AppliesToBaseline()
  {
    auto tensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    auto mask = itk::Image<short, 3>::New();
    mask->SetRegions(tensors->GetLargestPossibleRegion());
    mask->Allocate();
    mask->FillBuffer(0);

    auto dwi = TensorToDwi(tensors, mask);
    itk::ImageRegionConstIterator<itk::VectorImage<short, 3>> it(dwi, dwi->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      for (unsigned int i = 0; i < dwi->GetNumberOfComponentsPerPixel(); ++i)
        CPPUNIT_ASSERT_EQUAL(short(0), it.Get()[i]);
    }
  }

  // The GFA of an all-zero ODF is NaN, which was converted to unsigned char unclamped.
  void OdfToRgb_ZeroOdf_Black()
  {
    using OdfImageType = itk::Image<itk::Vector<float, ODF_SAMPLING_SIZE>, 3>;
    auto odfs = OdfImageType::New();
    OdfImageType::RegionType region;
    region.SetSize({{2, 1, 1}});
    odfs->SetRegions(region);
    odfs->Allocate();
    itk::Vector<float, ODF_SAMPLING_SIZE> odf;
    odf.Fill(0);
    odfs->SetPixel({{0, 0, 0}}, odf);
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = std::fabs(itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>::GetDirection(i)[0]);
    odfs->SetPixel({{1, 0, 0}}, odf);

    auto filter = itk::OdfToRgbImageFilter<OdfImageType>::New();
    filter->SetInput(odfs);
    filter->Update();

    const auto zero = filter->GetOutput()->GetPixel({{0, 0, 0}});
    CPPUNIT_ASSERT_EQUAL(0, int(zero.GetRed()) + int(zero.GetGreen()) + int(zero.GetBlue()) + int(zero.GetAlpha()));
    const auto xOdf = filter->GetOutput()->GetPixel({{1, 0, 0}});
    CPPUNIT_ASSERT(xOdf.GetRed() > xOdf.GetGreen() && xOdf.GetRed() > xOdf.GetBlue());
  }

  // ShOrder returned garbage for counts that belong to no order, GetValue read past too few
  // coefficients, and the importer copied any number of volumes into a fixed-size vector.
  void Sh_InvalidCoefficientCounts_Throw()
  {
    const std::pair<int, unsigned int> valid[] = {{1, 0}, {6, 2}, {15, 4}, {28, 6}, {45, 8}, {66, 10}, {91, 12}};
    for (const auto &[count, order] : valid)
      CPPUNIT_ASSERT_EQUAL(order, mitk::sh::ShOrder(count));
    for (int count : {0, 5, 7, 90})
      CPPUNIT_ASSERT_THROW(mitk::sh::ShOrder(count), mitk::Exception);

    vnl_vector<float> coefficients(14, 0.1f);
    vnl_vector_fixed<double, 3> dir(0, 0, 1);
    CPPUNIT_ASSERT_THROW(mitk::sh::GetValue(coefficients, 4, dir, true), mitk::Exception);
    coefficients.set_size(15);
    coefficients.fill(0.1f);
    CPPUNIT_ASSERT_NO_THROW(mitk::sh::GetValue(coefficients, 4, dir, true));

    using ImporterType = itk::ShCoefficientImageImporter<float, 4>;
    auto volumes = ImporterType::InputImageType::New();
    ImporterType::InputImageType::RegionType region;
    region.SetSize({{1, 1, 1, 28}});
    volumes->SetRegions(region);
    volumes->Allocate();
    volumes->FillBuffer(0);
    auto importer = ImporterType::New();
    importer->SetInputImage(volumes);
    CPPUNIT_ASSERT_THROW(importer->GenerateData(), itk::ExceptionObject);
  }

  // Voxels whose peak was removed by the FA threshold were counted as having one direction.
  void PrincipalDirection_FaThreshold_NoDirection()
  {
    auto tensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");
    for (float threshold : {0.0f, 1.1f})
    {
      auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
      filter->SetInput(tensors);
      filter->SetFaThreshold(threshold);
      filter->Update();

      unsigned int withDirection = 0;
      itk::ImageRegionConstIterator<MaskImageType> it(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
        withDirection += it.Get();
      if (threshold > 1)
        CPPUNIT_ASSERT_EQUAL(0u, withDirection);
      else
        CPPUNIT_ASSERT(withDirection > 0);
    }
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

  using ScalarImageType = itk::Image<float, 3>;

  ScalarImageType::Pointer LoadScalarImage(const std::string &path)
  {
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath(path));
    ScalarImageType::Pointer itkImage;
    mitk::CastToItkImage(image, itkImage);
    return itkImage;
  }

  static ScalarImageType::Pointer TensorMeasure(mitk::TensorImage::ItkTensorImageType *tensors, itk::TensorDerivedMeasurementsFilter<float>::Measure measure)
  {
    auto filter = itk::TensorDerivedMeasurementsFilter<float>::New();
    filter->SetInput(tensors);
    filter->SetMeasure(measure);
    filter->Update();
    return filter->GetOutput();
  }

  /** Largest absolute difference between two scalar images of the same size. */
  static double MaxDifference(const ScalarImageType *expected, const ScalarImageType *actual)
  {
    CPPUNIT_ASSERT(expected->GetLargestPossibleRegion().GetSize() == actual->GetLargestPossibleRegion().GetSize());
    itk::ImageRegionConstIterator<ScalarImageType> itExpected(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<ScalarImageType> itActual(actual, actual->GetLargestPossibleRegion());
    double maxDifference = 0;
    for (; !itExpected.IsAtEnd(); ++itExpected, ++itActual)
      maxDifference = std::max(maxDifference, std::fabs(static_cast<double>(itExpected.Get()) - itActual.Get()));
    return maxDifference;
  }

  // FA, AD and RD reproduce the MITK-Data references. The MD reference was made by an older
  // version and is not used; MD is checked against the trace, L2 and L3 against the eigenvalues.
  // The RA reference holds ITK's unit-dependent formula; RA is checked against the eigenvalues.
  void TensorMeasures_MatchReferences()
  {
    using FilterType = itk::TensorDerivedMeasurementsFilter<float>;
    auto tensors = LoadItkTensors("DiffusionImaging/tensorImage.dti");

    const std::pair<FilterType::Measure, std::string> references[] = {
      {FilterType::FA, "FA"}, {FilterType::AD, "AD"}, {FilterType::RD, "RD"}};
    for (const auto &[measure, name] : references)
    {
      auto reference = LoadScalarImage("DiffusionImaging/tensorImage_" + name + ".nrrd");
      CPPUNIT_ASSERT_MESSAGE(name, MaxDifference(reference, TensorMeasure(tensors, measure)) < 1e-6);
    }

    auto ca = TensorMeasure(tensors, FilterType::CA);
    auto md = TensorMeasure(tensors, FilterType::MD);
    auto l2 = TensorMeasure(tensors, FilterType::L2);
    auto l3 = TensorMeasure(tensors, FilterType::L3);
    auto ra = TensorMeasure(tensors, FilterType::RA);
    itk::ImageRegionConstIterator<mitk::TensorImage::ItkTensorImageType> it(tensors, tensors->GetLargestPossibleRegion());
    unsigned int anisotropic = 0;
    for (; !it.IsAtEnd(); ++it)
    {
      const auto index = it.GetIndex();
      FilterType::TensorType::EigenValuesArrayType ev;
      it.Get().ComputeEigenValues(ev);
      // RA = standard deviation / mean of the eigenvalues; the float eigenvalues limit the
      // accuracy relative to the largest one, also in the mean.
      const double mean = (static_cast<double>(ev[0]) + ev[1] + ev[2]) / 3.0;
      if (mean > 1e-6)
      {
        const double deviation = std::sqrt(((ev[0] - mean) * (ev[0] - mean) + (ev[1] - mean) * (ev[1] - mean) + (ev[2] - mean) * (ev[2] - mean)) / 3.0);
        const double largest = std::max({std::fabs(ev[0]), std::fabs(ev[1]), std::fabs(ev[2])});
        CPPUNIT_ASSERT_DOUBLES_EQUAL(deviation / mean, ra->GetPixel(index), 1e-6 * (1 + deviation / mean) * (1 + largest / mean));
      }
      // MD comes from the float eigenvalues, so its error scales with the largest eigenvalue.
      const double scale = std::max({std::fabs(ev[0]), std::fabs(ev[1]), std::fabs(ev[2])});
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it.Get().GetTrace() / 3.0, md->GetPixel(index), 1e-5 * scale + 1e-12);
      CPPUNIT_ASSERT_EQUAL(ev[1], l2->GetPixel(index));
      CPPUNIT_ASSERT_EQUAL(ev[0], l3->GetPixel(index));
      if (ev[2] > 0 && ev[0] > 0)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 - (ev[0] + ev[1]) / (2.0 * ev[2]), ca->GetPixel(index), 1e-6);
        ++anisotropic;
      }
    }
    CPPUNIT_ASSERT(anisotropic > 0);
  }

  static double RelativeAnisotropy(const itk::DiffusionTensor3D<float> &tensor)
  {
    auto image = TensorImageType::New();
    TensorImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(tensor);
    return TensorMeasure(image, itk::TensorDerivedMeasurementsFilter<float>::RA)->GetPixel({{0, 0, 0}});
  }

  static itk::DiffusionTensor3D<float> Tensor(double xx, double xy, double xz, double yy, double yz, double zz)
  {
    const float elements[6] = {float(xx), float(xy), float(xz), float(yy), float(yz), float(zz)};
    return itk::DiffusionTensor3D<float>(elements);
  }

  // RA was ITK's sqrt((3|D|^2 - tr^2) / (sqrt(3) tr)), which depends on the diffusivity unit
  // (0.0314 in mm^2/s, 0.992 in um^2/ms for the prolate tensor below). Basser & Pierpaoli 1996:
  // std / mean of the eigenvalues, 0 to sqrt(2).
  void TensorMeasures_RelativeAnisotropy_UnitIndependent()
  {
    // Eigenvalues (1.7, 0.3, 0.3), principal axis (1, 1, 0) / sqrt(2): RA = 0.860818.
    for (double scale : {1e-3, 1.0, 1e-9})
    {
      const auto prolate = Tensor(1.0 * scale, 0.7 * scale, 0, 1.0 * scale, 0, 0.3 * scale);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.860818, RelativeAnisotropy(prolate), 1e-5);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0), RelativeAnisotropy(Tensor(1e-3, 0, 0, 0, 0, 0)), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, RelativeAnisotropy(Tensor(1e-3, 0, 0, 1e-3, 0, 1e-3)), 1e-6);
    CPPUNIT_ASSERT_EQUAL(0.0, RelativeAnisotropy(Tensor(-1e-3, 0, 0, -1e-3, 0, 0.5e-3)));
    CPPUNIT_ASSERT_EQUAL(0.0, RelativeAnisotropy(Tensor(0, 0, 0, 0, 0, 0)));
  }

  // GFA of the ODFs in qBallImage.qbi as stored in MITK-Data.
  void Gfa_MatchesReference()
  {
    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/qBallImage.qbi"));
    auto filter = itk::DiffusionOdfGeneralizedFaImageFilter<float, float, ODF_SAMPLING_SIZE>::New();
    filter->SetInput(mitk::convert::GetItkOdfFromOdfImage(odfImage.GetPointer()));
    filter->SetComputationMethod(itk::DiffusionOdfGeneralizedFaImageFilter<float, float, ODF_SAMPLING_SIZE>::GFA_STANDARD);
    filter->Update();
    CPPUNIT_ASSERT(MaxDifference(LoadScalarImage("DiffusionImaging/qBallImage_GFA.nrrd"), filter->GetOutput()) < 1e-6);
  }

  // The MRtrix CSD coefficients (order 8) sampled on the ODF directions give the ODF file that
  // was made from them. The file was made from the coefficients read in MITK's legacy basis.
  void ShToOdf_MatchesReference()
  {
    mitk::IFileReader::Options options;
    options["SH convention"] = std::string("MITK legacy (MRtrix 0.2)");
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, options);
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image_sh.nii.gz"), &functor);
    CPPUNIT_ASSERT_EQUAL(8u, shImage->ShOrder());
    auto odfs = mitk::convert::GetItkOdfFromShImage(shImage.GetPointer());
    auto reference = mitk::convert::GetItkOdfFromOdfImage(mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image.qbi")).GetPointer());

    CPPUNIT_ASSERT(odfs->GetLargestPossibleRegion() == reference->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> it(odfs, odfs->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> itReference(reference, reference->GetLargestPossibleRegion());
    double maxDifference = 0;
    for (; !it.IsAtEnd(); ++it, ++itReference)
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
        maxDifference = std::max(maxDifference, std::fabs(static_cast<double>(it.Get()[i]) - itReference.Get()[i]));
    CPPUNIT_ASSERT(maxDifference < 2e-6);

    shImage->SetShConvention(mitk::ShImage::FSL);
    auto odfsFsl = mitk::convert::GetItkOdfFromShImage(shImage.GetPointer());
    CPPUNIT_ASSERT(odfsFsl->GetPixel({{30, 30, 1}}) != odfs->GetPixel({{30, 30, 1}}));
  }

  /** Unit vectors spread evenly over the sphere (Fibonacci lattice), as columns. */
  static vnl_matrix<double> SphereDirections(unsigned int count)
  {
    vnl_matrix<double> directions(3, count);
    const double golden = itk::Math::pi * (3.0 - std::sqrt(5.0));
    for (unsigned int i = 0; i < count; ++i)
    {
      const double z = 1.0 - (2.0 * i + 1.0) / count;
      const double r = std::sqrt(1.0 - z * z);
      directions(0, i) = r * std::cos(golden * i);
      directions(1, i) = r * std::sin(golden * i);
      directions(2, i) = z;
    }
    return directions;
  }

  // Pins both SH conventions: the FSL/Dipy basis is orthonormal; MITK's "MRtrix" basis (MRtrix
  // 0.2) is orthogonal with norm 1/2 for m != 0. GetValue, the basis matrix and the reconstruction
  // filters' basis (spherical input) give the same values.
  void ShBasis_ConventionsAndEvaluationAgree()
  {
    const unsigned int order = 8;
    const unsigned int numCoeffs = 45;
    const auto directions = SphereDirections(20000);
    for (bool mrtrix : {true, false})
    {
      const vnl_matrix<float> floatBasis = mitk::sh::CalcShBasisForDirections(order, directions, mrtrix);
      vnl_matrix<double> basis(floatBasis.rows(), floatBasis.cols());
      for (unsigned int r = 0; r < basis.rows(); ++r)
        for (unsigned int c = 0; c < basis.cols(); ++c)
          basis(r, c) = floatBasis(r, c);
      CPPUNIT_ASSERT_EQUAL(numCoeffs, basis.cols());
      const vnl_matrix<double> gram = basis.transpose() * basis * (4.0 * itk::Math::pi / directions.cols());
      for (unsigned int k = 0, j = 0; k <= order; k += 2)
      {
        for (int m = -static_cast<int>(k); m <= static_cast<int>(k); ++m, ++j)
        {
          const double expectedNorm = (mrtrix && m != 0) ? 0.5 : 1.0;
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedNorm, gram(j, j), 1e-3);
          for (unsigned int j2 = 0; j2 < numCoeffs; ++j2)
            if (j2 != j)
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, gram(j, j2), 1e-3);
        }
      }
    }

    vnl_vector<float> coefficients(numCoeffs);
    for (unsigned int i = 0; i < numCoeffs; ++i)
      coefficients[i] = std::sin(1.0f + i);
    const auto someDirections = SphereDirections(50);
    vnl_matrix<double> spherical(3, someDirections.cols());
    for (unsigned int i = 0; i < someDirections.cols(); ++i)
    {
      double sph[3];
      mitk::gradients::Cart2Sph(someDirections(0, i), someDirections(1, i), someDirections(2, i), sph);
      for (unsigned int r = 0; r < 3; ++r)
        spherical(r, i) = sph[r];
    }
    const vnl_matrix<double> reconstructionBasis = mitk::sh::ComputeSphericalHarmonicsBasis(spherical, order);
    for (bool mrtrix : {true, false})
    {
      const vnl_vector<float> values = mitk::sh::CalcShBasisForDirections(order, someDirections, mrtrix) * coefficients;
      for (unsigned int i = 0; i < someDirections.cols(); ++i)
      {
        const vnl_vector_fixed<double, 3> dir(someDirections(0, i), someDirections(1, i), someDirections(2, i));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(values[i], mitk::sh::GetValue(coefficients, order, dir, mrtrix), 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(values[i], mitk::sh::GetValue(coefficients, order, spherical(1, i), spherical(0, i), mrtrix), 1e-5);
      }
    }
    const vnl_matrix<float> mrtrixBasis = mitk::sh::CalcShBasisForDirections(order, someDirections, true);
    for (unsigned int i = 0; i < someDirections.cols(); ++i)
      for (unsigned int j = 0; j < numCoeffs; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(mrtrixBasis(i, j), reconstructionBasis(i, j), 1e-6);
  }

  void Convert_WrongImageType_Throws()
  {
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    auto odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi"));
    CPPUNIT_ASSERT_THROW(mitk::convert::GetItkOdfFromShImage(tensorImage.GetPointer()), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::convert::GetOdfFromShImage(odfImage.GetPointer()), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::convert::GetItkShFromShImage<15>(odfImage.GetPointer()), mitk::Exception);
  }

  /** FNV-1a hash of all RGBA bytes in buffer order. */
  static std::uint64_t Hash(const RgbaImageType *image)
  {
    std::uint64_t hash = 14695981039346656037ull;
    itk::ImageRegionConstIterator<RgbaImageType> it(image, image->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      for (unsigned int c = 0; c < 4; ++c)
        hash = (hash ^ it.Get()[c]) * 1099511628211ull;
    return hash;
  }

  // Colour maps of tensor, ODF and SH images as computed before the RGB code was shared.
  void RgbFilters_RecordedOutput()
  {
    auto tensorRgb = itk::TensorToRgbImageFilter<mitk::TensorImage::ItkTensorImageType>::New();
    tensorRgb->SetInput(LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    tensorRgb->Update();

    auto odfRgb = itk::OdfToRgbImageFilter<mitk::OdfImage::ItkOdfImageType>::New();
    odfRgb->SetInput(mitk::convert::GetItkOdfFromOdfImage(mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image.qbi")).GetPointer()));
    odfRgb->Update();

    // csd_image.qbi was made from the coefficients read in MITK's legacy basis.
    mitk::IFileReader::Options options;
    options["SH convention"] = std::string("MITK legacy (MRtrix 0.2)");
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, options);
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/StreamlineTractography/csd_image_sh.nii.gz"), &functor);
    auto shRgb = itk::ShToRgbImageFilter<itk::Image<itk::Vector<float, 45>, 3>, 8>::New();
    shRgb->SetInput(mitk::convert::GetItkShFromShImage<45>(shImage.GetPointer()));
    shRgb->SetShConvention(mitk::ShImage::MRTRIX);
    shRgb->Update();

    CPPUNIT_ASSERT_EQUAL(std::uint64_t(7043315695314250899ull), Hash(tensorRgb->GetOutput()));
    // csd_image.qbi was sampled from csd_image_sh.nii.gz, so both give the same colours.
    CPPUNIT_ASSERT_EQUAL(std::uint64_t(17028873552555350804ull), Hash(odfRgb->GetOutput()));
    CPPUNIT_ASSERT_EQUAL(std::uint64_t(17028873552555350804ull), Hash(shRgb->GetOutput()));
  }

  // ODFs of test_dti.dti, summarized by the sum and the sum of squares of all values.
  void TensorToOdf_RecordedOutput()
  {
    auto filter = itk::TensorImageToOdfImageFilter<float, float>::New();
    filter->SetInput(LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    filter->Update();
    double sum = 0;
    double sumOfSquares = 0;
    itk::ImageRegionConstIterator<mitk::OdfImage::ItkOdfImageType> it(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      {
        sum += it.Get()[i];
        sumOfSquares += it.Get()[i] * it.Get()[i];
      }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(362.99999671805926, sum, 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.4703429208501879, sumOfSquares, 1e-9);
  }

  // Principal directions of test_dti.dti, summarized per component.
  void PrincipalDirection_RecordedOutput()
  {
    auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    filter->SetInput(LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti"));
    filter->Update();

    double sums[3] = {0, 0, 0};
    using PeakImageType = itk::Image<float, 4>;
    itk::ImageRegionConstIteratorWithIndex<PeakImageType> it(filter->GetPeakImage(), filter->GetPeakImage()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      sums[it.GetIndex()[3]] += std::fabs(it.Get());
    unsigned int withDirection = 0;
    itk::ImageRegionConstIterator<MaskImageType> itNum(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
    for (; !itNum.IsAtEnd(); ++itNum)
      withDirection += itNum.Get();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(33.890375792980194, sums[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.773321121931076, sums[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(12.989389510825276, sums[2], 1e-6);
    CPPUNIT_ASSERT_EQUAL(47u, withDirection);
  }

  // Peaks of sh_image_test.nii.gz as extracted before the 3D/4D geometry code was shared,
  // summarized per component, with the number of peaks and the 4D geometry.
  void MaximaExtraction_RecordedOutput()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    auto coefficients = mitk::convert::GetItkShFromShImage<15>(shImage.GetPointer());

    using FilterType = itk::OdfMaximaExtractionFilter<float, 4, ODF_SAMPLING_SIZE>;
    auto filter = FilterType::New();
    filter->SetInput(coefficients);
    filter->SetMaxNumPeaks(3);
    filter->Update();

    auto peaks = filter->GetPeakImage();
    double sums[3] = {0, 0, 0};
    itk::ImageRegionConstIteratorWithIndex<FilterType::PeakImageType> it(peaks, peaks->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      sums[it.GetIndex()[3] % 3] += std::fabs(it.Get());
    unsigned int numPeaks = 0;
    itk::ImageRegionConstIterator<MaskImageType> itNum(filter->GetNumDirectionsImage(), filter->GetNumDirectionsImage()->GetLargestPossibleRegion());
    for (; !itNum.IsAtEnd(); ++itNum)
      numPeaks += itNum.Get();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(303.91917738690972, sums[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(505.74287556356285, sums[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(414.1981383588427, sums[2], 1e-6);
    CPPUNIT_ASSERT_EQUAL(809u, numPeaks);

    CPPUNIT_ASSERT(filter->GetNumDirectionsImage()->GetLargestPossibleRegion() == coefficients->GetLargestPossibleRegion());
    CPPUNIT_ASSERT(filter->GetNumDirectionsImage()->IsSameImageGeometryAs(coefficients));
    const itk::Size<4> peakSize = {{10, 10, 3, 9}};
    CPPUNIT_ASSERT(peaks->GetLargestPossibleRegion().GetSize() == peakSize);
    for (unsigned int r = 0; r < 4; ++r)
    {
      CPPUNIT_ASSERT_EQUAL(r < 3 ? coefficients->GetSpacing()[r] : 1.0, peaks->GetSpacing()[r]);
      CPPUNIT_ASSERT_EQUAL(r < 3 ? coefficients->GetOrigin()[r] : 0.0, peaks->GetOrigin()[r]);
      for (unsigned int c = 0; c < 4; ++c)
        CPPUNIT_ASSERT_EQUAL(r < 3 && c < 3 ? coefficients->GetDirection()[r][c] : (r == c ? 1.0 : 0.0), peaks->GetDirection()[r][c]);
    }
  }

  // Three filters replaced their primary output during the update, so an output taken before
  // Update(), as in a pipeline, stayed empty.
  void Filters_OutputTakenBeforeUpdate_GetsResult()
  {
    auto tensors = LoadItkTensors("DiffusionImaging/ImageReconstruction/test_dti.dti");

    auto toOdf = itk::TensorImageToOdfImageFilter<float, float>::New();
    toOdf->SetInput(tensors);
    auto odfs = toOdf->GetOutput();
    toOdf->Update();
    CPPUNIT_ASSERT(odfs == toOdf->GetOutput());
    CPPUNIT_ASSERT(odfs->GetLargestPossibleRegion() == tensors->GetLargestPossibleRegion());

    auto principal = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
    principal->SetInput(tensors);
    auto numDirections = principal->GetOutput();
    principal->Update();
    CPPUNIT_ASSERT(numDirections == principal->GetOutput());

    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    auto shImage = mitk::IOUtil::Load<mitk::ShImage>(GetTestDataFilePath("DiffusionImaging/sh_image_test.nii.gz"), &functor);
    auto maxima = itk::OdfMaximaExtractionFilter<float, 4, ODF_SAMPLING_SIZE>::New();
    maxima->SetInput(mitk::convert::GetItkShFromShImage<15>(shImage.GetPointer()));
    auto numPeaks = maxima->GetOutput();
    maxima->Update();
    CPPUNIT_ASSERT(numPeaks == maxima->GetNumDirectionsImage());
  }

  void FlipPeaks_NegatesChosenAxes()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    auto peakImage = mitk::IOUtil::Load<mitk::PeakImage>(GetTestDataFilePath("DiffusionImaging/peak_image_test.nii.gz"), &functor);
    auto peaks = peakImage->GetItkPeakImage();
    const bool flips[][3] = {{true, false, false}, {false, true, true}};
    for (const auto &flip : flips)
    {
      auto filter = itk::FlipPeaksFilter<float>::New();
      filter->SetInput(peaks);
      filter->SetFlipX(flip[0]);
      filter->SetFlipY(flip[1]);
      filter->SetFlipZ(flip[2]);
      filter->Update();
      CPPUNIT_ASSERT(filter->GetOutput()->GetLargestPossibleRegion() == peaks->GetLargestPossibleRegion());
      itk::ImageRegionConstIteratorWithIndex<mitk::PeakImage::ItkPeakImageType> it(peaks, peaks->GetLargestPossibleRegion());
      itk::ImageRegionConstIterator<mitk::PeakImage::ItkPeakImageType> itFlipped(filter->GetOutput(), peaks->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it, ++itFlipped)
        CPPUNIT_ASSERT_EQUAL(flip[it.GetIndex()[3] % 3] ? -it.Get() : it.Get(), itFlipped.Get());
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingAlgorithms)
