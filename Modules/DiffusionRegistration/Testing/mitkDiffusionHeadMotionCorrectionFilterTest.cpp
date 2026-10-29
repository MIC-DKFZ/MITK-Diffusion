/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkShrinkImageFilter.h>
#include <mitkDWIHeadMotionCorrectionFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

class mitkDiffusionHeadMotionCorrectionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionHeadMotionCorrectionFilterTestSuite);
  MITK_TEST(Correct_ReferenceNotFirst_KeepsVolumeOrder);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using FloatImageType = itk::Image<float, 3>;

  /** Mean absolute difference between channel a of image A and channel b of image B. */
  static double MeanDifference(const DwiType* imageA, unsigned int a, const DwiType* imageB, unsigned int b)
  {
    itk::ImageRegionConstIterator<DwiType> itA(imageA, imageA->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> itB(imageB, imageB->GetLargestPossibleRegion());
    double sum = 0;
    for (; !itA.IsAtEnd(); ++itA, ++itB)
      sum += std::abs(static_cast<double>(itA.Get()[a]) - static_cast<double>(itB.Get()[b]));
    return sum / imageA->GetLargestPossibleRegion().GetNumberOfPixels();
  }

public:
  void Correct_ReferenceNotFirst_KeepsVolumeOrder()
  {
    // A motion-free DWI with two weighted volumes and the b0 volume last, as Philips scanners store
    // it; the b0 volume is the reference. The volumes are scaled copies of a brain FA map, reduced
    // to 3 mm to keep the two affine registrations short. The test data's DWIs are too small for the
    // registration (3 slices or 11 voxels per axis).
    auto fa = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("mean_FA.nii.gz"));
    FloatImageType::Pointer itkFa;
    mitk::CastToItkImage(fa, itkFa);
    auto shrinker = itk::ShrinkImageFilter<FloatImageType, FloatImageType>::New();
    shrinker->SetInput(itkFa);
    shrinker->SetShrinkFactors(3);
    shrinker->Update();
    auto volume = shrinker->GetOutput();

    const std::vector<double> scales = {400, 600, 1000};
    const unsigned int reference = 2;
    auto dwi = DwiType::New();
    dwi->CopyInformation(volume);
    dwi->SetRegions(volume->GetLargestPossibleRegion());
    dwi->SetVectorLength(scales.size());
    dwi->Allocate();
    itk::ImageRegionConstIterator<FloatImageType> in(volume, volume->GetLargestPossibleRegion());
    itk::ImageRegionIterator<DwiType> out(dwi, dwi->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
    {
      auto pixel = out.Get();
      for (unsigned int c = 0; c < scales.size(); ++c)
        pixel[c] = static_cast<short>(scales[c] * in.Get());
      out.Set(pixel);
    }

    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    PropHelper::GradientDirectionType g;
    g.fill(0.0);
    g[0] = 1;
    gradients->push_back(g);
    g.fill(0.0);
    g[1] = 1;
    gradients->push_back(g);
    g.fill(0.0);
    gradients->push_back(g);

    auto image = mitk::GrabItkImageMemory(dwi.GetPointer());
    PropHelper::MeasurementFrameType identity;
    identity.set_identity();
    PropHelper::InitializeImage(image, gradients.GetPointer(), 1000, identity, false);

    auto filter = mitk::DWIHeadMotionCorrectionFilter::New();
    filter->SetInput(image);
    filter->Update();
    auto corrected = filter->GetCorrectedImage();
    auto correctedDwi = PropHelper::GetItkVectorImage(corrected);

    // The reference is copied unchanged to its own index.
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(scales.size()), correctedDwi->GetVectorLength());
    CPPUNIT_ASSERT_EQUAL(0.0, MeanDifference(correctedDwi, reference, dwi, reference));

    // Without motion the registered weighted volumes stay close to their inputs, much closer than
    // to the b0 volume, and their gradients barely turn.
    auto correctedGradients = PropHelper::GetGradientContainer(corrected);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, correctedGradients->ElementAt(reference).two_norm(), 1e-12);
    for (unsigned int c = 0; c < scales.size(); ++c)
    {
      if (c == reference)
        continue;
      CPPUNIT_ASSERT(MeanDifference(correctedDwi, c, dwi, c) < 0.5 * MeanDifference(correctedDwi, c, dwi, reference));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (correctedGradients->ElementAt(c) - gradients->ElementAt(c)).two_norm(), 0.1);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionHeadMotionCorrectionFilter)
