/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkIOUtil.h"
#include "mitkTestingMacros.h"
#include "mitkTestFixture.h"
#include "itkNonLocalMeansDenoisingFilter.h"
#include "mitkITKImageImport.h"
#include <mitkImageCast.h>
#include <mitkDiffusionPropertyHelper.h>
#include <itkImageRegionConstIterator.h>

class mitkNonLocalMeansDenoisingTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkNonLocalMeansDenoisingTestSuite);
  MITK_TEST(Denoise_NLMg_shouldReturnTrue);
  MITK_TEST(Denoise_NLMr_shouldReturnTrue);
  MITK_TEST(Denoise_NLMv_shouldReturnTrue);
  MITK_TEST(Denoise_NLMvr_shouldReturnTrue);
  MITK_TEST(Denoise_SeveralWorkUnits_SameResultAndCountsEveryVoxel);
  MITK_TEST(Denoise_LargeDifferences_KeepsValues);
  MITK_TEST(Denoise_ZeroVariance_Throws);
  MITK_TEST(Denoise_EmptyMask_ReturnsZeros);
  MITK_TEST(Denoise_MaskOfOtherSize_Throws);
  MITK_TEST(Denoise_NewInputWithoutMask_DenoisesEveryVoxel);
  CPPUNIT_TEST_SUITE_END();

private:

  typedef itk::VectorImage<short,3> VectorImagetType;

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/
  mitk::Image::Pointer m_Image;
  mitk::Image::Pointer m_ReferenceImage;
  mitk::Image::Pointer m_DenoisedImage;
  itk::Image<short, 3>::Pointer m_ImageMask;
  itk::NonLocalMeansDenoisingFilter<short>::Pointer m_DenoisingFilter;

  /** A row of voxels with one channel. */
  static VectorImagetType::Pointer CreateRow(const std::vector<short>& values)
  {
    auto image = VectorImagetType::New();
    VectorImagetType::RegionType region;
    region.SetSize({values.size(), 1, 1});
    image->SetRegions(region);
    image->SetVectorLength(1);
    image->Allocate();
    VectorImagetType::PixelType pixel(1);
    for (unsigned int x = 0; x < values.size(); ++x)
    {
      pixel[0] = values[x];
      image->SetPixel({x, 0, 0}, pixel);
    }
    return image;
  }

public:

  /**
* @brief Setup Always call this method before each Test-case to ensure correct and new intialization of the used members for a new test case. (If the members are not used in a test, the method does not need to be called).
*/
  void setUp() override
  {
    //generate test images
    std::string imagePath = GetTestDataFilePath("DiffusionImaging/Denoising/test_multi.dwi");

    m_Image =  mitk::IOUtil::Load<mitk::Image>(imagePath);
    // Also makes the driver link MitkDiffusionImage, which --as-needed drops otherwise; its
    // activator registers the DWI reader that the line above needs.
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(m_Image));
    m_ReferenceImage = nullptr;
    m_DenoisedImage = mitk::Image::New();

    //initialise Filter
    m_DenoisingFilter = itk::NonLocalMeansDenoisingFilter<short>::New();
    VectorImagetType::Pointer vectorImage;
    mitk::CastToItkImage(m_Image,vectorImage);
    m_DenoisingFilter->SetInputImage(vectorImage);
    m_DenoisingFilter->SetNumberOfWorkUnits(1);
    m_DenoisingFilter->SetComparisonRadius(1);
    m_DenoisingFilter->SetSearchRadius(1);
    m_DenoisingFilter->SetVariance(500);
  }

  void tearDown() override
  {
    m_Image = nullptr;
    m_ImageMask = nullptr;
    m_ReferenceImage = nullptr;
    m_DenoisingFilter = nullptr;
    m_DenoisedImage = nullptr;
  }

  void Denoise_NLMg_shouldReturnTrue()
  {
    std::string referenceImagePath = GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMg.dwi");
    m_ReferenceImage =  mitk::IOUtil::Load<mitk::Image>(referenceImagePath);

    m_DenoisingFilter->SetUseRicianAdaption(false);
    m_DenoisingFilter->SetUseJointInformation(false);
    try
    {
      m_DenoisingFilter->Update();
    }
    catch(std::exception& e)
    {
      MITK_ERROR << e.what();
    }

    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(),m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());

    MITK_ASSERT_EQUAL( m_DenoisedImage, m_ReferenceImage, "NLMg should always return the same result.");
  }

  void Denoise_NLMr_shouldReturnTrue()
  {
    std::string referenceImagePath = GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMr.dwi");
    m_ReferenceImage = mitk::IOUtil::Load<mitk::Image>(referenceImagePath);

    m_DenoisingFilter->SetUseRicianAdaption(true);
    m_DenoisingFilter->SetUseJointInformation(false);
    try
    {
      m_DenoisingFilter->Update();
    }
    catch(std::exception& e)
    {
      MITK_ERROR << e.what();
    }

    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(),m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());

    MITK_ASSERT_EQUAL( m_DenoisedImage, m_ReferenceImage, "NLMr should always return the same result.");
  }

  void Denoise_NLMv_shouldReturnTrue()
  {
    std::string referenceImagePath = GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMv.dwi");
    m_ReferenceImage = mitk::IOUtil::Load<mitk::Image>(referenceImagePath);
    m_DenoisingFilter->SetUseRicianAdaption(false);
    m_DenoisingFilter->SetUseJointInformation(true);
    try
    {
      m_DenoisingFilter->Update();
    }
    catch(std::exception& e)
    {
      MITK_ERROR << e.what();
    }

    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(),m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());

    MITK_ASSERT_EQUAL( m_DenoisedImage, m_ReferenceImage, "NLMv should always return the same result.");
  }

  void Denoise_NLMvr_shouldReturnTrue()
  {
    std::string referenceImagePath = GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMvr.dwi");
    m_ReferenceImage = mitk::IOUtil::Load<mitk::Image>(referenceImagePath);

    m_DenoisingFilter->SetUseRicianAdaption(true);
    m_DenoisingFilter->SetUseJointInformation(true);
    try
    {
      m_DenoisingFilter->Update();
    }
    catch(std::exception& e)
    {
      MITK_ERROR << e.what();
    }

    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(),m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());

    MITK_ASSERT_EQUAL( m_DenoisedImage, m_ReferenceImage, "NLMvr should always return the same result.");
  }

  void Denoise_SeveralWorkUnits_SameResultAndCountsEveryVoxel()
  {
    m_ReferenceImage = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMg.dwi"));

    // All work units increment the progress counter.
    m_DenoisingFilter->SetNumberOfWorkUnits(4);
    m_DenoisingFilter->SetUseRicianAdaption(false);
    m_DenoisingFilter->SetUseJointInformation(false);
    m_DenoisingFilter->Update();

    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(m_DenoisingFilter->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels()), m_DenoisingFilter->GetCurrentVoxelCount());
    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(), m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());
    MITK_ASSERT_EQUAL(m_DenoisedImage, m_ReferenceImage, "NLMg with several work units should return the single-unit result.");
  }

  void Denoise_LargeDifferences_KeepsValues()
  {
    // Neighbours differ by 60000, so their weights are exp(-3.6e9) = 0 and every voxel keeps its
    // value (negative results are clipped to 0); the squared difference does not fit an int.
    auto filter = itk::NonLocalMeansDenoisingFilter<short>::New();
    filter->SetInputImage(CreateRow({-30000, 30000, -30000}));
    filter->SetSearchRadius(1);
    filter->SetComparisonRadius(0);
    filter->SetVariance(1);
    filter->Update();
    const std::vector<short> expected = {0, 30000, 0};
    for (unsigned int x = 0; x < expected.size(); ++x)
      CPPUNIT_ASSERT_EQUAL(expected[x], filter->GetOutput()->GetPixel({x, 0, 0})[0]);
  }

  void Denoise_ZeroVariance_Throws()
  {
    m_DenoisingFilter->SetVariance(0);
    CPPUNIT_ASSERT_THROW(m_DenoisingFilter->Update(), itk::ExceptionObject);
  }

  void Denoise_EmptyMask_ReturnsZeros()
  {
    VectorImagetType::Pointer vectorImage;
    mitk::CastToItkImage(m_Image, vectorImage);
    auto mask = itk::Image<short, 3>::New();
    mask->SetRegions(vectorImage->GetLargestPossibleRegion());
    mask->Allocate();
    mask->FillBuffer(0);
    m_DenoisingFilter->SetInputMask(mask);
    m_DenoisingFilter->Update();

    itk::ImageRegionConstIterator<VectorImagetType> it(m_DenoisingFilter->GetOutput(), m_DenoisingFilter->GetOutput()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      for (unsigned int c = 0; c < it.Get().GetSize(); ++c)
        CPPUNIT_ASSERT_EQUAL(short(0), it.Get()[c]);
  }

  void Denoise_MaskOfOtherSize_Throws()
  {
    auto mask = itk::Image<short, 3>::New();
    itk::Image<short, 3>::RegionType region;
    region.SetSize({2, 2, 2});
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(1);
    m_DenoisingFilter->SetInputMask(mask);
    CPPUNIT_ASSERT_THROW(m_DenoisingFilter->Update(), itk::ExceptionObject);
  }

  void Denoise_NewInputWithoutMask_DenoisesEveryVoxel()
  {
    // The full mask made for the first, smaller input must not be reused for the second.
    m_DenoisingFilter->SetInputImage(CreateRow({100, 100, 100}));
    m_DenoisingFilter->Update();

    VectorImagetType::Pointer vectorImage;
    mitk::CastToItkImage(m_Image, vectorImage);
    m_DenoisingFilter->SetInputImage(vectorImage);
    m_DenoisingFilter->UpdateLargestPossibleRegion();

    m_ReferenceImage = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/Denoising/test_multi_NLMg.dwi"));
    mitk::GrabItkImageMemory(m_DenoisingFilter->GetOutput(), m_DenoisedImage);
    m_DenoisedImage->SetPropertyList(m_Image->GetPropertyList()->Clone());
    MITK_ASSERT_EQUAL(m_DenoisedImage, m_ReferenceImage, "NLMg after a run on another input should return the reference.");
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkNonLocalMeansDenoising)
