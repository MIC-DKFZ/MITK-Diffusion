/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>

class mitkDiffusionPropertyHelperTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionPropertyHelperTestSuite);
  MITK_TEST(GetMeasurementFrame_Missing_ReturnsIdentity);
  MITK_TEST(InitializeImage_NoMeasurementFrame_KeepsGradients);
  MITK_TEST(Getters_NoDiffusionProperties_Throw);
  MITK_TEST(AverageRedundantGradients_NegativeValues_AveragesSigned);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;

  mitk::Image::Pointer m_Dwi;

  static PropHelper::GradientDirectionType Vec(double x, double y, double z)
  {
    PropHelper::GradientDirectionType v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
  }

  /** One voxel, one component per gradient. */
  static mitk::Image::Pointer CreateDwi(const std::vector<short>& values, const std::vector<PropHelper::GradientDirectionType>& gradients)
  {
    auto itkImage = PropHelper::ImageType::New();
    PropHelper::ImageType::RegionType region;
    region.SetSize({1, 1, 1});
    itkImage->SetRegions(region);
    itkImage->SetVectorLength(values.size());
    itkImage->Allocate();
    PropHelper::ImageType::PixelType pixel(values.size());
    for (unsigned int i = 0; i < values.size(); ++i)
      pixel[i] = values[i];
    itkImage->SetPixel({0, 0, 0}, pixel);

    auto image = mitk::GrabItkImageMemory(itkImage);
    auto container = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : gradients)
      container->push_back(g);
    PropHelper::SetOriginalGradientContainer(image, container);
    PropHelper::SetReferenceBValue(image, 1000);
    PropHelper::InitializeImage(image);
    return image;
  }

public:
  void setUp() override
  {
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void GetMeasurementFrame_Missing_ReturnsIdentity()
  {
    m_Dwi->GetPropertyList()->RemoveProperty(PropHelper::GetMeasurementFramePropertyName());
    PropHelper::MeasurementFrameType identity;
    identity.set_identity();
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(m_Dwi) == identity);
  }

  void InitializeImage_NoMeasurementFrame_KeepsGradients()
  {
    // Filters set only the modified gradients and no measurement frame before InitializeImage;
    // it then derives the original gradients while reading the (missing) frame.
    auto expected = PropHelper::GetGradientContainer(m_Dwi);
    auto image = m_Dwi->Clone();
    image->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    image->GetPropertyList()->RemoveProperty(PropHelper::GetMeasurementFramePropertyName());
    PropHelper::InitializeImage(image);

    auto original = PropHelper::GetOriginalGradientContainer(image);
    CPPUNIT_ASSERT_EQUAL(expected->Size(), original->Size());
    for (unsigned int i = 0; i < expected->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (expected->ElementAt(i) - original->ElementAt(i)).two_norm(), 1e-12);
  }

  void Getters_NoDiffusionProperties_Throw()
  {
    auto image = mitk::Image::New();
    image->Initialize(m_Dwi);
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT_THROW(PropHelper::GetGradientContainer(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetOriginalGradientContainer(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetReferenceBValue(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetBValueMap(image), mitk::Exception);
  }

  void AverageRedundantGradients_NegativeValues_AveragesSigned()
  {
    auto image = CreateDwi({-10, -20, 5}, {Vec(1, 0, 0), Vec(1, 0, 0), Vec(0, 1, 0)});
    PropHelper::AverageRedundantGradients(image, 0.001);

    CPPUNIT_ASSERT_EQUAL(size_t(2), image->GetPixelType().GetNumberOfComponents());
    CPPUNIT_ASSERT_EQUAL(2u, PropHelper::GetGradientContainer(image)->Size());
    auto itkImage = PropHelper::GetItkVectorImage(image);
    auto pixel = itkImage->GetPixel({0, 0, 0});
    CPPUNIT_ASSERT_EQUAL(short(-15), pixel[0]);
    CPPUNIT_ASSERT_EQUAL(short(5), pixel[1]);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionPropertyHelper)
