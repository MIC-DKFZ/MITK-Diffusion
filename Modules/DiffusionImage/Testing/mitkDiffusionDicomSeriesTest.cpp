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

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <array>
#include <cmath>

/**
 * Reads the synthetic DWI DICOM series of mitk-diffusion-test-data (DicomDwi/, one per vendor
 * parser; DicomDwi/README.md describes how they were made) with the DWI DICOM reader and checks
 * b-values, gradients, geometry and pixel values against the values written into the headers.
 */
class mitkDiffusionDicomSeriesTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionDicomSeriesTestSuite);
  MITK_TEST(SiemensMosaic_CsaHeader_GradientsAndSlices);
  MITK_TEST(Siemens_CsaHeader_GradientsAndSlices);
  MITK_TEST(GE_PrivateTags_GradientsAndSlices);
  MITK_TEST(Philips_PrivateAndStandardTags_GradientsAndSlices);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;

  // The four volumes of every series: b-value and unit gradient as written in the headers.
  struct Volume
  {
    double b;
    std::array<double, 3> g;
  };
  static constexpr std::array<Volume, 4> m_Volumes = {{{0, {0, 0, 0}},
                                                       {1000, {0.6, 0, -0.8}},
                                                       {1000, {0, 0.8, 0.6}},
                                                       {2000, {-0.48, 0.6, 0.64}}}};

  static mitk::Image::Pointer Load(const std::string& folder)
  {
    // The reader loads every DICOM file in the folder of the given file.
    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({"Diffusion Weighted Images"}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    const std::string file = mitk::GetDiffusionTestDataFilePath("DicomDwi/" + folder + "/" + (folder == "SiemensMosaic" ? "v0.dcm" : "v0_z0.dcm"));
    auto data = mitk::IOUtil::Load(file, &functor);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("one image per series", std::size_t(1), data.size());
    mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.front().GetPointer());
    CPPUNIT_ASSERT(image.IsNotNull());
    return image;
  }

  /** Checks the DWI properties and the pixel value 1000 v + 100 z + 10 y + x of every voxel and volume. */
  static void Check(const std::string& folder, unsigned int numSlices)
  {
    mitk::Image::Pointer image = Load(folder);
    CPPUNIT_ASSERT_MESSAGE(folder, PropHelper::IsDiffusionWeightedImage(image));
    // DICOM pads odd-length strings with a space, which is not part of the value.
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::string("synthetic DWI " + folder), image->GetProperty("name")->GetValueAsString());
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::string("Synthetic^DWI"), image->GetProperty("DICOM.patient_name")->GetValueAsString());

    // MITK scales the unit gradients with sqrt(b / reference b-value), the largest b-value.
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, 2000.0f, PropHelper::GetReferenceBValue(image));
    auto gradients = PropHelper::GetGradientContainer(image);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, m_Volumes.size(), std::size_t(gradients->Size()));
    for (unsigned int v = 0; v < m_Volumes.size(); ++v)
    {
      const double scale = std::sqrt(m_Volumes[v].b / 2000.0);
      for (unsigned int k = 0; k < 3; ++k)
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(folder + " volume " + std::to_string(v), scale * m_Volumes[v].g[k], gradients->ElementAt(v)[k], 1e-6);
    }
    const auto& bValues = PropHelper::GetBValueMap(image);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::size_t(3), bValues.size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::size_t(2), bValues.at(1000).size());

    auto itkImage = PropHelper::GetItkVectorImage(image);
    const auto size = itkImage->GetLargestPossibleRegion().GetSize();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::size_t(3), std::size_t(size[0]));
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::size_t(2), std::size_t(size[1]));
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, std::size_t(numSlices), std::size_t(size[2]));
    CPPUNIT_ASSERT_EQUAL_MESSAGE(folder, 4u, itkImage->GetNumberOfComponentsPerPixel());
    const auto spacing = itkImage->GetSpacing();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, spacing[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, spacing[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, spacing[2], 1e-6);

    for (unsigned int z = 0; z < numSlices; ++z)
      for (unsigned int y = 0; y < 2; ++y)
        for (unsigned int x = 0; x < 3; ++x)
        {
          const auto pixel = itkImage->GetPixel({{x, y, z}});
          for (unsigned int v = 0; v < m_Volumes.size(); ++v)
            CPPUNIT_ASSERT_EQUAL_MESSAGE(folder + " voxel " + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z),
                                         static_cast<int>(1000 * v + 100 * z + 10 * y + x), static_cast<int>(pixel[v]));
        }
  }

public:
  void SiemensMosaic_CsaHeader_GradientsAndSlices() { Check("SiemensMosaic", 4); }

  void Siemens_CsaHeader_GradientsAndSlices() { Check("Siemens", 2); }

  void GE_PrivateTags_GradientsAndSlices() { Check("GE", 2); }

  void Philips_PrivateAndStandardTags_GradientsAndSlices() { Check("Philips", 2); }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionDicomSeries)
