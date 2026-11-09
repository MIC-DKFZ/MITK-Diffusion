/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionDICOMFileReader.h>
#include <mitkDiffusionHeaderGEDICOMFileReader.h>
#include <mitkDiffusionHeaderPhilipsDICOMFileReader.h>
#include <mitkDiffusionHeaderSiemensDICOMFileReader.h>
#include <mitkDiffusionHeaderSiemensMosaicDICOMFileReader.h>

/** Parts of the diffusion DICOM import that work without DICOM files; mitkDiffusionDicomSeriesTest reads DICOM series. */
class mitkDiffusionDicomHeaderTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionDicomHeaderTestSuite);
  MITK_TEST(CreateHeaderReader_Vendors_MatchingReader);
  MITK_TEST(RevealBinaryTagC_WrongSize_CopiesNothing);
  CPPUNIT_TEST_SUITE_END();

  template <class TReader>
  static bool Creates(const std::string& vendor, const std::string& imageType, bool expectMosaic)
  {
    bool isMosaic = !expectMosaic;
    auto reader = mitk::DiffusionDICOMFileReader::CreateHeaderReader(vendor, imageType, isMosaic);
    return reader.IsNotNull() && typeid(*reader) == typeid(TReader) && isMosaic == expectMosaic;
  }

public:
  void CreateHeaderReader_Vendors_MatchingReader()
  {
    const std::string normal = "ORIGINAL\\PRIMARY\\DIFFUSION\\NONE";
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderSiemensDICOMFileReader>("SIEMENS", normal, false));
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderSiemensMosaicDICOMFileReader>("SIEMENS", "ORIGINAL\\PRIMARY\\DIFFUSION\\NONE\\ND\\MOSAIC", true));
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderSiemensDICOMFileReader>("Siemens HealthCare GmbH", normal, false));
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderSiemensDICOMFileReader>("Siemens Healthineers", normal, false));
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderGEDICOMFileReader>("GE MEDICAL SYSTEMS", normal, false));
    CPPUNIT_ASSERT(Creates<mitk::DiffusionHeaderPhilipsDICOMFileReader>("Philips Medical Systems", normal, false));

    bool isMosaic = true;
    CPPUNIT_ASSERT(mitk::DiffusionDICOMFileReader::CreateHeaderReader("Other Vendor", normal, isMosaic).IsNull());
    CPPUNIT_ASSERT(!isMosaic);
  }

  void RevealBinaryTagC_WrongSize_CopiesNothing()
  {
    const gdcm::Tag bvalueTag(0x0018, 0x9087);
    const double values[2] = {1000, 2000};
    gdcm::DataElement element(bvalueTag);
    element.SetByteValue(reinterpret_cast<const char*>(values), sizeof(values));
    gdcm::DataSet dataset;
    dataset.Insert(element);

    double bvalue = -1;
    CPPUNIT_ASSERT(!mitk::RevealBinaryTagC(bvalueTag, dataset, reinterpret_cast<char*>(&bvalue), sizeof(bvalue)));
    CPPUNIT_ASSERT_EQUAL(-1.0, bvalue);

    double both[2] = {0, 0};
    CPPUNIT_ASSERT(mitk::RevealBinaryTagC(bvalueTag, dataset, reinterpret_cast<char*>(both), sizeof(both)));
    CPPUNIT_ASSERT_EQUAL(2000.0, both[1]);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionDicomHeader)
