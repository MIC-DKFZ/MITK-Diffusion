/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkIOMimeTypes.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itksys/SystemTools.hxx>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <set>
#include <string>

/** Regression tests for the fiber bundle readers and writers. Files are written to a temporary
 *  directory, nothing is added to MITK-Data. */
class mitkFiberBundleIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleIOTestSuite);
  MITK_TEST(VtkWriter_AllOptionCombinations_LeaveBundleUnchanged);
  MITK_TEST(DicomWriter_ReferencedImageUids_TractogramGetsOwnUids);
  CPPUNIT_TEST_SUITE_END();

  std::string m_TempDir;

  mitk::FiberBundle::Pointer LoadFibers(const std::string &relativePath)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/" + relativePath));
  }

  /** original.fib with a different weight per fiber */
  mitk::FiberBundle::Pointer LoadWeightedFibers()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      fib->SetFiberWeight(i, 0.5f + static_cast<float>(i) / fib->GetNumFibers());
    return fib;
  }

  static std::set<std::string> ArrayNames(vtkPolyData *polyData)
  {
    std::set<std::string> names;
    for (int i = 0; i < polyData->GetPointData()->GetNumberOfArrays(); ++i)
      names.insert(std::string("point ") + polyData->GetPointData()->GetArrayName(i));
    for (int i = 0; i < polyData->GetCellData()->GetNumberOfArrays(); ++i)
      names.insert(std::string("cell ") + polyData->GetCellData()->GetArrayName(i));
    return names;
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("FiberBundleIOTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void VtkWriter_AllOptionCombinations_LeaveBundleUnchanged()
  {
    auto fib = LoadWeightedFibers();
    const auto arraysBefore = ArrayNames(fib->GetFiberPolyData());
    const std::string path = m_TempDir + "/fibers.fib";

    for (bool saveWeights : {true, false})
    {
      for (bool saveColors : {true, false})
      {
        mitk::IFileWriter::Options options;
        options["Save fiber weights"] = saveWeights;
        options["Save color information"] = saveColors;
        mitk::IOUtil::Save(fib, path, options);

        CPPUNIT_ASSERT_MESSAGE("Saving does not change the arrays of the bundle",
                               arraysBefore == ArrayNames(fib->GetFiberPolyData()));

        auto reloaded = mitk::IOUtil::Load<mitk::FiberBundle>(path);
        CPPUNIT_ASSERT(reloaded->Equals(fib));
        for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
        {
          const float expected = saveWeights ? fib->GetFiberWeight(i) : 1.0f;
          CPPUNIT_ASSERT_EQUAL(expected, reloaded->GetFiberWeight(i));
        }
      }
    }
  }

  static std::string StringProperty(const mitk::FiberBundle *fib, const std::string &name)
  {
    std::string value;
    fib->GetPropertyList()->GetStringProperty(name.c_str(), value);
    return value;
  }

  void DicomWriter_ReferencedImageUids_TractogramGetsOwnUids()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    const std::string mrImageStorage = "1.2.840.10008.5.1.4.1.1.4";
    fib->GetPropertyList()->SetStringProperty("DICOM.patient_id", "TEST0001");
    fib->GetPropertyList()->SetStringProperty("DICOM.patient_name", "Test^Tractogram");
    fib->GetPropertyList()->SetStringProperty("DICOM.study_instance_uid", "1.2.826.0.1.3680043.2.1125.1");
    fib->GetPropertyList()->SetStringProperty("DICOM.series_instance_uid", "1.2.826.0.1.3680043.2.1125.2");
    fib->GetPropertyList()->SetStringProperty("DICOM.sop_instance_uid", "1.2.826.0.1.3680043.2.1125.3");
    fib->GetPropertyList()->SetStringProperty("DICOM.sop_class_uid", mrImageStorage.c_str());
    fib->GetPropertyList()->SetStringProperty("DICOM.frame_of_reference_uid", "1.2.826.0.1.3680043.2.1125.4");

    // The mime type helper class of the module is not exported; this is its DICOM name.
    const std::string mimeType = mitk::IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.dcm";
    const std::string path1 = m_TempDir + "/tractogram1.dcm";
    mitk::IOUtil::Save(fib, mimeType, path1);
    auto first = mitk::IOUtil::Load<mitk::FiberBundle>(path1);
    CPPUNIT_ASSERT(first->Equals(fib, 0.001));

    // The properties name the referenced image, as before saving
    for (const std::string name : {"DICOM.patient_id", "DICOM.study_instance_uid", "DICOM.series_instance_uid",
                                   "DICOM.sop_instance_uid", "DICOM.sop_class_uid", "DICOM.frame_of_reference_uid"})
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, StringProperty(fib, name), StringProperty(first, name));

    // The tractogram's own instance is neither the image nor its SOP class
    const auto ownInstance = StringProperty(first, "DICOM.tractogram.sop_instance_uid");
    const auto ownSeries = StringProperty(first, "DICOM.tractogram.series_instance_uid");
    CPPUNIT_ASSERT(!ownInstance.empty() && ownInstance != "0");
    CPPUNIT_ASSERT(ownInstance != StringProperty(fib, "DICOM.sop_instance_uid"));
    CPPUNIT_ASSERT(ownInstance != "1.2.840.10008.5.1.4.1.1.66.6");
    CPPUNIT_ASSERT(ownSeries != StringProperty(fib, "DICOM.series_instance_uid"));

    // Saving the loaded bundle again refers to the same image under a new instance UID
    const std::string path2 = m_TempDir + "/tractogram2.dcm";
    mitk::IOUtil::Save(first, mimeType, path2);
    auto second = mitk::IOUtil::Load<mitk::FiberBundle>(path2);
    CPPUNIT_ASSERT_EQUAL(StringProperty(fib, "DICOM.sop_instance_uid"), StringProperty(second, "DICOM.sop_instance_uid"));
    CPPUNIT_ASSERT(StringProperty(second, "DICOM.tractogram.sop_instance_uid") != ownInstance);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleIO)
