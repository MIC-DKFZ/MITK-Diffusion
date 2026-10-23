/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
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
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleIO)
