/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <set>
#include <string>
#include <vector>

/** Regression tests for the operations of mitk::FiberBundle (selection, removal, combination). */
class mitkFiberBundleOperationsTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleOperationsTestSuite);
  MITK_TEST(ExtractFiberIdSubset_NotWithOverlappingChildren_EqualsSetDifference);
  CPPUNIT_TEST_SUITE_END();

  mitk::FiberBundle::Pointer LoadFibers(const std::string &relativePath)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/" + relativePath));
  }

public:
  void ExtractFiberIdSubset_NotWithOverlappingChildren_EqualsSetDifference()
  {
    auto fib = LoadFibers("fiberBundleX.fib");

    auto storage = mitk::StandaloneDataStorage::New();
    auto composite = mitk::PlanarFigureComposite::New();
    composite->setOperationType(mitk::PlanarFigureComposite::NOT);
    auto compositeNode = mitk::DataNode::New();
    compositeNode->SetData(composite);
    storage->Add(compositeNode);

    // ROI1 twice: the second child overlaps the first completely, which made the difference
    // larger than the buffer it was written to.
    std::set<unsigned int> excluded;
    for (const std::string name : {"ROI1.pf", "ROI1.pf", "ROI2.pf", "ROI3.pf"})
    {
      auto node = mitk::DataNode::New();
      node->SetData(mitk::IOUtil::Load(GetTestDataFilePath("DiffusionImaging/" + name)).at(0));
      storage->Add(node, compositeNode);
      auto ids = fib->ExtractFiberIdSubset(node, storage);
      excluded.insert(ids.begin(), ids.end());
    }
    CPPUNIT_ASSERT_MESSAGE("The ROIs hit fibers", !excluded.empty());

    std::vector<unsigned int> expected;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      if (excluded.count(i) == 0)
        expected.push_back(i);

    CPPUNIT_ASSERT_MESSAGE("NOT keeps exactly the fibers no child hits",
                           fib->ExtractFiberIdSubset(compositeNode, storage) == expected);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleOperations)
