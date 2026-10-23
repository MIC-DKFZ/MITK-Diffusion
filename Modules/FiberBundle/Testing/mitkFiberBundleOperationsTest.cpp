/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberBundle.h>
#include <mitkImage.h>
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
  MITK_TEST(RemoveFibersOutside_MaskMissesAllFibers_ReturnsEmptyBundle);
  MITK_TEST(ExtractFiberSubset_RoiIsNoPlanarFigure_Throws);
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

  void RemoveFibersOutside_MaskMissesAllFibers_ReturnsEmptyBundle()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");

    // A small mask far away from every fiber
    auto mask = mitk::FiberBundle::ItkUcharImgType::New();
    mitk::FiberBundle::ItkUcharImgType::RegionType region;
    region.SetSize({{2, 2, 2}});
    mask->SetRegions(region);
    mitk::FiberBundle::ItkUcharImgType::PointType origin;
    origin.Fill(10000.0);
    mask->SetOrigin(origin);
    mask->Allocate();
    mask->FillBuffer(1);

    auto inside = fib->RemoveFibersOutside(mask);
    CPPUNIT_ASSERT_MESSAGE("Nothing inside the mask gives an empty bundle", inside.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(0u, inside->GetNumFibers());

    auto outside = fib->RemoveFibersOutside(mask, true);
    CPPUNIT_ASSERT(outside.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(fib->GetNumFibers(), outside->GetNumFibers());
  }

  void ExtractFiberSubset_RoiIsNoPlanarFigure_Throws()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    auto storage = mitk::StandaloneDataStorage::New();
    auto node = mitk::DataNode::New();
    node->SetData(mitk::Image::New());
    CPPUNIT_ASSERT_THROW(fib->ExtractFiberSubset(node, storage), mitk::Exception);
    CPPUNIT_ASSERT_THROW(fib->ExtractFiberSubset(nullptr, storage), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleOperations)
