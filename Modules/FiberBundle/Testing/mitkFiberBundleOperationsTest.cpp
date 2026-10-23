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

#include <itkMath.h>
#include <omp.h>
#include <vtkCell.h>
#include <vtkPoints.h>

#include <array>
#include <cmath>
#include <map>
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
  MITK_TEST(Weights_RemovingOrSplittingFibers_FollowTheirFibers);
  MITK_TEST(CompressAndSubtract_EightThreads_EqualOneThread);
  CPPUNIT_TEST_SUITE_END();

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

  using PointKey = std::array<long long, 3>;

  static PointKey Key(const double *p)
  {
    return {std::llround(p[0] * 1000), std::llround(p[1] * 1000), std::llround(p[2] * 1000)};
  }

  /** Fibers of a bundle by the positions of their points, to find the source of a kept or split fiber. */
  static std::map<PointKey, std::set<unsigned int>> FibersByPoint(mitk::FiberBundle *fib)
  {
    std::map<PointKey, std::set<unsigned int>> owners;
    auto polyData = fib->GetFiberPolyData();
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = polyData->GetCell(i);
      for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
        owners[Key(cell->GetPoints()->GetPoint(j))].insert(i);
    }
    return owners;
  }

  /** Every fiber of result carries the weight of the source fiber it comes from. */
  static void AssertWeightsFollowFibers(mitk::FiberBundle *source, mitk::FiberBundle *result, const std::string &operation)
  {
    CPPUNIT_ASSERT_MESSAGE(operation + " keeps fibers", result->GetNumFibers() > 0);
    auto owners = FibersByPoint(source);
    for (unsigned int i = 0; i < result->GetNumFibers(); ++i)
    {
      auto it = owners.find(Key(result->GetFiberPolyData()->GetCell(i)->GetPoints()->GetPoint(0)));
      CPPUNIT_ASSERT_MESSAGE(operation + ": fiber comes from the source", it != owners.end());
      bool found = false;
      for (auto owner : it->second)
        found = found || source->GetFiberWeight(owner) == result->GetFiberWeight(i);
      CPPUNIT_ASSERT_MESSAGE(operation + ": weight of fiber " + std::to_string(i), found);
    }
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

  void Weights_RemovingOrSplittingFibers_FollowTheirFibers()
  {
    auto source = LoadWeightedFibers();
    const float medianLength = source->GetMedianFiberLength();

    auto fib = source->GetDeepCopy();
    CPPUNIT_ASSERT(fib->RemoveShortFibers(medianLength));
    CPPUNIT_ASSERT(fib->GetNumFibers() < source->GetNumFibers());
    AssertWeightsFollowFibers(source, fib, "RemoveShortFibers");

    fib = source->GetDeepCopy();
    CPPUNIT_ASSERT(fib->RemoveLongFibers(medianLength));
    CPPUNIT_ASSERT(fib->GetNumFibers() < source->GetNumFibers());
    AssertWeightsFollowFibers(source, fib, "RemoveLongFibers");

    fib = source->GetDeepCopy();
    vnl_vector_fixed<double, 3> dir(0.0, 1.0, 0.0);
    fib->RemoveDir(dir, std::cos(5.0 * itk::Math::pi / 180));
    AssertWeightsFollowFibers(source, fib, "RemoveDir");

    auto shortFibers = source->GetDeepCopy();
    shortFibers->RemoveLongFibers(medianLength);
    auto subtracted = source->SubtractBundle(shortFibers);
    CPPUNIT_ASSERT(subtracted->GetNumFibers() < source->GetNumFibers());
    AssertWeightsFollowFibers(source, subtracted, "SubtractBundle");
  }

  /** Same fibers in the same order with the same points and weights */
  static void AssertIdentical(mitk::FiberBundle *expected, mitk::FiberBundle *actual, const std::string &message)
  {
    CPPUNIT_ASSERT_MESSAGE(message, expected->Equals(actual, 0.0));
    for (unsigned int i = 0; i < expected->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(message, expected->GetFiberWeight(i), actual->GetFiberWeight(i));
  }

  void CompressAndSubtract_EightThreads_EqualOneThread()
  {
    auto source = LoadWeightedFibers();
    auto shortFibers = source->GetDeepCopy();
    shortFibers->RemoveLongFibers(source->GetMedianFiberLength());

    const int threads = omp_get_max_threads();
    omp_set_num_threads(1);
    auto compressedSingle = source->GetDeepCopy();
    compressedSingle->Compress(0.1f);
    auto subtractedSingle = source->SubtractBundle(shortFibers);

    omp_set_num_threads(8);
    for (int run = 0; run < 3; ++run)
    {
      auto compressed = source->GetDeepCopy();
      compressed->Compress(0.1f);
      AssertIdentical(compressedSingle, compressed, "Compress, run " + std::to_string(run));
      AssertIdentical(subtractedSingle, source->SubtractBundle(shortFibers), "SubtractBundle, run " + std::to_string(run));
    }
    omp_set_num_threads(threads);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleOperations)
