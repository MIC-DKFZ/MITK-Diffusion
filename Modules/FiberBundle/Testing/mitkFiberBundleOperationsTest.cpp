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
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

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
  MITK_TEST(SubsampleFibers_WithoutRandomSeed_IsReproducible);
  MITK_TEST(ResampleToNumPoints_FibersWithoutLength_Terminate);
  MITK_TEST(DerivedBundles_KeepRasFlagAndHeader);
  MITK_TEST(ApplyCurvatureThreshold_KeepsTwoPointFibersAndNoSinglePoints);
  MITK_TEST(TransformFibers_NoTransform_Throws);
  MITK_TEST(SetFiberPolyData_LeavesOldPolyDataAndPointDataBehind);
  MITK_TEST(SetSingleFiberColor_UsesRgbOrder);
  MITK_TEST(SetFiberColors_ChecksTheArray);
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

  /** Bundle of the given fibers, built point by point */
  static mitk::FiberBundle::Pointer MakeBundle(const std::vector<std::vector<std::array<double, 3>>> &fibers)
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const auto &fiber : fibers)
    {
      auto line = vtkSmartPointer<vtkPolyLine>::New();
      for (const auto &p : fiber)
        line->GetPointIds()->InsertNextId(points->InsertNextPoint(p.data()));
      lines->InsertNextCell(line);
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    return mitk::FiberBundle::New(polyData);
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

    fib = source->GetDeepCopy();
    CPPUNIT_ASSERT(fib->ApplyCurvatureThreshold(5.0f, false));
    AssertWeightsFollowFibers(source, fib, "ApplyCurvatureThreshold");

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

  void SubsampleFibers_WithoutRandomSeed_IsReproducible()
  {
    auto source = LoadWeightedFibers();
    auto first = source->SubsampleFibers(0.5f, false);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(std::round(source->GetNumFibers() * 0.5f)), first->GetNumFibers());
    AssertIdentical(first, source->SubsampleFibers(0.5f, false), "second subsample");
  }

  void ResampleToNumPoints_FibersWithoutLength_Terminate()
  {
    auto fib = MakeBundle({{{0, 0, 0}, {2, 0, 0}, {5, 0, 0}, {7, 0, 0}, {10, 0, 0}},
                           {{5, 5, 5}, {5, 5, 5}, {5, 5, 5}},
                           {{7, 7, 7}}});
    CPPUNIT_ASSERT_EQUAL(3u, fib->GetNumFibers());
    fib->ResampleToNumPoints(4);
    CPPUNIT_ASSERT_EQUAL(3u, fib->GetNumFibers());
    bool foundLine = false;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(4), cell->GetNumberOfPoints());
      const double *first = cell->GetPoints()->GetPoint(0);
      if (first[0] == 0.0 && first[1] == 0.0)
      {
        foundLine = true;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, cell->GetPoints()->GetPoint(3)[0], 1e-4);
      }
    }
    CPPUNIT_ASSERT_MESSAGE("The line keeps its end points", foundLine);
  }

  void DerivedBundles_KeepRasFlagAndHeader()
  {
    auto source = LoadWeightedFibers();
    source->setIsRAS(true);
    auto header = source->GetTrackVisHeader();
    header.dim[0] = 77;
    header.voxel_size[0] = 1.25f;
    source->SetTrackVisHeader(header);
    auto other = LoadWeightedFibers();
    auto shortFibers = source->GetDeepCopy();
    shortFibers->RemoveLongFibers(source->GetMedianFiberLength());

    auto check = [](mitk::FiberBundle *fib, const std::string &operation)
    {
      CPPUNIT_ASSERT_MESSAGE(operation + " keeps the RAS flag", fib->IsRAS());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 77, static_cast<int>(fib->GetTrackVisHeader().dim[0]));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 1.25f, fib->GetTrackVisHeader().voxel_size[0]);
    };
    check(source->GetDeepCopy(), "GetDeepCopy");
    check(source->AddBundle(other), "AddBundle");
    check(source->AddBundles({other}), "AddBundles");
    check(source->SubtractBundle(shortFibers), "SubtractBundle");
    check(source->SubtractBundle(source), "SubtractBundle without remaining fibers");
    check(source->FilterByWeights(1.0f), "FilterByWeights");
    check(source->SubsampleFibers(0.5f, false), "SubsampleFibers");
  }

  void ApplyCurvatureThreshold_KeepsTwoPointFibersAndNoSinglePoints()
  {
    // A two-point fiber, a straight fiber and a zigzag whose every corner is sharper than 5 mm
    auto fib = MakeBundle({{{0, 0, 10}, {3, 0, 10}},
                           {{0, 5, 0}, {2, 5, 0}, {4, 5, 0}, {6, 5, 0}, {8, 5, 0}},
                           {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 2, 0}, {1, 2, 0}, {1, 3, 0}}});
    CPPUNIT_ASSERT(fib->ApplyCurvatureThreshold(5.0f, false));

    bool twoPointFiber = false;
    bool straightFiber = false;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_MESSAGE("No fiber with a single point", cell->GetNumberOfPoints() > 1);
      const double *first = cell->GetPoints()->GetPoint(0);
      twoPointFiber = twoPointFiber || (first[2] == 10.0 && cell->GetNumberOfPoints() == 2);
      straightFiber = straightFiber || (first[1] == 5.0 && cell->GetNumberOfPoints() == 5);
    }
    CPPUNIT_ASSERT_MESSAGE("The two-point fiber is kept", twoPointFiber);
    CPPUNIT_ASSERT_MESSAGE("The straight fiber is kept whole", straightFiber);
  }

  void TransformFibers_NoTransform_Throws()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    CPPUNIT_ASSERT_THROW(fib->TransformFibers(itk::ScalableAffineTransform<mitk::ScalarType>::Pointer()), mitk::Exception);

    auto transform = itk::ScalableAffineTransform<mitk::ScalarType>::New();
    auto translation = transform->GetTranslation();
    translation[0] = 5.0;
    transform->SetTranslation(translation);
    const double x = fib->GetFiberPolyData()->GetPoint(0)[0];
    fib->TransformFibers(transform);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x + 5.0, fib->GetFiberPolyData()->GetPoint(0)[0], 1e-4);
  }

  void SetFiberPolyData_LeavesOldPolyDataAndPointDataBehind()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    auto oldPolyData = fib->GetFiberPolyData();
    const auto oldNumPoints = oldPolyData->GetNumberOfPoints();
    auto stale = vtkSmartPointer<vtkFloatArray>::New();
    stale->SetName("STALE");
    stale->SetNumberOfValues(oldNumPoints);
    oldPolyData->GetPointData()->AddArray(stale);

    auto other = LoadFibers("fiberBundleX.fib");
    fib->SetFiberPolyData(other->GetFiberPolyData());

    CPPUNIT_ASSERT(fib->Equals(other));
    CPPUNIT_ASSERT_MESSAGE("Point data of the old fibers is not carried over",
                           !fib->GetFiberPolyData()->GetPointData()->HasArray("STALE"));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("The old polydata is not changed", oldNumPoints, oldPolyData->GetNumberOfPoints());
  }

  void SetSingleFiberColor_UsesRgbOrder()
  {
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}}, {{0, 5, 0}, {1, 5, 0}, {2, 5, 0}}});
    fib->SetSingleFiberColor(10, 20, 30, 1);
    auto colours = fib->GetFiberColors();
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
      {
        const vtkIdType id = cell->GetPointId(j);
        CPPUNIT_ASSERT_EQUAL(i == 1 ? 10.0 : 255.0, colours->GetComponent(id, 0));
        CPPUNIT_ASSERT_EQUAL(i == 1 ? 20.0 : 255.0, colours->GetComponent(id, 1));
        CPPUNIT_ASSERT_EQUAL(i == 1 ? 30.0 : 255.0, colours->GetComponent(id, 2));
      }
    }
  }

  void SetFiberColors_ChecksTheArray()
  {
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}});
    auto rgb = vtkSmartPointer<vtkUnsignedCharArray>::New();
    rgb->SetNumberOfComponents(3);
    for (int i = 0; i < 3; ++i)
    {
      const unsigned char value[3] = {10, 20, 30};
      rgb->InsertNextTypedTuple(value);
    }
    fib->SetFiberColors(rgb);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("RGB gets an opaque alpha", 255.0, fib->GetFiberColors()->GetComponent(2, 3));
    CPPUNIT_ASSERT_EQUAL(20.0, fib->GetFiberColors()->GetComponent(2, 1));

    auto tooShort = vtkSmartPointer<vtkUnsignedCharArray>::New();
    tooShort->SetNumberOfComponents(4);
    tooShort->SetNumberOfTuples(2);
    CPPUNIT_ASSERT_THROW(fib->SetFiberColors(tooShort), mitk::Exception);
    CPPUNIT_ASSERT_THROW(fib->SetFiberColors(vtkSmartPointer<vtkUnsignedCharArray>()), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleOperations)
