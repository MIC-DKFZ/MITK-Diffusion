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
#include <mitkITKImageImport.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkDiffusionTestData.h>
#include <mitkOmpThreadGuard.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkMath.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <omp.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/** Regression tests for the operations of mitk::FiberBundle (selection, removal, combination). */
class mitkFiberBundleOperationsTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleOperationsTestSuite);
  MITK_TEST(ExtractFiberIdSubset_NotWithOverlappingChildren_EqualsSetDifference);
  MITK_TEST(RemoveFibersOutside_MaskMissesAllFibers_ReturnsEmptyBundle);
  MITK_TEST(RemoveFibersOutside_SegmentsCrossingTheMask_CutAtTheMaskBorder);
  MITK_TEST(ExtractFiberSubset_RoiIsNoPlanarFigure_Throws);
  MITK_TEST(Weights_RemovingOrSplittingFibers_FollowTheirFibers);
  MITK_TEST(CompressAndSubtract_EightThreads_EqualOneThread);
  MITK_TEST(OmpThreadGuard_FailingCase_RestoresThreadCount);
  MITK_TEST(SubsampleFibers_WithoutRandomSeed_IsReproducible);
  MITK_TEST(ResampleToNumPoints_FibersWithoutLength_Terminate);
  MITK_TEST(DerivedBundles_KeepRasFlagAndHeader);
  MITK_TEST(ApplyCurvatureThreshold_KeepsTwoPointFibersAndNoSinglePoints);
  MITK_TEST(TransformFibers_NoTransform_Throws);
  MITK_TEST(SetFiberPolyData_LeavesOldPolyDataAndPointDataBehind);
  MITK_TEST(SetFiberPolyData_EqualFiberCount_NewFibersGetWeightOne);
  MITK_TEST(SetSingleFiberColor_UsesRgbOrder);
  MITK_TEST(SetFiberColors_ChecksTheArray);
  MITK_TEST(ColorFibersByOrientation_OnePointFibers_GetColours);
  MITK_TEST(ColorFibersByCurvature_UsesFirstSegment);
  MITK_TEST(ResampleLinear_NoStep_LeavesFibersUnchanged);
  MITK_TEST(PrintSelf_EmptyBundle_Prints);
  MITK_TEST(PrintStatistics_WithoutObjectAddresses);
  MITK_TEST(GeneratePolyDataByIds_InvalidId_Throws);
  MITK_TEST(ColorFibersByFiberWeights_EqualWeightsAboveOne_OpaqueColours);
  MITK_TEST(AddBundlesAndCompositeWithoutStorage_NoCrash);
  MITK_TEST(Clone_CopiesFibersAndMetadataIndependently);
  MITK_TEST(FiberWeights_WrongSizeOrIndex_Throw);
  MITK_TEST(LengthStatistics_OriginalFib_AsRecorded);
  MITK_TEST(ResampleToNumPoints_RealFibers_GiveNPointsAndKeepEndpoints);
  MITK_TEST(FilterByWeights_BothDirections_SplitAtThreshold);
  MITK_TEST(AddBundles_EqualsRepeatedAddBundle);
  MITK_TEST(RotateAroundAxis_EqualsTransformFibersWithoutTranslation);
  MITK_TEST(Clone_KeepsDataAndIsIndependent);
  MITK_TEST(EmptyBundle_EveryOperation_StaysEmpty);
  MITK_TEST(ResampleLinear_SyntheticFibers_EqualChordSpacingOnTheFiber);
  MITK_TEST(ResampleSpline_StraightFiber_PointsOnTheLine);
  MITK_TEST(ResampleToNumPoints_StraightFiber_EqualSpacing);
  MITK_TEST(ResampleToNumPoints_CurvedFibers_EqualArcLengthSpacing);
  MITK_TEST(Compress_SyntheticFibers_RemovesOnlyPointsWithinTheError);
  MITK_TEST(ApplyCurvatureThreshold_CircularArc_RadiusDecides);
  MITK_TEST(ColorFibersByOrientation_ThreePointFiber_SegmentDirections);
  MITK_TEST(ColorFibersByScalarMap_Mask_ColoursAndWeightsFromTheMask);
  MITK_TEST(TransformFibers_ItkTransform_MovesEveryPoint);
  MITK_TEST(ColorFibersByScalarMap_SignedIntegerImage_InterpolatedLikeFloat);
  MITK_TEST(TrackVisHeader_JoinedAndMovedBundles_ContainEveryPoint);
  MITK_TEST(Weights_OnePointFibers_KeepPlaceAndWeight);
  CPPUNIT_TEST_SUITE_END();

  mitk::FiberBundle::Pointer LoadFibers(const std::string &relativePath)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath(relativePath));
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
      node->SetData(mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(name)).at(0));
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

  /** Points of fiber i of fib */
  static std::vector<std::array<double, 3>> FiberPoints(mitk::FiberBundle *fib, unsigned int i)
  {
    vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
    std::vector<std::array<double, 3>> points;
    for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
    {
      double p[3];
      cell->GetPoints()->GetPoint(j, p);
      points.push_back({p[0], p[1], p[2]});
    }
    return points;
  }

  static void AssertPiece(mitk::FiberBundle *fib, unsigned int i, const std::array<double, 3> &first, const std::array<double, 3> &last,
                          std::size_t numPoints, float weight, const std::string &message)
  {
    const auto points = FiberPoints(fib, i);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": number of points", numPoints, points.size());
    for (int c = 0; c < 3; ++c)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message + ": first point", first[c], points.front()[c], 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message + ": last point", last[c], points.back()[c], 1e-5);
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": weight", weight, fib->GetFiberWeight(i));
  }

  // #1: the cut tested only the fiber points, so a segment that crosses the mask without a point
  // inside it was dropped by "remove outside" and kept whole by "remove inside", and the kept
  // pieces ended at their last point on the kept side instead of at the mask border.
  void RemoveFibersOutside_SegmentsCrossingTheMask_CutAtTheMaskBorder()
  {
    // 10 x 10 x 10 voxels of 1 mm, origin 0, set where the x index is 4 or 5: x from 3.5 to 5.5 mm
    auto mask = mitk::FiberBundle::ItkUcharImgType::New();
    mitk::FiberBundle::ItkUcharImgType::RegionType region;
    region.SetSize({{10, 10, 10}});
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(0);
    for (itk::IndexValueType y = 0; y < 10; ++y)
      for (itk::IndexValueType z = 0; z < 10; ++z)
      {
        mask->SetPixel({{4, y, z}}, 1);
        mask->SetPixel({{5, y, z}}, 1);
      }

    // A: two points, both outside the mask; B: the same line in steps of 0.5 mm; C: inside the
    // mask, across several voxel borders; D: diagonal in the x-y plane, both points outside.
    std::vector<std::array<double, 3>> fiberB;
    for (int j = 0; j <= 18; ++j)
      fiberB.push_back({0.5 * j, 5, 5});
    const std::vector<std::array<double, 3>> fiberC = {{3.6, 1, 1}, {4.2, 3.3, 2}, {5.4, 8, 8.7}};
    auto fib = MakeBundle({{{0, 5, 5}, {9, 5, 5}}, fiberB, fiberC, {{2, 2, 5}, {7, 7, 5}}});
    for (unsigned int i = 0; i < 4; ++i)
      fib->SetFiberWeight(i, 2.0f + i);

    auto inside = fib->RemoveFibersOutside(mask);
    CPPUNIT_ASSERT_EQUAL(4u, inside->GetNumFibers());
    AssertPiece(inside, 0, {3.5, 5, 5}, {5.5, 5, 5}, 2, 2.0f, "outside removed, A");
    AssertPiece(inside, 1, {3.5, 5, 5}, {5.5, 5, 5}, 5, 3.0f, "outside removed, B");
    CPPUNIT_ASSERT_MESSAGE("a fiber inside the mask keeps its points", FiberPoints(inside, 2) == FiberPoints(fib, 2));
    CPPUNIT_ASSERT_EQUAL(4.0f, inside->GetFiberWeight(2));
    AssertPiece(inside, 3, {3.5, 3.5, 5}, {5.5, 5.5, 5}, 2, 5.0f, "outside removed, D");

    auto outside = fib->RemoveFibersOutside(mask, true);
    CPPUNIT_ASSERT_EQUAL(6u, outside->GetNumFibers());
    AssertPiece(outside, 0, {0, 5, 5}, {3.5, 5, 5}, 2, 2.0f, "inside removed, A before the mask");
    AssertPiece(outside, 1, {5.5, 5, 5}, {9, 5, 5}, 2, 2.0f, "inside removed, A after the mask");
    AssertPiece(outside, 2, {0, 5, 5}, {3.5, 5, 5}, 8, 3.0f, "inside removed, B before the mask");
    AssertPiece(outside, 3, {5.5, 5, 5}, {9, 5, 5}, 8, 3.0f, "inside removed, B after the mask");
    AssertPiece(outside, 4, {2, 2, 5}, {3.5, 3.5, 5}, 2, 5.0f, "inside removed, D before the mask");
    AssertPiece(outside, 5, {5.5, 5.5, 5}, {7, 7, 5}, 2, 5.0f, "inside removed, D after the mask");
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

    auto fib = source->Clone();
    CPPUNIT_ASSERT(fib->RemoveShortFibers(medianLength));
    CPPUNIT_ASSERT(fib->GetNumFibers() < source->GetNumFibers());
    AssertWeightsFollowFibers(source, fib, "RemoveShortFibers");

    fib = source->Clone();
    CPPUNIT_ASSERT(fib->RemoveLongFibers(medianLength));
    CPPUNIT_ASSERT(fib->GetNumFibers() < source->GetNumFibers());
    AssertWeightsFollowFibers(source, fib, "RemoveLongFibers");

    fib = source->Clone();
    vnl_vector_fixed<double, 3> dir(0.0, 1.0, 0.0);
    fib->RemoveDir(dir, std::cos(5.0 * itk::Math::pi / 180));
    AssertWeightsFollowFibers(source, fib, "RemoveDir");

    fib = source->Clone();
    CPPUNIT_ASSERT(fib->ApplyCurvatureThreshold(5.0f, false));
    AssertWeightsFollowFibers(source, fib, "ApplyCurvatureThreshold");

    auto shortFibers = source->Clone();
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
    auto shortFibers = source->Clone();
    shortFibers->RemoveLongFibers(source->GetMedianFiberLength());

    mitk::OmpThreadGuard threads(1);
    auto compressedSingle = source->Clone();
    compressedSingle->Compress(0.1f);
    auto subtractedSingle = source->SubtractBundle(shortFibers);

    threads.SetNumberOfThreads(8);
    for (int run = 0; run < 3; ++run)
    {
      auto compressed = source->Clone();
      compressed->Compress(0.1f);
      AssertIdentical(compressedSingle, compressed, "Compress, run " + std::to_string(run));
      AssertIdentical(subtractedSingle, source->SubtractBundle(shortFibers), "SubtractBundle, run " + std::to_string(run));
    }
  }

  // #32: the thread count a failing case set stayed for all later cases of the driver.
  void OmpThreadGuard_FailingCase_RestoresThreadCount()
  {
    const int before = omp_get_max_threads();
    const int other = before == 3 ? 2 : 3;
    try
    {
      mitk::OmpThreadGuard threads(other);
      CPPUNIT_ASSERT_EQUAL(other, omp_get_max_threads());
      threads.SetNumberOfThreads(5);
      CPPUNIT_ASSERT_EQUAL(5, omp_get_max_threads());
      throw std::runtime_error("a failing assertion");
    }
    catch (const std::runtime_error &)
    {
    }
    CPPUNIT_ASSERT_EQUAL(before, omp_get_max_threads());
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
    auto shortFibers = source->Clone();
    shortFibers->RemoveLongFibers(source->GetMedianFiberLength());

    auto check = [](mitk::FiberBundle *fib, const std::string &operation)
    {
      CPPUNIT_ASSERT_MESSAGE(operation + " keeps the RAS flag", fib->IsRAS());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 77, static_cast<int>(fib->GetTrackVisHeader().dim[0]));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 1.25f, fib->GetTrackVisHeader().voxel_size[0]);
    };
    check(source->Clone(), "Clone");
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

  // #32: with as many new fibers as old ones, the new fibers kept the weights of the old ones.
  void SetFiberPolyData_EqualFiberCount_NewFibersGetWeightOne()
  {
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {1, 1, 0}}});
    fib->SetFiberWeight(0, 3.0f);
    fib->SetFiberWeight(1, 4.0f);
    auto other = MakeBundle({{{5, 5, 5}, {6, 5, 5}}, {{5, 6, 5}, {6, 6, 5}}});
    auto oldWeights = fib->GetFiberWeights();

    fib->SetFiberPolyData(other->GetFiberPolyData());
    CPPUNIT_ASSERT_EQUAL(2u, fib->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1.0f, fib->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(1.0f, fib->GetFiberWeight(1));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("a weight array taken before is not changed", 4.0f, oldWeights->GetValue(1));
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

  void ColorFibersByOrientation_OnePointFibers_GetColours()
  {
    // the point of the one-point fiber has no orientation
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}}, {{5, 5, 5}}, {{0, 3, 0}, {0, 4, 0}}});
    auto colours = fib->GetFiberColors();
    CPPUNIT_ASSERT_EQUAL(fib->GetFiberPolyData()->GetNumberOfPoints(), colours->GetNumberOfTuples());
    bool found = false;
    for (vtkIdType i = 0; i < fib->GetFiberPolyData()->GetNumberOfPoints(); ++i)
    {
      const double *p = fib->GetFiberPolyData()->GetPoint(i);
      if (p[0] == 5.0)
      {
        found = true;
        for (int c = 0; c < 4; ++c)
          CPPUNIT_ASSERT_EQUAL(255.0, colours->GetComponent(i, c));
      }
    }
    CPPUNIT_ASSERT(found);
    // a copy can take the colours over
    CPPUNIT_ASSERT(fib->Clone().IsNotNull());
  }

  void ColorFibersByCurvature_UsesFirstSegment()
  {
    // A right angle at the second point: its backward window holds the first segment, so its
    // curvature colour differs from that of a point on the straight part.
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 2, 0}, {1, 3, 0}, {1, 4, 0}, {1, 5, 0}, {1, 6, 0}}});
    fib->ColorFibersByCurvature(false, false, mitk::LookupTable::GRAYSCALE);
    auto colours = fib->GetFiberColors();
    vtkCell *cell = fib->GetFiberPolyData()->GetCell(0);
    const vtkIdType corner = cell->GetPointId(1);
    const vtkIdType straight = cell->GetPointId(6);
    CPPUNIT_ASSERT_MESSAGE("The corner point sees the bend",
                           colours->GetComponent(corner, 0) != colours->GetComponent(straight, 0));
  }

  void ResampleLinear_NoStep_LeavesFibersUnchanged()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    auto copy = fib->Clone();
    fib->ResampleLinear(0.0);
    CPPUNIT_ASSERT(copy->Equals(fib, 0.0));
  }

  void PrintSelf_EmptyBundle_Prints()
  {
    std::ostringstream out;
    mitk::FiberBundle::New()->Print(out);
    CPPUNIT_ASSERT(out.str().find("Number of fibers: 0") != std::string::npos);
    CPPUNIT_ASSERT(out.str().find("No fiber weights") != std::string::npos);
  }

  void PrintStatistics_WithoutObjectAddresses()
  {
    // MitkPrintFiberStatistics writes these statistics to a file, which must be the same on
    // every run: no object addresses as in the BaseData part of PrintSelf.
    auto fib = LoadWeightedFibers();
    std::ostringstream statistics;
    fib->PrintStatistics(statistics, 2);
    CPPUNIT_ASSERT(statistics.str().find("Number of fibers: 71") != std::string::npos);
    CPPUNIT_ASSERT(statistics.str().find("Fiber weight statistics") != std::string::npos);
    CPPUNIT_ASSERT(statistics.str().find("0x") == std::string::npos);

    std::ostringstream all;
    fib->Print(all);
    CPPUNIT_ASSERT(all.str().find("Geometry3D (0x") != std::string::npos);
  }

  void GeneratePolyDataByIds_InvalidId_Throws()
  {
    auto fib = LoadWeightedFibers();
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    auto polyData = fib->GeneratePolyDataByIds({0, 2}, weights);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(2), polyData->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(fib->GetFiberWeight(2), weights->GetValue(1));
    CPPUNIT_ASSERT_THROW(fib->GeneratePolyDataByIds({0, fib->GetNumFibers()}, weights), mitk::Exception);
  }

  void ColorFibersByFiberWeights_EqualWeightsAboveOne_OpaqueColours()
  {
    // equal weights force the colour range to [0, 1]; weight/max is then 3
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}}, {{0, 2, 0}, {1, 2, 0}}});
    fib->SetFiberWeights(3.0f);
    fib->ColorFibersByFiberWeights(true, mitk::LookupTable::JET);
    auto colours = fib->GetFiberColors();
    for (vtkIdType i = 0; i < colours->GetNumberOfTuples(); ++i)
      CPPUNIT_ASSERT_EQUAL(255.0, colours->GetComponent(i, 3));
  }

  void AddBundlesAndCompositeWithoutStorage_NoCrash()
  {
    auto fib = LoadWeightedFibers();
    auto other = LoadWeightedFibers();
    auto joined = fib->AddBundles({nullptr, other, nullptr});
    CPPUNIT_ASSERT_EQUAL(2 * fib->GetNumFibers(), joined->GetNumFibers());

    auto compositeNode = mitk::DataNode::New();
    compositeNode->SetData(mitk::PlanarFigureComposite::New());
    CPPUNIT_ASSERT_THROW(fib->ExtractFiberIdSubset(compositeNode, nullptr), mitk::Exception);
  }

  void Clone_CopiesFibersAndMetadataIndependently()
  {
    auto fib = LoadWeightedFibers();
    fib->SetSingleFiberColor(10, 20, 30, 1);
    fib->setIsRAS(true);
    auto header = fib->GetTrackVisHeader();
    header.dim[0] = 77;
    fib->SetTrackVisHeader(header);

    mitk::FiberBundle::Pointer clone = fib->Clone();
    CPPUNIT_ASSERT(clone.GetPointer() != fib.GetPointer());
    CPPUNIT_ASSERT(clone->GetFiberPolyData().GetPointer() != fib->GetFiberPolyData().GetPointer());
    CPPUNIT_ASSERT(fib->Equals(clone, 0.0));
    CPPUNIT_ASSERT_EQUAL(fib->GetNumberOfPoints(), clone->GetNumberOfPoints());
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberWeight(i), clone->GetFiberWeight(i));
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberLength(i), clone->GetFiberLength(i));
    }
    CPPUNIT_ASSERT_EQUAL(fib->GetFiberColors()->GetNumberOfTuples(), clone->GetFiberColors()->GetNumberOfTuples());
    for (vtkIdType i = 0; i < fib->GetFiberColors()->GetNumberOfValues(); ++i)
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberColors()->GetValue(i), clone->GetFiberColors()->GetValue(i));
    CPPUNIT_ASSERT_EQUAL(fib->GetMeanFiberLength(), clone->GetMeanFiberLength());
    CPPUNIT_ASSERT(clone->IsRAS());
    CPPUNIT_ASSERT_EQUAL(77, static_cast<int>(clone->GetTrackVisHeader().dim[0]));

    // Changing the arrays of the clone in place leaves the original as it was
    const double x = fib->GetFiberPolyData()->GetPoint(0)[0];
    const float weight = fib->GetFiberWeight(0);
    const unsigned char red = fib->GetFiberColors()->GetValue(0);
    const double moved[3] = {x + 5, 0, 0};
    clone->GetFiberPolyData()->GetPoints()->SetPoint(0, moved);
    clone->SetFiberWeight(0, weight + 1);
    clone->GetFiberColors()->SetValue(0, static_cast<unsigned char>(red + 1));
    CPPUNIT_ASSERT_EQUAL(x, fib->GetFiberPolyData()->GetPoint(0)[0]);
    CPPUNIT_ASSERT_EQUAL(weight, fib->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(red, fib->GetFiberColors()->GetValue(0));
    // The clone stays usable for operations that use the fiber ids
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(2), clone->GeneratePolyDataByIds({0, 1}, weights)->GetNumberOfLines());

    auto composite = mitk::PlanarFigureComposite::New();
    composite->setOperationType(mitk::PlanarFigureComposite::NOT);
    CPPUNIT_ASSERT_EQUAL(mitk::PlanarFigureComposite::NOT, composite->Clone()->getOperationType());
  }

  void FiberWeights_WrongSizeOrIndex_Throw()
  {
    auto fib = LoadWeightedFibers();
    const unsigned int numFibers = fib->GetNumFibers();
    const float lastWeight = fib->GetFiberWeight(numFibers - 1);

    auto tooShort = vtkSmartPointer<vtkFloatArray>::New();
    tooShort->SetNumberOfValues(numFibers - 1);
    tooShort->FillValue(2.0f);
    CPPUNIT_ASSERT_THROW(fib->SetFiberWeights(tooShort), mitk::Exception);
    CPPUNIT_ASSERT_THROW(fib->SetFiberWeights(vtkSmartPointer<vtkFloatArray>()), mitk::Exception);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("A rejected array leaves the weights", lastWeight, fib->GetFiberWeight(numFibers - 1));

    CPPUNIT_ASSERT_THROW(fib->GetFiberWeight(numFibers), mitk::Exception);
    CPPUNIT_ASSERT_THROW(fib->SetFiberWeight(numFibers, 1.0f), mitk::Exception);

    auto matching = vtkSmartPointer<vtkFloatArray>::New();
    matching->SetNumberOfValues(numFibers);
    matching->FillValue(2.0f);
    fib->SetFiberWeights(matching);
    CPPUNIT_ASSERT_EQUAL(2.0f, fib->GetFiberWeight(numFibers - 1));
  }

  // Safety net for the restructuring of FiberBundle: the results the operations give today.

  void LengthStatistics_OriginalFib_AsRecorded()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    CPPUNIT_ASSERT_EQUAL(71u, fib->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(19066u, fib->GetNumberOfPoints());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20.4, fib->GetMinFiberLength(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(67.8, fib->GetMaxFiberLength(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(53.50709, fib->GetMeanFiberLength(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(62.80009, fib->GetMedianFiberLength(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(16.07929, fib->GetLengthStDev(), 1e-3);
    double sum = 0;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      sum += fib->GetFiberLength(i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fib->GetMeanFiberLength() * fib->GetNumFibers(), sum, 1e-2);
  }

  void ResampleToNumPoints_RealFibers_GiveNPointsAndKeepEndpoints()
  {
    auto source = LoadWeightedFibers();
    auto fib = source->Clone();
    fib->ResampleToNumPoints(12);
    CPPUNIT_ASSERT_EQUAL(source->GetNumFibers(), fib->GetNumFibers());
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      vtkCell *sourceCell = source->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(12), cell->GetNumberOfPoints());
      for (int c = 0; c < 3; ++c)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceCell->GetPoints()->GetPoint(0)[c], cell->GetPoints()->GetPoint(0)[c], 1e-4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceCell->GetPoints()->GetPoint(sourceCell->GetNumberOfPoints() - 1)[c],
                                     cell->GetPoints()->GetPoint(11)[c], 1e-3);
      }
      CPPUNIT_ASSERT_EQUAL(source->GetFiberWeight(i), fib->GetFiberWeight(i));
    }
  }

  void FilterByWeights_BothDirections_SplitAtThreshold()
  {
    auto source = LoadWeightedFibers();
    auto heavy = source->FilterByWeights(1.0f);
    auto light = source->FilterByWeights(1.0f, true);
    CPPUNIT_ASSERT(heavy->GetNumFibers() > 0 && light->GetNumFibers() > 0);
    CPPUNIT_ASSERT_EQUAL(source->GetNumFibers(), heavy->GetNumFibers() + light->GetNumFibers());
    for (unsigned int i = 0; i < heavy->GetNumFibers(); ++i)
      CPPUNIT_ASSERT(heavy->GetFiberWeight(i) > 1.0f);
    for (unsigned int i = 0; i < light->GetNumFibers(); ++i)
      CPPUNIT_ASSERT(light->GetFiberWeight(i) <= 1.0f);
    AssertWeightsFollowFibers(source, heavy, "FilterByWeights");
    AssertWeightsFollowFibers(source, light, "FilterByWeights inverted");
  }

  void AddBundles_EqualsRepeatedAddBundle()
  {
    auto a = LoadWeightedFibers();
    auto b = LoadFibers("fiberBundleX.fib");
    auto c = LoadWeightedFibers();
    c->TranslateFibers(3, 0, 0);
    auto joined = a->AddBundles({b, c});
    auto stepwise = a->AddBundle(b)->AddBundle(c);
    CPPUNIT_ASSERT_EQUAL(a->GetNumFibers() + b->GetNumFibers() + c->GetNumFibers(), joined->GetNumFibers());
    AssertIdentical(stepwise, joined, "AddBundles");
    // the first fibers are those of a, in order
    for (unsigned int i = 0; i < a->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(a->GetFiberWeight(i), joined->GetFiberWeight(i));
  }

  void RotateAroundAxis_EqualsTransformFibersWithoutTranslation()
  {
    auto rotated = LoadWeightedFibers();
    auto transformed = rotated->Clone();
    const auto colours = rotated->GetFiberColors();
    std::vector<unsigned char> before(colours->GetPointer(0), colours->GetPointer(0) + colours->GetNumberOfValues());
    rotated->RotateAroundAxis(10, 20, 30);
    transformed->TransformFibers(10, 20, 30, 0, 0, 0);
    CPPUNIT_ASSERT(rotated->Equals(transformed, 1e-4));
    for (unsigned int i = 0; i < rotated->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(transformed->GetFiberWeight(i), rotated->GetFiberWeight(i));
    // RotateAroundAxis keeps the colours of the points
    const auto after = rotated->GetFiberColors();
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(before.size()), after->GetNumberOfValues());
    for (vtkIdType i = 0; i < after->GetNumberOfValues(); ++i)
      CPPUNIT_ASSERT_EQUAL(before[i], after->GetValue(i));
  }

  void Clone_KeepsDataAndIsIndependent()
  {
    auto fib = LoadWeightedFibers();
    fib->ColorFibersByLength(false, false, mitk::LookupTable::JET);
    fib->setIsRAS(true);
    auto header = fib->GetTrackVisHeader();
    header.dim[1] = 55;
    fib->SetTrackVisHeader(header);

    auto copy = fib->Clone();
    CPPUNIT_ASSERT(copy->GetFiberPolyData().GetPointer() != fib->GetFiberPolyData().GetPointer());
    AssertIdentical(fib, copy, "Clone");
    CPPUNIT_ASSERT_EQUAL(fib->GetFiberColors()->GetNumberOfValues(), copy->GetFiberColors()->GetNumberOfValues());
    for (vtkIdType i = 0; i < fib->GetFiberColors()->GetNumberOfValues(); ++i)
      CPPUNIT_ASSERT_EQUAL(fib->GetFiberColors()->GetValue(i), copy->GetFiberColors()->GetValue(i));
    CPPUNIT_ASSERT(copy->IsRAS());
    CPPUNIT_ASSERT_EQUAL(55, static_cast<int>(copy->GetTrackVisHeader().dim[1]));

    copy->TranslateFibers(1, 0, 0);
    copy->SetFiberWeights(7.0f);
    CPPUNIT_ASSERT(!fib->Equals(copy, 0.5));
    CPPUNIT_ASSERT(fib->GetFiberWeight(0) != 7.0f);
  }

  void EmptyBundle_EveryOperation_StaysEmpty()
  {
    auto empty = []() { return mitk::FiberBundle::New(); };
    auto check = [](mitk::FiberBundle *fib, const std::string &operation)
    {
      CPPUNIT_ASSERT_MESSAGE(operation, fib != nullptr);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 0u, fib->GetNumFibers());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation, 0u, fib->GetNumberOfPoints());
    };
    auto fib = empty();
    fib->ResampleSpline(1.0f);
    check(fib, "ResampleSpline");
    fib->ResampleLinear(1.0);
    check(fib, "ResampleLinear");
    fib->ResampleToNumPoints(5);
    check(fib, "ResampleToNumPoints");
    fib->Compress(0.1f);
    check(fib, "Compress");
    // "no fiber would be left": every length is above the maximum 0
    CPPUNIT_ASSERT(!fib->RemoveShortFibers(10));
    CPPUNIT_ASSERT(fib->RemoveLongFibers(10));
    CPPUNIT_ASSERT(!fib->ApplyCurvatureThreshold(5, false));
    vnl_vector_fixed<double, 3> dir(0.0, 1.0, 0.0);
    fib->RemoveDir(dir, 0.5);
    fib->MirrorFibers(1);
    fib->RotateAroundAxis(10, 0, 0);
    fib->TranslateFibers(1, 2, 3);
    fib->ScaleFibers(2, 2, 2);
    fib->TransformFibers(1, 2, 3, 4, 5, 6);
    fib->TransformFibers(itk::ScalableAffineTransform<mitk::ScalarType>::New());
    check(fib, "transforms");
    check(fib->FilterByWeights(0.5f), "FilterByWeights");
    check(fib->SubsampleFibers(0.5f, false), "SubsampleFibers");
    check(fib->AddBundle(empty()), "AddBundle");
    check(fib->AddBundles({empty(), empty()}), "AddBundles");
    check(fib->SubtractBundle(empty()), "SubtractBundle");
    check(fib->Clone(), "Clone");
    fib->ColorFibersByLength(true, true, mitk::LookupTable::JET);
    fib->ColorFibersByCurvature(true, true, mitk::LookupTable::JET);
    fib->ColorFibersByFiberWeights(true, mitk::LookupTable::JET);
    fib->ColorFibersByOrientation();
    fib->SetFiberColors(1, 2, 3);
    check(fib, "colouring");
    CPPUNIT_ASSERT(fib->Equals(empty()));
    std::ostringstream out;
    fib->Print(out);
    CPPUNIT_ASSERT(out.str().find("Number of fibers: 0") != std::string::npos);
  }

  // Method checks: the geometric operations against their definitions on synthetic fibers.

  static double Distance(const double *a, const double *b)
  {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
  }

  /** Distance of p from the polyline through the given points */
  static double DistanceToPolyline(const double *p, const std::vector<std::array<double, 3>> &line)
  {
    double best = std::numeric_limits<double>::max();
    for (std::size_t k = 0; k + 1 < line.size(); ++k)
    {
      double d[3], w[3];
      double dd = 0, t = 0;
      for (int c = 0; c < 3; ++c)
      {
        d[c] = line[k + 1][c] - line[k][c];
        w[c] = p[c] - line[k][c];
        dd += d[c] * d[c];
        t += w[c] * d[c];
      }
      t = dd > 0 ? std::min(1.0, std::max(0.0, t / dd)) : 0.0;
      double q[3];
      for (int c = 0; c < 3; ++c)
        q[c] = line[k][c] + t * d[c];
      best = std::min(best, Distance(p, q));
    }
    return best;
  }

  void ResampleLinear_SyntheticFibers_EqualChordSpacingOnTheFiber()
  {
    // Linear resampling places the points on the fiber, each at the given straight-line distance
    // from the previous one; the end point is kept.
    const std::vector<std::vector<std::array<double, 3>>> fibers = {
      {{0, 0, 0}, {0.3, 0, 0}, {2.5, 0, 0}, {4, 0, 0}, {10, 0, 0}},
      {{0, 5, 0}, {4, 5, 0}, {4, 9, 0}, {4, 9, 3.3}}};
    auto fib = MakeBundle(fibers);
    const double step = 0.7;
    fib->ResampleLinear(step);
    CPPUNIT_ASSERT_EQUAL(2u, fib->GetNumFibers());
    for (unsigned int i = 0; i < 2; ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      vtkPoints *points = cell->GetPoints();
      const auto n = cell->GetNumberOfPoints();
      CPPUNIT_ASSERT(n > 10);
      for (vtkIdType j = 0; j < n; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, DistanceToPolyline(points->GetPoint(j), fibers[i]), 1e-5);
      for (vtkIdType j = 1; j < n - 1; ++j)
      {
        double a[3], b[3];
        points->GetPoint(j - 1, a);
        points->GetPoint(j, b);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(step, Distance(a, b), 1e-5);
      }
      double a[3], b[3];
      points->GetPoint(n - 2, a);
      points->GetPoint(n - 1, b);
      CPPUNIT_ASSERT(Distance(a, b) <= step + 1e-5);
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(fibers[i].back()[c], b[c], 1e-6);
    }
    // the straight fiber: 10 mm in steps of 0.7 mm
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(16), fib->GetFiberPolyData()->GetCell(0)->GetNumberOfPoints());
  }

  void ResampleSpline_StraightFiber_PointsOnTheLine()
  {
    // The spline through collinear points is the line; ceil(length/distance) segments.
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}, {2.5, 0, 0}, {6, 0, 0}, {10, 0, 0}}});
    fib->ResampleSpline(1.0f);
    vtkCell *cell = fib->GetFiberPolyData()->GetCell(0);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(11), cell->GetNumberOfPoints());
    for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
    {
      const double *p = cell->GetPoints()->GetPoint(j);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, p[1], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, p[2], 1e-6);
      CPPUNIT_ASSERT(p[0] >= -1e-6 && p[0] <= 10 + 1e-6);
      if (j > 0)
        CPPUNIT_ASSERT_MESSAGE("monotonic along the line", p[0] > cell->GetPoints()->GetPoint(j - 1)[0]);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, cell->GetPoints()->GetPoint(0)[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, cell->GetPoints()->GetPoint(10)[0], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, fib->GetFiberLength(0), 1e-4);
  }

  void ResampleToNumPoints_StraightFiber_EqualSpacing()
  {
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}, {2.5, 0, 0}, {6, 0, 0}, {10, 0, 0}}});
    fib->ResampleToNumPoints(6);
    vtkCell *cell = fib->GetFiberPolyData()->GetCell(0);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(6), cell->GetNumberOfPoints());
    for (vtkIdType j = 0; j < 6; ++j)
    {
      const double *p = cell->GetPoints()->GetPoint(j);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0 * j, p[0], 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, p[1], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, p[2], 1e-6);
    }
  }

  /** Arc length along the polyline at which p lies on it, searched from segment *segment on
   *  (points come in order along the line); -1 if p is not on the polyline. */
  static double ArcLengthOnPolyline(const double *p, const std::vector<std::array<double, 3>> &line, std::size_t *segment)
  {
    double arc = 0;
    for (std::size_t k = 0; k < *segment; ++k)
      arc += Distance(line[k].data(), line[k + 1].data());
    for (std::size_t k = *segment; k + 1 < line.size(); ++k)
    {
      const double length = Distance(line[k].data(), line[k + 1].data());
      if (DistanceToPolyline(p, {line[k], line[k + 1]}) < 1e-4)
      {
        *segment = k;
        return arc + Distance(line[k].data(), p);
      }
      arc += length;
    }
    return -1;
  }

  // #16: the points were placed at equal straight-line (chord) distance, which on curved fibers
  // is not equal arc length, and fibers that already had the requested number of points were
  // kept with their uneven spacing.
  void ResampleToNumPoints_CurvedFibers_EqualArcLengthSpacing()
  {
    // A half circle of radius 10 mm and two turns of a helix, both with uneven point steps, and a
    // straight fiber that already has 7 unevenly spaced points.
    std::vector<std::array<double, 3>> arc, helix;
    for (int j = 0; j <= 60; ++j)
    {
      const double u = j / 60.0;
      const double phi = itk::Math::pi * u * u;
      arc.push_back({10 * std::cos(phi), 10 * std::sin(phi), 0});
    }
    for (int j = 0; j <= 200; ++j)
    {
      const double u = j / 200.0;
      const double phi = 4 * itk::Math::pi * std::sqrt(u);
      helix.push_back({5 * std::cos(phi), 5 * std::sin(phi), 20 * std::sqrt(u)});
    }
    const std::vector<std::array<double, 3>> line = {{0, 0, 0}, {0.5, 0, 0}, {1, 0, 0}, {5, 0, 0}, {6, 0, 0}, {11, 0, 0}, {12, 0, 0}};
    const std::vector<std::vector<std::array<double, 3>>> fibers = {arc, helix, line};
    auto fib = MakeBundle(fibers);
    for (unsigned int i = 0; i < 3; ++i)
      fib->SetFiberWeight(i, 1.5f + i);
    // the input as stored (float), against which the result is measured
    std::vector<std::vector<std::array<double, 3>>> stored;
    std::vector<double> lengths;
    for (unsigned int i = 0; i < 3; ++i)
    {
      stored.push_back(FiberPoints(fib, i));
      lengths.push_back(fib->GetFiberLength(i));
    }

    const unsigned int n = 7;
    fib->ResampleToNumPoints(n);
    CPPUNIT_ASSERT_EQUAL(3u, fib->GetNumFibers());
    for (unsigned int i = 0; i < 3; ++i)
    {
      const std::string message = "fiber " + std::to_string(i);
      const auto points = FiberPoints(fib, i);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(message, static_cast<std::size_t>(n), points.size());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": weight", 1.5f + i, fib->GetFiberWeight(i));
      for (int c = 0; c < 3; ++c)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message + ": first point", stored[i].front()[c], points.front()[c], 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message + ": last point", stored[i].back()[c], points.back()[c], 1e-5);
      }
      std::size_t segment = 0;
      for (unsigned int k = 0; k < n; ++k)
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message + ": arc length of point " + std::to_string(k),
                                             lengths[i] * k / (n - 1),
                                             ArcLengthOnPolyline(points[k].data(), stored[i], &segment), 1e-3);
    }
  }

  void Compress_SyntheticFibers_RemovesOnlyPointsWithinTheError()
  {
    // A point is removed while its distance from the line through its kept neighbours is below
    // the error: collinear points and a 0.05 mm wobble go, corners and a 1 mm zigzag stay.
    auto fib = MakeBundle({{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0}, {5, 0, 0}},
                           {{0, 5, 0}, {1, 5, 0}, {2, 5, 0}, {2, 6, 0}, {2, 7, 0}, {2, 8, 0}},
                           {{0, 10, 0}, {1, 11, 0}, {2, 10, 0}, {3, 11, 0}, {4, 10, 0}},
                           {{0, 15, 0}, {1, 15.05, 0}, {2, 15, 0}, {3, 14.95, 0}, {4, 15, 0}}});
    fib->Compress(0.1f);
    CPPUNIT_ASSERT_EQUAL(4u, fib->GetNumFibers());
    std::map<long long, vtkIdType> pointsByFirstY;
    for (unsigned int i = 0; i < 4; ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      pointsByFirstY[std::llround(cell->GetPoints()->GetPoint(0)[1])] = cell->GetNumberOfPoints();
    }
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(2), pointsByFirstY[0]);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(3), pointsByFirstY[5]);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(5), pointsByFirstY[10]);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(2), pointsByFirstY[15]);
  }

  void ApplyCurvatureThreshold_CircularArc_RadiusDecides()
  {
    // Three consecutive points of an arc of radius 10 have the circumradius 10.
    std::vector<std::array<double, 3>> arc;
    for (int a = 0; a <= 90; a += 10)
      arc.push_back({10 * std::cos(a * itk::Math::pi / 180), 10 * std::sin(a * itk::Math::pi / 180), 0});
    const std::vector<std::array<double, 3>> line = {{0, 20, 0}, {1, 20, 0}, {2, 20, 0}, {3, 20, 0}};

    auto fib = MakeBundle({arc, line});
    CPPUNIT_ASSERT(fib->ApplyCurvatureThreshold(9.5f, true));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("radius 10 is above 9.5 mm", 2u, fib->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(arc.size()), fib->GetFiberPolyData()->GetCell(0)->GetNumberOfPoints());

    fib = MakeBundle({arc, line});
    CPPUNIT_ASSERT(fib->ApplyCurvatureThreshold(10.5f, true));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("radius 10 is below 10.5 mm: the arc goes", 1u, fib->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0, fib->GetFiberPolyData()->GetCell(0)->GetPoints()->GetPoint(0)[1], 1e-6);
  }

  void ColorFibersByOrientation_ThreePointFiber_SegmentDirections()
  {
    // End points take the direction of their segment, inner points the mean of both segments.
    auto fib = MakeBundle({{{0, 0, 0}, {2, 0, 0}, {2, 2, 0}}});
    auto colours = fib->GetFiberColors();
    vtkCell *cell = fib->GetFiberPolyData()->GetCell(0);
    const std::array<std::array<double, 4>, 3> expected = {{{255, 0, 0, 255},
                                                             {static_cast<double>(static_cast<unsigned char>(255.0 * std::sqrt(0.5))),
                                                              static_cast<double>(static_cast<unsigned char>(255.0 * std::sqrt(0.5))), 0, 255},
                                                             {0, 255, 0, 255}}};
    for (vtkIdType j = 0; j < 3; ++j)
      for (int c = 0; c < 4; ++c)
        CPPUNIT_ASSERT_EQUAL(expected[j][c], colours->GetComponent(cell->GetPointId(j), c));
  }

  void ColorFibersByScalarMap_Mask_ColoursAndWeightsFromTheMask()
  {
    auto fib = LoadFibers("FiberProcessing/original.fib");
    auto mask = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("FiberProcessing/MASK.nrrd"));
    fib->ColorFibersByScalarMap(mask, true, true, mitk::LookupTable::JET, 1.0, false);
    auto colours = fib->GetFiberColors();
    CPPUNIT_ASSERT_EQUAL(fib->GetFiberPolyData()->GetNumberOfPoints(), colours->GetNumberOfTuples());
    // without interpolation the opacity is 255 * mask value, i.e. 0 or 255
    std::set<double> alphas;
    for (vtkIdType i = 0; i < colours->GetNumberOfTuples(); ++i)
      alphas.insert(colours->GetComponent(i, 3));
    CPPUNIT_ASSERT(alphas == std::set<double>({0.0, 255.0}));
    // the fiber weight is the mean mask value along the fiber
    bool inside = false;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      CPPUNIT_ASSERT(fib->GetFiberWeight(i) >= 0.0f && fib->GetFiberWeight(i) <= 1.0f);
      inside = inside || fib->GetFiberWeight(i) > 0.0f;
    }
    CPPUNIT_ASSERT(inside);
  }

  /** A 4x1x1 image with the values 0, 1, 1, 0 along x (voxel i at x = i) */
  template <typename TPixel>
  static mitk::Image::Pointer EdgeImage()
  {
    using ImageType = itk::Image<TPixel, 3>;
    typename ImageType::RegionType region;
    region.SetSize({4, 1, 1});
    auto image = ImageType::New();
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0);
    image->SetPixel({1, 0, 0}, 1);
    image->SetPixel({2, 0, 0}, 1);
    return mitk::GrabItkImageMemory(image.GetPointer());
  }

  // C-L16: interpolated values of signed integer images were truncated to int, so a point
  // halfway across a mask edge counted as 0.
  void ColorFibersByScalarMap_SignedIntegerImage_InterpolatedLikeFloat()
  {
    auto fromShort = MakeBundle({{{0.5, 0, 0}, {1.0, 0, 0}, {1.5, 0, 0}}});
    auto fromFloat = fromShort->Clone();
    fromShort->ColorFibersByScalarMap(EdgeImage<short>(), true, true, mitk::LookupTable::JET, 1.0, true);
    fromFloat->ColorFibersByScalarMap(EdgeImage<float>(), true, true, mitk::LookupTable::JET, 1.0, true);

    CPPUNIT_ASSERT_DOUBLES_EQUAL((0.5 + 1.0 + 1.0) / 3, fromFloat->GetFiberWeight(0), 1e-6);
    CPPUNIT_ASSERT_EQUAL(fromFloat->GetFiberWeight(0), fromShort->GetFiberWeight(0));
    auto expected = fromFloat->GetFiberColors();
    auto actual = fromShort->GetFiberColors();
    CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfValues(), actual->GetNumberOfValues());
    for (vtkIdType i = 0; i < expected->GetNumberOfValues(); ++i)
      CPPUNIT_ASSERT_EQUAL(static_cast<int>(expected->GetValue(i)), static_cast<int>(actual->GetValue(i)));
    // the half-way point is half opaque
    CPPUNIT_ASSERT_EQUAL(127, static_cast<int>(actual->GetComponent(0, 3)));
  }

  void TransformFibers_ItkTransform_MovesEveryPoint()
  {
    auto source = LoadWeightedFibers();
    auto fib = source->Clone();
    auto transform = itk::ScalableAffineTransform<mitk::ScalarType>::New();
    itk::ScalableAffineTransform<mitk::ScalarType>::OutputVectorType axis;
    axis[0] = 0.2;
    axis[1] = 0.5;
    axis[2] = 0.8;
    transform->Rotate3D(axis, 0.3);
    transform->Scale(1.1);
    auto translation = transform->GetTranslation();
    translation[0] = 1;
    translation[1] = -2;
    translation[2] = 3;
    transform->SetTranslation(translation);
    fib->TransformFibers(transform);
    CPPUNIT_ASSERT_EQUAL(source->GetNumFibers(), fib->GetNumFibers());
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
    {
      vtkCell *cell = fib->GetFiberPolyData()->GetCell(i);
      vtkCell *sourceCell = source->GetFiberPolyData()->GetCell(i);
      CPPUNIT_ASSERT_EQUAL(sourceCell->GetNumberOfPoints(), cell->GetNumberOfPoints());
      for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j)
      {
        auto expected = transform->TransformPoint(itk::Point<double, 3>(sourceCell->GetPoints()->GetPoint(j)));
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[c], cell->GetPoints()->GetPoint(j)[c], 1e-3);
      }
      CPPUNIT_ASSERT_EQUAL(source->GetFiberWeight(i), fib->GetFiberWeight(i));
    }
  }

  /** Every point lies in the volume the TrackVis header declares: its corner-based voxel
   *  coordinate (vox_to_ras^-1 * RAS + 0.5) is within [0, dim], as dipy's bounding box check asks. */
  static void AssertHeaderContainsAllPoints(mitk::FiberBundle *fib, const std::string &message)
  {
    const auto header = fib->GetTrackVisHeader();
    vnl_matrix_fixed<double, 3, 3> a;
    vnl_vector_fixed<double, 3> t;
    for (int r = 0; r < 3; ++r)
    {
      t[r] = header.vox_to_ras[r][3];
      for (int c = 0; c < 3; ++c)
        a[r][c] = header.vox_to_ras[r][c];
    }
    const auto inverse = vnl_inverse(a);
    auto points = fib->GetFiberPolyData()->GetPoints();
    CPPUNIT_ASSERT_MESSAGE(message, points->GetNumberOfPoints() > 0);
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    {
      const double *p = points->GetPoint(i);
      const vnl_vector_fixed<double, 3> ras(-p[0], -p[1], p[2]);
      const auto index = inverse * (ras - t);
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_MESSAGE(message + ": point " + std::to_string(i),
                               index[c] + 0.5 >= -1e-3 && index[c] + 0.5 <= header.dim[c] + 1e-3);
    }
  }

  void TrackVisHeader_JoinedAndMovedBundles_ContainEveryPoint()
  {
    auto a = LoadFibers("FiberProcessing/original.fib");
    auto b = LoadFibers("fiberBundleX.fib");
    CPPUNIT_ASSERT_MESSAGE("A .fib file sets no header", !a->HasTrackVisHeader());

    // The joins of FiberJoin, FitFibersToImage and FiberClustering, and joins in the views
    AssertHeaderContainsAllPoints(mitk::FiberBundle::New()->AddBundles({a, b}), "New()->AddBundles");
    AssertHeaderContainsAllPoints(a->AddBundle(b), "AddBundle");
    AssertHeaderContainsAllPoints(b->AddBundles({a}), "AddBundles");

    auto moved = a->Clone();
    moved->TranslateFibers(100, 0, 0);
    AssertHeaderContainsAllPoints(moved, "TranslateFibers");
    moved->RotateAroundAxis(0, 0, 90);
    AssertHeaderContainsAllPoints(moved, "RotateAroundAxis");

    // The first bundle with fibers gives the join its coordinate convention and a set header
    auto withReference = a->Clone();
    withReference->setIsRAS(true);
    auto header = withReference->GetTrackVisHeader();
    header.dim[0] = 77;
    withReference->SetTrackVisHeader(header);
    CPPUNIT_ASSERT(withReference->HasTrackVisHeader());
    auto joined = mitk::FiberBundle::New()->AddBundles({nullptr, withReference, b});
    CPPUNIT_ASSERT_MESSAGE("RAS flag of the first bundle with fibers", joined->IsRAS());
    CPPUNIT_ASSERT(joined->HasTrackVisHeader());
    CPPUNIT_ASSERT_EQUAL(77, static_cast<int>(joined->GetTrackVisHeader().dim[0]));
    CPPUNIT_ASSERT_MESSAGE("A derived header is not copied", !a->AddBundle(withReference)->HasTrackVisHeader());
  }

  /** Every fiber of result has the weight that the fiber starting at the same point has in one of
   *  the sources. */
  static void AssertWeightsByFirstPoint(const std::vector<mitk::FiberBundle *> &sources, mitk::FiberBundle *result, const std::string &operation)
  {
    std::map<PointKey, float> weights;
    for (auto source : sources)
      for (unsigned int i = 0; i < source->GetNumFibers(); ++i)
        weights[Key(source->GetFiberPolyData()->GetCell(i)->GetPoints()->GetPoint(0))] = source->GetFiberWeight(i);
    CPPUNIT_ASSERT_MESSAGE(operation, result->GetNumFibers() > 0);
    for (unsigned int i = 0; i < result->GetNumFibers(); ++i)
    {
      auto it = weights.find(Key(result->GetFiberPolyData()->GetCell(i)->GetPoints()->GetPoint(0)));
      CPPUNIT_ASSERT_MESSAGE(operation + ": fiber comes from a source", it != weights.end());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(operation + ": weight of fiber " + std::to_string(i), it->second, result->GetFiberWeight(i));
    }
  }

  void Weights_OnePointFibers_KeepPlaceAndWeight()
  {
    // Fibers of one point keep their place among the others (#20) and their weights.
    auto a = MakeBundle({{{0, 0, 0}, {0, 1, 0}}, {{1, 0, 0}, {1, 1, 0}}});
    a->SetFiberWeight(0, 1.0f);
    a->SetFiberWeight(1, 2.0f);
    auto b = MakeBundle({{{9, 9, 9}}, {{5, 0, 0}, {5, 1, 0}}});
    b->SetFiberWeight(0, 3.0f);
    b->SetFiberWeight(1, 4.0f);

    auto joined = a->AddBundle(b);
    CPPUNIT_ASSERT_EQUAL(4u, joined->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("the one-point fiber of b follows the fibers of a", static_cast<vtkIdType>(1),
                                 joined->GetFiberPolyData()->GetCell(2)->GetNumberOfPoints());
    AssertWeightsByFirstPoint({a, b}, joined, "AddBundle");
    AssertWeightsByFirstPoint({a, b}, a->AddBundles({b}), "AddBundles");

    for (unsigned int i = 0; i < joined->GetNumFibers(); ++i)
      joined->SetFiberWeight(i, 10.0f + i);
    AssertWeightsByFirstPoint({joined}, joined->SubsampleFibers(1.0f, false), "SubsampleFibers");
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    auto polyData = joined->GeneratePolyDataByIds({3, 2, 1, 0}, weights);
    auto subset = mitk::FiberBundle::New(polyData);
    subset->SetFiberWeights(weights);
    AssertWeightsByFirstPoint({joined}, subset, "GeneratePolyDataByIds");

    // A two-point fiber without length becomes a fiber of one point
    auto resampled = MakeBundle({{{0, 0, 0}, {0, 3, 0}}, {{7, 7, 7}, {7, 7, 7}}});
    resampled->SetFiberWeight(0, 5.0f);
    resampled->SetFiberWeight(1, 6.0f);
    auto source = resampled->Clone();
    resampled->ResampleLinear(1);
    AssertWeightsByFirstPoint({source}, resampled, "ResampleLinear");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleOperations)
