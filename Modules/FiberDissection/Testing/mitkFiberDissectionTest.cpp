/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkSphereInteractor.h>
#include <mitkStreamlineFeatureExtractor.h>
#include <mitkStreamlineInteractorBrush.h>

#include <mitkDataNode.h>
#include <mitkException.h>
#include <mitkFiberBundle.h>
#include <mitkRenderingManagerFactory.h>
#include <mitkSurface.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <set>
#include <vector>

namespace
{
  using Point = std::array<double, 3>;
  using Fiber = std::vector<Point>;

  /** A straight fiber from a to b with the point count the classifier needs. */
  Fiber Line(const Point &a, const Point &b)
  {
    const unsigned int n = mitk::StreamlineFeatureExtractor::NumPoints;
    Fiber fiber;
    for (unsigned int i = 0; i < n; ++i)
    {
      const double t = static_cast<double>(i) / (n - 1);
      fiber.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])});
    }
    return fiber;
  }

  mitk::FiberBundle::Pointer Bundle(const std::vector<Fiber> &fibers, const std::vector<float> &weights = {})
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
    auto bundle = mitk::FiberBundle::New(polyData);
    if (!weights.empty())
    {
      auto weightArray = vtkSmartPointer<vtkFloatArray>::New();
      weightArray->SetNumberOfValues(weights.size());
      for (std::size_t i = 0; i < weights.size(); ++i)
        weightArray->SetValue(i, weights[i]);
      bundle->SetFiberWeights(weightArray);
    }
    return bundle;
  }

  Point FirstPoint(const mitk::FiberBundle *bundle, unsigned int fiber)
  {
    Point p;
    bundle->GetFiberPolyData()->GetCell(fiber)->GetPoints()->GetPoint(0, p.data());
    return p;
  }

  std::vector<unsigned int> AllIds(const mitk::FiberBundle *bundle)
  {
    std::vector<unsigned int> ids(bundle->GetNumFibers());
    std::iota(ids.begin(), ids.end(), 0);
    return ids;
  }

  /** Three well separated clusters of 60 fibers in random orientation: 0-59 along x (the target
   *  tract), 60-119 along y, 120-179 along z. */
  std::vector<Fiber> ThreeClusters()
  {
    std::mt19937 rng(1);
    std::normal_distribution<double> noise(0, 0.5);
    std::vector<Fiber> fibers;
    for (int i = 0; i < 60; ++i)
    {
      const Point s{0, noise(rng), noise(rng)}, e{40, noise(rng), noise(rng)};
      fibers.push_back(i % 2 ? Line(s, e) : Line(e, s));
    }
    for (int i = 0; i < 60; ++i)
    {
      const Point s{20 + noise(rng), -20, noise(rng)}, e{20 + noise(rng), 20, noise(rng)};
      fibers.push_back(i % 2 ? Line(s, e) : Line(e, s));
    }
    for (int i = 0; i < 60; ++i)
    {
      const Point s{noise(rng) - 30, noise(rng), -20}, e{noise(rng) - 30, noise(rng), 20};
      fibers.push_back(i % 2 ? Line(s, e) : Line(e, s));
    }
    return fibers;
  }

  /** Every sixth fiber as prototype, instead of the random subset farthest first selection. */
  mitk::FiberBundle::Pointer Prototypes(mitk::FiberBundle *pool)
  {
    std::vector<unsigned int> ids;
    for (unsigned int i = 0; i < pool->GetNumFibers(); i += 6)
      ids.push_back(i);
    return mitk::StreamlineFeatureExtractor::SubsetByIds(pool, ids);
  }

  /** A sphere ROI node as the SphereInteractor creates it: a sphere surface around the origin of
   *  its geometry, and its radius as property. */
  mitk::DataNode::Pointer SphereRoi(const Point &center, float radius)
  {
    // the interactor requests render updates; registered once for the whole test driver
    static mitk::TestingRenderingManagerFactory renderingManagerFactory;

    auto node = mitk::DataNode::New();
    node->SetData(mitk::Surface::New());
    mitk::Point3D origin;
    origin[0] = center[0];
    origin[1] = center[1];
    origin[2] = center[2];
    node->GetData()->GetGeometry()->SetOrigin(origin);
    node->SetFloatProperty(mitk::SphereInteractor::DATANODE_PROPERTY_SIZE, radius);
    mitk::SphereInteractor::UpdateSurface(node);
    return node;
  }

  /** The fibers of bundle that the sphere ROI selection of the Interactive Fiber Dissection view
   *  keeps, as their first points. */
  std::vector<Point> SelectBySphereRois(mitk::FiberBundle *bundle, const mitk::DataNode::Pointer &start, const mitk::DataNode::Pointer &end)
  {
    auto workingNode = mitk::DataNode::New();
    workingNode->SetData(bundle);
    auto reducedNode = mitk::DataNode::New();
    reducedNode->SetData(mitk::FiberBundle::New());

    auto interactor = mitk::SphereInteractor::New();
    interactor->workingBundleNode(workingNode, reducedNode);
    interactor->m_startDataNode = start;
    interactor->m_endDataNode = end;
    interactor->ExtractFibers();

    auto reduced = dynamic_cast<mitk::FiberBundle *>(reducedNode->GetData());
    std::vector<Point> firstPoints;
    for (unsigned int i = 0; i < reduced->GetNumFibers(); ++i)
      firstPoints.push_back(FirstPoint(reduced, i));
    return firstPoints;
  }
}

class mitkFiberDissectionTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberDissectionTestSuite);
  MITK_TEST(MoveFiber_KeepsFiberCountAndWeights);
  MITK_TEST(MoveFiber_LabelBundlesTrainAndPredict);
  MITK_TEST(SphereRoi_KeepsFibersOfBothOrientations);
  MITK_TEST(SphereRoi_IsASphere);
  CPPUNIT_TEST_SUITE_END();

public:

  /** Labeling and unlabeling move one fiber at a time between the bundle being labeled and the
   *  label bundles; every bundle stays consistent (fiber count, one weight per fiber) and the
   *  weights travel with their fibers. */
  void MoveFiber_KeepsFiberCountAndWeights()
  {
    std::vector<Fiber> fibers;
    for (int i = 0; i < 5; ++i)
      fibers.push_back(Line({0, static_cast<double>(i), 0}, {39, static_cast<double>(i), 0}));
    auto toLabel = Bundle(fibers, {1, 2, 3, 4, 5});
    auto positive = mitk::FiberBundle::New();
    auto negative = mitk::FiberBundle::New();

    mitk::StreamlineInteractorBrush::MoveFiber(toLabel, positive, 2);
    CPPUNIT_ASSERT_EQUAL(4u, toLabel->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1u, positive->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(3.0f, positive->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(2.0, FirstPoint(positive, 0)[1]);
    const std::vector<float> remaining{1, 2, 4, 5};
    for (unsigned int i = 0; i < 4; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(remaining[i], toLabel->GetFiberWeight(i));
      CPPUNIT_ASSERT_EQUAL(static_cast<double>(remaining[i] - 1), FirstPoint(toLabel, i)[1]);
    }

    mitk::StreamlineInteractorBrush::MoveFiber(toLabel, negative, 0);
    mitk::StreamlineInteractorBrush::MoveFiber(toLabel, positive, 0);
    CPPUNIT_ASSERT_EQUAL(2u, toLabel->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1u, negative->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1.0f, negative->GetFiberWeight(0));
    CPPUNIT_ASSERT_EQUAL(2u, positive->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(2.0f, positive->GetFiberWeight(1));

    // unlabel: back to the bundle being labeled
    mitk::StreamlineInteractorBrush::MoveFiber(positive, toLabel, 0);
    CPPUNIT_ASSERT_EQUAL(1u, positive->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(3u, toLabel->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(3.0f, toLabel->GetFiberWeight(2));

    // exclusion from a prediction: removed only
    mitk::StreamlineInteractorBrush::MoveFiber(toLabel, nullptr, 1);
    CPPUNIT_ASSERT_EQUAL(2u, toLabel->GetNumFibers());

    for (auto *bundle : {toLabel.GetPointer(), positive.GetPointer(), negative.GetPointer()})
    {
      CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(bundle->GetNumFibers()), bundle->GetFiberWeights()->GetNumberOfValues());
      CPPUNIT_ASSERT_EQUAL(bundle->GetNumFibers(), mitk::StreamlineFeatureExtractor::SubsetByIds(bundle, AllIds(bundle))->GetNumFibers());
    }
    CPPUNIT_ASSERT_EQUAL(3u, positive->AddBundle(toLabel)->GetNumFibers());

    CPPUNIT_ASSERT_THROW(mitk::StreamlineInteractorBrush::MoveFiber(toLabel, positive, 2), mitk::Exception);
  }

  /** Label bundles filled by the brush train the classifier, and the prediction node is built from
   *  the predicted pool fibers and the positive labels, as the Interactive Fiber Dissection view
   *  does (CreatePredictionNode). */
  void MoveFiber_LabelBundlesTrainAndPredict()
  {
    auto pool = Bundle(ThreeClusters());
    auto toLabel = mitk::StreamlineFeatureExtractor::SubsetByIds(pool, {0, 1, 2, 3, 4, 5, 60, 61, 62, 120, 121, 122});
    auto positive = mitk::FiberBundle::New();
    auto negative = mitk::FiberBundle::New();
    for (int i = 0; i < 6; ++i)
      mitk::StreamlineInteractorBrush::MoveFiber(toLabel, positive, 0);
    while (toLabel->GetNumFibers() > 0)
      mitk::StreamlineInteractorBrush::MoveFiber(toLabel, negative, 0);

    mitk::StreamlineFeatureExtractor classifier;
    classifier.SetTractogramTest(pool);
    classifier.SetTractogramPrototypes(Prototypes(pool));
    classifier.SetTractogramPlus(positive);
    classifier.SetTractogramMinus(negative);
    classifier.SetActiveCycle(0);
    classifier.Update();

    auto prediction = classifier.ExtractFibers(classifier.GetPredictedIds());
    prediction = prediction->AddBundle(mitk::StreamlineFeatureExtractor::SubsetByIds(positive, AllIds(positive)));

    unsigned int targetFibers = 0;
    for (auto id : classifier.GetPredictedIds())
      targetFibers += id < 60 ? 1 : 0;
    CPPUNIT_ASSERT(targetFibers >= 50);
    CPPUNIT_ASSERT(classifier.GetPredictedIds().size() - targetFibers <= 3);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(classifier.GetPredictedIds().size()) + 6, prediction->GetNumFibers());
  }

  /** A streamline connecting the two ROIs is kept whichever end it starts at. */
  void SphereRoi_KeepsFibersOfBothOrientations()
  {
    auto bundle = Bundle({Line({0, 0, 0}, {40, 0, 0}),     // start ROI -> end ROI
                          Line({40, 0, 0}, {0, 0, 0}),     // end ROI -> start ROI
                          Line({0, 0, 20}, {40, 0, 0})});  // misses the start ROI
    const auto kept = SelectBySphereRois(bundle, SphereRoi({0, 0, 0}, 5), SphereRoi({40, 0, 0}, 5));

    CPPUNIT_ASSERT_EQUAL(std::size_t(2), kept.size());
    std::set<double> startX;
    for (const auto &p : kept)
      startX.insert(p[0]);
    CPPUNIT_ASSERT(startX.count(0) == 1 && startX.count(40) == 1);
  }

  /** A ROI of radius 5 keeps an end point 4.9 mm from its centre and rejects one in the corner of
   *  its bounding box, 6.4 mm from the centre. */
  void SphereRoi_IsASphere()
  {
    auto bundle = Bundle({Line({4.9, 0, 0}, {40, 0, 0}),
                          Line({4.5, 4.5, 0}, {40, 0, 0})});
    const auto kept = SelectBySphereRois(bundle, SphereRoi({0, 0, 0}, 5), SphereRoi({40, 0, 0}, 5));

    CPPUNIT_ASSERT_EQUAL(std::size_t(1), kept.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, kept[0][1], 1e-6);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberDissection)
