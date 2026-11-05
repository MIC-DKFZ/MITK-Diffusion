/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkActiveLearningSession.h>
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

#include <algorithm>
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

  /** A straight fiber from a to b, by default with the point count the classifier needs. */
  Fiber Line(const Point &a, const Point &b, unsigned int n = mitk::StreamlineFeatureExtractor::NumPoints)
  {
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

  /** The tractogram of the session tests, fibers with other point counts than the classifier's:
   *  0-59 the target tract along x (25 points), 60-119 a bundle crossing it in its middle along y
   *  (31 points), 120-179 a bundle that starts where the target starts and turns off to y and z
   *  (57 points). Each fiber has a weight of its own. */
  mitk::FiberBundle::Pointer SessionTractogram()
  {
    std::mt19937 rng(2);
    std::normal_distribution<double> noise(0, 0.7);
    std::vector<Fiber> fibers;
    for (int i = 0; i < 60; ++i)
    {
      const Point s{0, noise(rng), noise(rng)}, e{40, noise(rng), noise(rng)};
      fibers.push_back(i % 2 ? Line(s, e, 25) : Line(e, s, 25));
    }
    for (int i = 0; i < 60; ++i)
    {
      const Point s{20 + noise(rng), -20, noise(rng)}, e{20 + noise(rng), 20, noise(rng)};
      fibers.push_back(i % 2 ? Line(s, e, 31) : Line(e, s, 31));
    }
    for (int i = 0; i < 60; ++i)
    {
      const Point s{noise(rng), noise(rng), noise(rng)}, e{10 + noise(rng), 30 + noise(rng), 15 + noise(rng)};
      fibers.push_back(i % 2 ? Line(s, e, 57) : Line(e, s, 57));
    }
    std::vector<float> weights(fibers.size());
    for (std::size_t i = 0; i < weights.size(); ++i)
      weights[i] = 1.0f + 0.01f * i;
    return Bundle(fibers, weights);
  }

  bool IsTarget(unsigned int originalId) { return originalId < 60; }

  /** Labels the fibers as an oracle that knows the target tract. */
  void LabelByOracle(mitk::ActiveLearningSession &session, const std::vector<unsigned int> &ids)
  {
    std::vector<unsigned int> positive, negative;
    const auto original = session.GetOriginalIds(ids);
    for (std::size_t i = 0; i < ids.size(); ++i)
      (IsTarget(original[i]) ? positive : negative).push_back(ids[i]);
    session.Label(positive, negative);
  }

  /** The tractogram's fibers of the session's prediction. */
  std::set<unsigned int> PredictedOriginals(const mitk::ActiveLearningSession &session)
  {
    const auto ids = session.GetOriginalIds(session.GetPrediction());
    return std::set<unsigned int>(ids.begin(), ids.end());
  }

  /** The pool ids of the tractogram's fibers in [first, last). */
  std::vector<unsigned int> PoolIdsOfOriginals(const mitk::ActiveLearningSession &session, unsigned int first, unsigned int last)
  {
    std::vector<unsigned int> all(session.GetNumberOfPoolFibers());
    std::iota(all.begin(), all.end(), 0);
    const auto original = session.GetOriginalIds(all);
    std::vector<unsigned int> ids;
    for (auto id : all)
    {
      if (original[id] >= first && original[id] < last)
        ids.push_back(id);
    }
    return ids;
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
  MITK_TEST(Classifier_LabeledFibersAreNoCandidates);
  MITK_TEST(Session_OracleLabels_RecoverTheTargetBundle);
  MITK_TEST(Session_FixedSeed_SameSamplesAndPrediction);
  MITK_TEST(Session_FinishedTract_HasTheOriginalFibers);
  MITK_TEST(Session_InitialSample_HalfFromEnrichment);
  MITK_TEST(Session_PoolIdsOf_FindsExtractedFibers);
  MITK_TEST(Session_InvalidLabels_Throw);
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

  /** Labeled pool fibers are not classified, so the prediction node (predicted pool fibers plus
   *  positive labels) holds each fiber once; fibers handed out for labeling but left unlabeled
   *  are classified again in the next cycle. */
  void Classifier_LabeledFibersAreNoCandidates()
  {
    auto pool = Bundle(ThreeClusters());
    const std::vector<unsigned int> plusIds{0, 1, 2, 3, 4, 5};
    const std::vector<unsigned int> minusIds{60, 61, 62, 120, 121, 122};

    mitk::StreamlineFeatureExtractor classifier;
    classifier.SetTractogramTest(pool);
    classifier.SetTractogramPrototypes(Prototypes(pool));
    classifier.SetTractogramPlus(mitk::StreamlineFeatureExtractor::SubsetByIds(pool, plusIds));
    classifier.SetTractogramMinus(mitk::StreamlineFeatureExtractor::SubsetByIds(pool, minusIds));
    classifier.SetActiveCycle(0);
    classifier.Update();

    const std::set<unsigned int> candidates(classifier.GetCandidateIds().begin(), classifier.GetCandidateIds().end());
    const std::set<unsigned int> predicted(classifier.GetPredictedIds().begin(), classifier.GetPredictedIds().end());
    CPPUNIT_ASSERT_EQUAL(pool->GetNumFibers() - 12, static_cast<unsigned int>(candidates.size()));
    for (auto id : plusIds)
    {
      CPPUNIT_ASSERT(candidates.count(id) == 0 && predicted.count(id) == 0);
      CPPUNIT_ASSERT(!classifier.IsAvailableForSampling(id));
    }
    for (auto id : minusIds)
      CPPUNIT_ASSERT(candidates.count(id) == 0);
    CPPUNIT_ASSERT(predicted.count(6) == 1);

    // hand out uncertain fibers, the user labels none of them
    const auto handedOut = classifier.GetDiverseUncertainIds(0.0f, 5);
    CPPUNIT_ASSERT_EQUAL(std::size_t(5), handedOut.size());
    classifier.ExcludeFromSampling(handedOut);
    for (auto id : handedOut)
      CPPUNIT_ASSERT(!classifier.IsAvailableForSampling(id));
    for (auto id : classifier.GetDiverseUncertainIds(0.0f, 20))
      CPPUNIT_ASSERT(std::find(handedOut.begin(), handedOut.end(), id) == handedOut.end());

    classifier.SetActiveCycle(1);
    classifier.Update();
    const std::set<unsigned int> nextCandidates(classifier.GetCandidateIds().begin(), classifier.GetCandidateIds().end());
    CPPUNIT_ASSERT_EQUAL(candidates.size(), nextCandidates.size());
    for (auto id : handedOut)
    {
      CPPUNIT_ASSERT(nextCandidates.count(id) == 1);
      CPPUNIT_ASSERT(classifier.IsAvailableForSampling(id));
    }

    // an unlabeled fiber is a candidate again
    classifier.SetTractogramPlus(mitk::StreamlineFeatureExtractor::SubsetByIds(pool, {0, 1, 2, 3, 4}));
    classifier.SetActiveCycle(2);
    classifier.Update();
    const std::set<unsigned int> afterUnlabel(classifier.GetCandidateIds().begin(), classifier.GetCandidateIds().end());
    CPPUNIT_ASSERT(afterUnlabel.count(5) == 1 && afterUnlabel.count(4) == 0);
  }

  /** A session whose samples are labeled by an oracle finds the target tract, and nothing else,
   *  within three cycles: an initial sample, then the most uncertain fibers. */
  void Session_OracleLabels_RecoverTheTargetBundle()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession session(tractogram, 20, 1);
    CPPUNIT_ASSERT_EQUAL(180u, session.GetNumberOfPoolFibers());
    CPPUNIT_ASSERT_EQUAL(20u, session.GetNumberOfPrototypes());
    CPPUNIT_ASSERT_EQUAL(mitk::StreamlineFeatureExtractor::NumPoints, static_cast<unsigned int>(session.GetPoolFiberPoints(0).size()));

    std::set<unsigned int> target;
    for (unsigned int i = 0; i < 60; ++i)
      target.insert(i);

    auto sample = session.SampleInitial(20);
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), sample.size());
    bool found = false;
    for (unsigned int cycle = 0; cycle < 3 && !found; ++cycle)
    {
      LabelByOracle(session, sample);
      const auto metrics = session.Train();
      CPPUNIT_ASSERT_EQUAL(cycle, metrics.Cycle);
      CPPUNIT_ASSERT_EQUAL(session.GetPrediction().size(), metrics.Predicted);
      CPPUNIT_ASSERT_EQUAL(session.GetPositiveIds().size(), metrics.Positives);
      CPPUNIT_ASSERT_EQUAL(session.GetNegativeIds().size(), metrics.Negatives);
      if (cycle == 0)
        CPPUNIT_ASSERT_EQUAL(1.0, metrics.ChangedFraction);
      CPPUNIT_ASSERT(metrics.ChangedFraction >= 0 && metrics.ChangedFraction <= 1);
      found = PredictedOriginals(session) == target;
      sample = session.SampleUncertain(10);
    }
    CPPUNIT_ASSERT(found);

    // Once the prediction holds, another cycle changes nothing.
    LabelByOracle(session, sample);
    const auto metrics = session.Train();
    CPPUNIT_ASSERT(PredictedOriginals(session) == target);
    CPPUNIT_ASSERT_EQUAL(0.0, metrics.ChangedFraction);
  }

  /** Two sessions with one seed give the same pool, samples and predictions. */
  void Session_FixedSeed_SameSamplesAndPrediction()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession a(tractogram, 20, 7);
    mitk::ActiveLearningSession b(tractogram, 20, 7);

    const auto all = PoolIdsOfOriginals(a, 0, 180);
    CPPUNIT_ASSERT(a.GetOriginalIds(all) == b.GetOriginalIds(all));

    const auto enrichment = PoolIdsOfOriginals(a, 120, 180);
    const auto sampleA = a.SampleInitial(20, enrichment);
    CPPUNIT_ASSERT(sampleA == b.SampleInitial(20, enrichment));

    for (int cycle = 0; cycle < 2; ++cycle)
    {
      LabelByOracle(a, cycle == 0 ? sampleA : a.SampleUncertain(10));
      LabelByOracle(b, cycle == 0 ? sampleA : b.SampleUncertain(10));
      const auto metricsA = a.Train();
      const auto metricsB = b.Train();
      CPPUNIT_ASSERT(a.GetPositiveIds() == b.GetPositiveIds() && a.GetNegativeIds() == b.GetNegativeIds());
      CPPUNIT_ASSERT(a.GetPrediction() == b.GetPrediction());
      CPPUNIT_ASSERT_EQUAL(metricsA.Uncertain, metricsB.Uncertain);
      CPPUNIT_ASSERT_EQUAL(metricsA.ChangedFraction, metricsB.ChangedFraction);
    }
    CPPUNIT_ASSERT(a.SampleDiverse(0.0f, 5) == b.SampleDiverse(0.0f, 5));

    // another seed shuffles the pool otherwise
    mitk::ActiveLearningSession c(tractogram, 20, 8);
    CPPUNIT_ASSERT(a.GetOriginalIds(all) != c.GetOriginalIds(all));
  }

  /** The finished tract holds the predicted fibers as the tractogram holds them: point counts,
   *  coordinates and weights, in the tractogram's order. */
  void Session_FinishedTract_HasTheOriginalFibers()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession session(tractogram, 20, 3);
    LabelByOracle(session, session.SampleInitial(20));
    session.Train();

    auto tract = session.GetFinishedTract();
    auto expected = session.GetOriginalIds(session.GetPrediction());
    std::sort(expected.begin(), expected.end());
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(expected.size()), tract->GetNumFibers());
    CPPUNIT_ASSERT(tract->GetNumFibers() > 0);

    auto original = tractogram->GetFiberPolyData();
    auto result = tract->GetFiberPolyData();
    for (unsigned int i = 0; i < tract->GetNumFibers(); ++i)
    {
      vtkCell *a = original->GetCell(expected[i]);
      const auto numPoints = a->GetNumberOfPoints();
      std::vector<Point> points(numPoints);
      for (vtkIdType j = 0; j < numPoints; ++j)
        a->GetPoints()->GetPoint(j, points[j].data());
      vtkCell *b = result->GetCell(i);
      CPPUNIT_ASSERT_EQUAL(numPoints, b->GetNumberOfPoints());
      CPPUNIT_ASSERT(numPoints != static_cast<vtkIdType>(mitk::StreamlineFeatureExtractor::NumPoints));
      for (vtkIdType j = 0; j < numPoints; ++j)
      {
        Point p;
        b->GetPoints()->GetPoint(j, p.data());
        CPPUNIT_ASSERT(p == points[j]);
      }
      CPPUNIT_ASSERT_EQUAL(tractogram->GetFiberWeight(expected[i]), tract->GetFiberWeight(i));
    }
  }

  /** Half of an initial sample comes from the enrichment, the rest from the whole pool; fibers
   *  handed out are not drawn again before the next training. */
  void Session_InitialSample_HalfFromEnrichment()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession session(tractogram, 20, 4);
    const auto enrichment = PoolIdsOfOriginals(session, 120, 180);
    const std::set<unsigned int> enriched(enrichment.begin(), enrichment.end());

    const auto first = session.SampleInitial(20, enrichment);
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), first.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), std::set<unsigned int>(first.begin(), first.end()).size());
    for (std::size_t i = 0; i < 10; ++i)
      CPPUNIT_ASSERT(enriched.count(first[i]) == 1);

    // a small enrichment: all of it, the rest spread over the pool
    const std::vector<unsigned int> small{enrichment[20], enrichment[21], enrichment[22]};
    const auto second = session.SampleInitial(20, small);
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), second.size());
    CPPUNIT_ASSERT(std::set<unsigned int>(second.begin(), second.begin() + 3) == std::set<unsigned int>(small.begin(), small.end()));
    for (auto id : second)
      CPPUNIT_ASSERT(std::find(first.begin(), first.end(), id) == first.end());

    // without enrichment, the sample covers all three bundles
    const auto third = session.SampleInitial(20);
    std::set<unsigned int> bundles;
    for (auto id : session.GetOriginalIds(third))
      bundles.insert(id / 60);
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), bundles.size());
  }

  /** Fibers extracted from the pool are found in it by their geometry. */
  void Session_PoolIdsOf_FindsExtractedFibers()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession session(tractogram, 20, 5);
    const std::vector<unsigned int> ids{3, 17, 42, 99, 179};
    CPPUNIT_ASSERT(ids == session.GetPoolIdsOf(session.GetPoolSubset(ids)));
    // the original fibers have other point counts, so they are no pool fibers
    CPPUNIT_ASSERT(session.GetPoolIdsOf(tractogram).empty());
  }

  /** Wrong ids and contradicting lists are refused without changing a label; training needs both
   *  classes, sampling uncertain fibers a training. */
  void Session_InvalidLabels_Throw()
  {
    auto tractogram = SessionTractogram();
    mitk::ActiveLearningSession session(tractogram, 20, 6);

    CPPUNIT_ASSERT_THROW(session.Label({180}, {}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(session.Label({1, 2}, {2}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(session.Label({1}, {}, {1}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(session.Label({}, {3}, {3}), mitk::Exception);
    CPPUNIT_ASSERT(session.GetPositiveIds().empty() && session.GetNegativeIds().empty());

    CPPUNIT_ASSERT_THROW(session.SampleUncertain(10), mitk::Exception);
    CPPUNIT_ASSERT_THROW(session.SampleDiverse(0.95f, 10), mitk::Exception);

    session.Label({1, 2}, {});
    CPPUNIT_ASSERT_THROW(session.Train(), mitk::Exception);
    CPPUNIT_ASSERT_EQUAL(0u, session.GetNumberOfTrainings());

    // labeled again: negative, then no label
    session.Label({}, {2});
    CPPUNIT_ASSERT(session.GetPositiveIds() == std::set<unsigned int>{1});
    CPPUNIT_ASSERT(session.GetNegativeIds() == std::set<unsigned int>{2});
    session.Label({}, {}, {2});
    CPPUNIT_ASSERT(session.GetNegativeIds().empty());

    CPPUNIT_ASSERT_THROW(session.GetOriginalIds({180}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(session.GetPoolFiberPoints(180), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::ActiveLearningSession(mitk::FiberBundle::New(), 20), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberDissection)
