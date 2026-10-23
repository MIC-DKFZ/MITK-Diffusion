/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricEuclideanMean.h>
#include <mitkClusteringMetricInnerAngles.h>
#include <mitkClusteringMetricLength.h>
#include <mitkClusteringMetricScalarMap.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <random>

/** Fiber clustering on synthetic bundles whose partition is known. */
class mitkFiberClusteringTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberClusteringTestSuite);
  MITK_TEST(Clusters_KeepInputWeights);
  MITK_TEST(Clusters_SameResultForAnyThreadCount);
  MITK_TEST(InputCentroids_SameResultForAnyThreadCount);
  MITK_TEST(InputCentroids_OutputKBelongsToCentroidK);
  MITK_TEST(Metrics_FibersWithFewerThanFourPoints);
  CPPUNIT_TEST_SUITE_END();

  /** Three fibers along x at y = 0, 0.5, 1 mm and, 40 mm away, two at y = 40, 41 mm, with
   * weights 0.2, 0.3, 0.4, 2 and 3. */
  static mitk::FiberBundle::Pointer TwoGroups()
  {
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({0, 0, 0}, {50, 0, 0}, 20),
                             Line({0, 0.5, 0}, {50, 0.5, 0}, 20),
                             Line({0, 40, 0}, {50, 40, 0}, 20),
                             Line({0, 1, 0}, {50, 1, 0}, 20),
                             Line({0, 41, 0}, {50, 41, 0}, 20)});
    const float weights[] = {0.2f, 0.3f, 2.0f, 0.4f, 3.0f};
    for (unsigned int i = 0; i < 5; ++i)
      fib->SetFiberWeight(i, weights[i]);
    return fib;
  }

  /** Six groups 100 mm apart, each of two sub-bundles 25 mm apart with ten fibers each, half of
   * them reversed, with fixed random offsets and in shuffled order: the hierarchical clustering
   * splits them on several levels. */
  static mitk::FiberBundle::Pointer ManyGroups()
  {
    using namespace mitk::FiberProcessingTestHelper;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> jitter(-1.0, 1.0);
    std::vector<Polyline> fibers;
    for (int group = 0; group < 6; ++group)
      for (int sub = 0; sub < 2; ++sub)
        for (int f = 0; f < 10; ++f)
        {
          const double y = 100.0 * group + 25.0 * sub + jitter(rng);
          const double z = jitter(rng);
          auto line = Line({jitter(rng), y, z}, {40 + jitter(rng), y + jitter(rng), z + jitter(rng)}, 15);
          if (f % 2 == 1)
            std::reverse(line.begin(), line.end());
          fibers.push_back(line);
        }
    std::shuffle(fibers.begin(), fibers.end(), rng);
    return CreateBundle(fibers);
  }

  struct Result
  {
    std::vector<std::vector<unsigned int>> indices;
    std::vector<mitk::FiberBundle::Pointer> centroids;
  };

  static Result Cluster(mitk::FiberBundle::Pointer fib, mitk::FiberBundle::Pointer inCentroids, int threads)
  {
    // The cluster sizes MitkFiberClustering uses for a 10 mm cluster size on this bundle.
    std::vector<float> distances;
    for (int i = 1; i <= 26; ++i)
      distances.push_back(10.0f * i);

    const int previousThreads = omp_get_max_threads();
    omp_set_num_threads(threads);
    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(fib);
    if (inCentroids.IsNotNull())
      clusterer.SetInCentroids(inCentroids);
    clusterer.SetDistances(distances);
    clusterer.SetMetrics({new mitk::ClusteringMetricEuclideanMean()});
    clusterer.Update();
    omp_set_num_threads(previousThreads);
    return {clusterer.GetOutFiberIndices(), clusterer.GetOutCentroids()};
  }

  static void AssertEqualResults(const Result& expected, const Result& actual)
  {
    CPPUNIT_ASSERT(expected.indices == actual.indices);
    CPPUNIT_ASSERT_EQUAL(expected.centroids.size(), actual.centroids.size());
    for (std::size_t c = 0; c < expected.centroids.size(); ++c)
      CPPUNIT_ASSERT(expected.centroids[c]->Equals(actual.centroids[c], 1e-5));
  }

public:
  void Clusters_SameResultForAnyThreadCount()
  {
    auto fib = ManyGroups();
    const auto single = Cluster(fib, nullptr, 1);
    CPPUNIT_ASSERT_EQUAL(std::size_t(12), single.indices.size());
    for (int run = 0; run < 3; ++run)
      AssertEqualResults(single, Cluster(fib, nullptr, 8));
  }

  void InputCentroids_SameResultForAnyThreadCount()
  {
    auto fib = ManyGroups();
    auto centroids = mitk::FiberBundle::New()->AddBundles(Cluster(fib, nullptr, 1).centroids);
    const auto single = Cluster(fib, centroids, 1);
    for (int run = 0; run < 3; ++run)
      AssertEqualResults(single, Cluster(fib, centroids, 8));
  }

  void InputCentroids_OutputKBelongsToCentroidK()
  {
    // Centroid 0 lies on the group at y = 40, centroid 1 far from every fiber, centroid 2 on the
    // group at y = 0; a sixth fiber at y = 200 is far from all of them.
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = TwoGroups();
    fib = fib->AddBundle(CreateBundle({Line({0, 200, 0}, {50, 200, 0}, 20)}));
    auto centroids = CreateBundle({Line({0, 40.5, 0}, {50, 40.5, 0}, 12), Line({0, -100, 0}, {50, -100, 0}, 12), Line({50, 0.5, 0}, {0, 0.5, 0}, 12)});

    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(fib);
    clusterer.SetInCentroids(centroids);
    clusterer.SetDistances({10});
    clusterer.SetMetrics({new mitk::ClusteringMetricEuclideanMean()});
    clusterer.Update();

    const auto indices = clusterer.GetOutFiberIndices();
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), clusterer.GetOutTractograms().size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), clusterer.GetOutCentroids().size());
    CPPUNIT_ASSERT(indices.at(0) == std::vector<unsigned int>({2, 4}));
    CPPUNIT_ASSERT(indices.at(1).empty());
    CPPUNIT_ASSERT_EQUAL(0u, clusterer.GetOutTractograms().at(1)->GetNumFibers());
    CPPUNIT_ASSERT(indices.at(2) == std::vector<unsigned int>({0, 1, 3}));
    CPPUNIT_ASSERT_EQUAL(1u, clusterer.GetOutUnassigned()->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-100.0, clusterer.GetOutCentroids().at(1)->GetFiberPolyData()->GetPoint(0)[1], 1e-4);
  }

  void Metrics_FibersWithFewerThanFourPoints()
  {
    // The metrics step through every fourth point; with fewer than four points the step was 0
    // and the loop never ended.
    vnl_matrix<float> s(3, 3);
    vnl_matrix<float> t(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
    {
      s.set_column(i, vnl_vector_fixed<float, 3>(2.0f * i, 0, 0).as_vector());
      t.set_column(i, vnl_vector_fixed<float, 3>(2.0f * i, 1, 0).as_vector());
    }

    auto map = itk::Image<float, 3>::New();
    itk::Image<float, 3>::SizeType size = {{6, 3, 3}};
    map->SetRegions(itk::Image<float, 3>::RegionType(size));
    map->Allocate();
    map->FillBuffer(1);

    mitk::ClusteringMetricLength length;
    mitk::ClusteringMetricInnerAngles angles;
    mitk::ClusteringMetricScalarMap scalarMap;
    scalarMap.SetImages({map});
    bool flipped = true;
    CPPUNIT_ASSERT(std::isfinite(length.CalculateDistance(s, t, flipped)));
    CPPUNIT_ASSERT(!flipped);
    CPPUNIT_ASSERT(std::isfinite(angles.CalculateDistance(s, t, flipped)));
    CPPUNIT_ASSERT(std::isfinite(scalarMap.CalculateDistance(s, t, flipped)));
  }

  void Clusters_KeepInputWeights()
  {
    auto input = TwoGroups();
    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(input);
    clusterer.SetDistances({10});
    clusterer.SetMetrics({new mitk::ClusteringMetricEuclideanMean()});
    clusterer.Update();

    const auto tracts = clusterer.GetOutTractograms();
    const auto indices = clusterer.GetOutFiberIndices();
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), tracts.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), indices.at(0).size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), indices.at(1).size());
    for (std::size_t c = 0; c < tracts.size(); ++c)
    {
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(indices.at(c).size()), tracts.at(c)->GetNumFibers());
      for (unsigned int i = 0; i < tracts.at(c)->GetNumFibers(); ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(input->GetFiberWeight(indices.at(c).at(i)), tracts.at(c)->GetFiberWeight(i), 1e-6);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberClustering)
