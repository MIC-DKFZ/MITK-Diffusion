/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricAnatomic.h>
#include <mitkClusteringMetricEuclideanMax.h>
#include <mitkClusteringMetricEuclideanMean.h>
#include <mitkClusteringMetricEuclideanStd.h>
#include <mitkClusteringMetricInnerAngles.h>
#include <mitkClusteringMetricLength.h>
#include <mitkClusteringMetricScalarMap.h>
#include <mitkIOUtil.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkExceptionMacro.h>

#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>

#include "mitkFiberProcessingTestHelper.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <string>

/** Fiber clustering on synthetic bundles whose partition is known. */
class mitkFiberClusteringTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberClusteringTestSuite);
  MITK_TEST(Clusters_KeepInputWeights);
  MITK_TEST(Clusters_SameResultForAnyThreadCount);
  MITK_TEST(InputCentroids_SameResultForAnyThreadCount);
  MITK_TEST(InputCentroids_OutputKBelongsToCentroidK);
  MITK_TEST(Metrics_FibersWithFewerThanFourPoints);
  MITK_TEST(Resampling_EqualArcLength);
  MITK_TEST(ScalarMap_SeveralMapsInParallel);
  MITK_TEST(Anatomic_PointOutsideFirstParcellation);
  MITK_TEST(Update_TwiceGivesSameOutputs);
  MITK_TEST(Metrics_Definitions);
  MITK_TEST(Metrics_SameValuesAsVnlFormulation);
  MITK_TEST(Outputs_KeepRasFlagAndSetHeaderOnly);
  MITK_TEST(FiberFitBundles_PartitionAsQuickBundles);
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
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
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
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
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
    // The metrics step through every fourth point; fewer than four points must not give a step
    // of 0, which would never end the loop.
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

  void Resampling_EqualArcLength()
  {
    // The same 100 mm line, once with 12 points crowded at its start, once with 50 even points:
    // after resampling to 12 points at equal arc length both are the same fiber.
    using namespace mitk::FiberProcessingTestHelper;
    Polyline crowded;
    for (int i = 0; i <= 10; ++i)
      crowded.push_back({static_cast<double>(i), 0, 0});
    crowded.push_back({100, 0, 0});
    auto fib = CreateBundle({crowded, Line({0, 0.5, 0}, {100, 0.5, 0}, 50)});

    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(fib);
    clusterer.SetDistances({10});
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
    clusterer.Update();

    CPPUNIT_ASSERT_EQUAL(std::size_t(1), clusterer.GetOutTractograms().size());
    auto centroid = clusterer.GetOutCentroids().at(0)->GetFiberPolyData();
    CPPUNIT_ASSERT_EQUAL(vtkIdType(12), centroid->GetNumberOfPoints());
    for (vtkIdType k = 0; k < 12; ++k)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(100.0 * k / 11, centroid->GetPoint(k)[0], 1e-3);
  }

  void ScalarMap_SeveralMapsInParallel()
  {
    // Two maps (value = x and value = 10 y); distances computed by 8 threads must equal the
    // sequential ones. A shared interpolator switched between the maps gave other values.
    using ImageType = itk::Image<float, 3>;
    std::vector<ImageType::Pointer> maps;
    for (int m = 0; m < 2; ++m)
    {
      auto map = ImageType::New();
      ImageType::SizeType size = {{20, 20, 3}};
      map->SetRegions(ImageType::RegionType(size));
      map->Allocate();
      itk::ImageRegionIteratorWithIndex<ImageType> it(map, map->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
        it.Set(m == 0 ? it.GetIndex()[0] : 10.0f * it.GetIndex()[1]);
      maps.push_back(map);
    }
    mitk::ClusteringMetricScalarMap metric;
    metric.SetImages(maps);

    std::mt19937 rng(7);
    std::uniform_real_distribution<float> coordinate(0.0f, 19.0f);
    const int numPairs = 2000;
    std::vector<vnl_matrix<float>> s(numPairs, vnl_matrix<float>(3, 12)), t(numPairs, vnl_matrix<float>(3, 12));
    for (int i = 0; i < numPairs; ++i)
      for (unsigned int k = 0; k < 12; ++k)
      {
        s[i](0, k) = coordinate(rng); s[i](1, k) = coordinate(rng); s[i](2, k) = 1;
        t[i](0, k) = coordinate(rng); t[i](1, k) = coordinate(rng); t[i](2, k) = 1;
      }

    std::vector<float> sequential(numPairs), parallel(numPairs);
    bool flipped = false;
    for (int i = 0; i < numPairs; ++i)
      sequential[i] = metric.CalculateDistance(s[i], t[i], flipped);

    const int previousThreads = omp_get_max_threads();
    omp_set_num_threads(8);
#pragma omp parallel for
    for (int i = 0; i < numPairs; ++i)
    {
      bool f = false;
      parallel[i] = metric.CalculateDistance(s[i], t[i], f);
    }
    omp_set_num_threads(previousThreads);
    for (int i = 0; i < numPairs; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sequential[i], parallel[i], 1e-3);
  }

  void Anatomic_PointOutsideFirstParcellation()
  {
    // Parcellation A (labels 1, 2) covers [0, 2) mm, parcellation B (labels 5, 6, 7) starts at
    // 10 mm. Fibers inside B only must be looked up in B's labels, not in A's.
    using ImageType = itk::Image<short, 3>;
    const auto createParcellation = [](double origin, unsigned int size, short firstLabel, short numLabels) {
      auto parc = ImageType::New();
      ImageType::SizeType imageSize = {{size, size, size}};
      parc->SetRegions(ImageType::RegionType(imageSize));
      parc->SetOrigin(itk::Point<double, 3>(origin));
      parc->Allocate();
      itk::ImageRegionIteratorWithIndex<ImageType> it(parc, parc->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
        it.Set(firstLabel + it.GetIndex()[0] % numLabels);
      return parc;
    };

    mitk::ClusteringMetricAnatomic metric;
    metric.SetParcellations({createParcellation(0, 2, 1, 2), createParcellation(10, 6, 5, 3)});

    vnl_matrix<float> s(3, 12);
    vnl_matrix<float> t(3, 12);
    for (unsigned int k = 0; k < 12; ++k)
    {
      s.set_column(k, vnl_vector_fixed<float, 3>(10.5f + 0.3f * k, 12, 12).as_vector());
      t.set_column(k, vnl_vector_fixed<float, 3>(10.5f + 0.3f * k, 12.2f, 12).as_vector());
    }
    bool flipped = true;
    float distance = 0;
    CPPUNIT_ASSERT_NO_THROW(distance = metric.CalculateDistance(s, t, flipped));
    CPPUNIT_ASSERT(distance < 9999);
  }

  void Update_TwiceGivesSameOutputs()
  {
    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(TwoGroups());
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
    CPPUNIT_ASSERT_THROW(clusterer.Update(), mitk::Exception);

    clusterer.SetDistances({10});
    clusterer.Update();
    const auto first = clusterer.GetOutFiberIndices();
    clusterer.ClearMetrics();
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
    clusterer.Update();
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), clusterer.GetOutTractograms().size());
    CPPUNIT_ASSERT(first == clusterer.GetOutFiberIndices());
  }

  static vnl_matrix<float> Points(const std::vector<std::array<float, 3>>& points)
  {
    vnl_matrix<float> m(3, points.size());
    for (unsigned int k = 0; k < points.size(); ++k)
      for (unsigned int d = 0; d < 3; ++d)
        m(d, k) = points[k][d];
    return m;
  }

  void Metrics_Definitions()
  {
    bool flipped = true;

    // EU_STD: point distances 1, 2, 3, 4 mm: (mean 2.5 + standard deviation sqrt(1.25)) / 2.
    auto s = Points({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}});
    auto t = Points({{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 4, 0}});
    mitk::ClusteringMetricEuclideanStd std;
    CPPUNIT_ASSERT_DOUBLES_EQUAL((2.5 + std::sqrt(1.25)) / 2, std.CalculateDistance(s, t, flipped), 1e-5);
    CPPUNIT_ASSERT(!flipped);

    // LENGTH: straight fibers of 11 and 16.5 mm with 12 points.
    std::vector<std::array<float, 3>> a, b;
    for (int k = 0; k < 12; ++k)
    {
      a.push_back({static_cast<float>(k), 0, 0});
      b.push_back({1.5f * k, 1, 0});
    }
    auto la = Points(a);
    auto lb = Points(b);
    mitk::ClusteringMetricLength length;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.5, length.CalculateDistance(la, lb, flipped), 1e-4);

    // ANGLES: straight fibers in random directions have no NaN angle differences.
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> coordinate(-50.0f, 50.0f);
    mitk::ClusteringMetricInnerAngles angles;
    for (int f = 0; f < 200; ++f)
    {
      std::array<float, 3> start = {coordinate(rng), coordinate(rng), coordinate(rng)};
      std::array<float, 3> end = {coordinate(rng), coordinate(rng), coordinate(rng)};
      std::vector<std::array<float, 3>> line;
      for (int k = 0; k < 12; ++k)
        line.push_back({start[0] + (end[0] - start[0]) * k / 11.0f, start[1] + (end[1] - start[1]) * k / 11.0f, start[2] + (end[2] - start[2]) * k / 11.0f});
      auto m = Points(line);
      CPPUNIT_ASSERT(std::isfinite(angles.CalculateDistance(m, m, flipped)));
    }
  }

  // #35: the metrics compute without per-comparison allocations, the anatomic and scalar map
  // metrics from per-fiber data prepared once. These are the former vnl formulations; the
  // clustering results depend on the metric values bit for bit.
  static float VnlDistance(const vnl_matrix<float>& s, unsigned int i, const vnl_matrix<float>& t, unsigned int j)
  {
    return (s.get_column(i) - t.get_column(j)).magnitude();
  }

  static float VnlEuclideanMean(float scale, const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool& flipped)
  {
    float d_direct = 0, d_flipped = 0;
    for (unsigned int i = 0; i < s.cols(); ++i)
    {
      d_direct += VnlDistance(s, i, t, i);
      d_flipped += VnlDistance(s, i, t, s.cols() - i - 1);
    }
    flipped = d_direct > d_flipped;
    return scale * (flipped ? d_flipped : d_direct) / s.cols();
  }

  static float VnlEuclideanStdOrMax(float scale, const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool max, bool& flipped)
  {
    float d_direct = 0, d_flipped = 0;
    vnl_vector<float> dists_d(s.cols()), dists_f(s.cols());
    for (unsigned int i = 0; i < s.cols(); ++i)
    {
      dists_d[i] = VnlDistance(s, i, t, i);
      d_direct += dists_d[i];
      dists_f[i] = VnlDistance(s, i, t, s.cols() - i - 1);
      d_flipped += dists_f[i];
    }
    flipped = d_direct > d_flipped;
    vnl_vector<float>& dists = flipped ? dists_f : dists_d;
    if (max)
      return scale * dists.max_value();
    float d = (flipped ? d_flipped : d_direct) / s.cols();
    dists -= d;
    d += dists.magnitude() / std::sqrt(static_cast<float>(s.cols()));
    return scale * d / 2;
  }

  static float VnlLength(float scale, const vnl_matrix<float>& s, const vnl_matrix<float>& t)
  {
    float l1 = 0, l2 = 0;
    for (unsigned int i = 0; i + 1 < s.cols(); ++i)
    {
      l1 += (s.get_column(i + 1) - s.get_column(i)).magnitude();
      l2 += (t.get_column(i + 1) - t.get_column(i)).magnitude();
    }
    return scale * std::fabs(l1 - l2);
  }

  static float VnlAnatomic(const mitk::ClusteringMetricAnatomic& metric, unsigned int numLabels, const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool& flipped)
  {
    vnl_matrix<float> hist1(27, numLabels, 0.0f), hist2(27, numLabels, 0.0f);
    float d_direct = 0, d_flipped = 0;
    for (unsigned int i = 0; i < s.cols(); ++i)
    {
      metric.GetHistogramAtPoint(itk::Point<float, 3>(s.get_column(i).data_block()), hist1);
      metric.GetHistogramAtPoint(itk::Point<float, 3>(t.get_column(i).data_block()), hist2);
      d_direct += VnlDistance(s, i, t, i);
      d_flipped += VnlDistance(s, i, t, s.cols() - i - 1);
    }
    flipped = d_direct > d_flipped;

    float label_intersection = 0;
    for (unsigned int c = 0; c < hist1.cols(); ++c)
    {
      bool l1 = false, l2 = false;
      for (unsigned int r = 0; r < hist1.rows(); ++r)
      {
        l1 = l1 || hist1[r][c] > 0;
        l2 = l2 || hist2[r][c] > 0;
        if (l1 && l2)
        {
          label_intersection += 1;
          break;
        }
      }
    }
    float similarity = 0;
    if (label_intersection > 0)
    {
      hist1.normalize_rows();
      hist2.normalize_rows();
      label_intersection /= static_cast<short>(numLabels);
      for (unsigned int l = 0; l < hist1.rows(); ++l)
        similarity += dot_product(hist1.get_row(l), hist2.get_row(l));
      similarity /= hist1.rows();
      similarity *= label_intersection;
    }
    return similarity > 0 ? metric.GetScale() * 0.2 / similarity : 9999;
  }

  using FloatInterpolator = itk::LinearInterpolateImageFunction<itk::Image<float, 3>, float>;

  static float VnlScalarMap(const std::vector<FloatInterpolator::Pointer>& interpolators, float scale, const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool& flipped)
  {
    auto values = [&interpolators](const vnl_matrix<float>& m, unsigned int i) {
      vnl_vector<float> vals(static_cast<unsigned int>(interpolators.size()));
      for (unsigned int c = 0; c < interpolators.size(); ++c)
        vals[c] = mitk::imv::GetImageValue<float>(itk::Point<float, 3>(m.get_column(i).data_block()), true, interpolators[c]);
      return vals;
    };
    float d_direct = 0, d_flipped = 0;
    const unsigned int inc = s.cols() >= 4 ? s.cols() / 4 : 1;
    for (unsigned int i = 0; i < s.cols(); i += inc)
    {
      d_direct += VnlDistance(s, i, t, i);
      d_flipped += VnlDistance(s, i, t, s.cols() - i - 1);
    }
    flipped = d_direct > d_flipped;
    float map_distance = 0;
    for (unsigned int i = 0; i < s.cols(); ++i)
      map_distance += (values(s, i) - values(t, flipped ? s.cols() - i - 1 : i)).magnitude();
    return scale * map_distance;
  }

  void Metrics_SameValuesAsVnlFormulation()
  {
    // Two parcellations (labels 0..6 and 20..23, 1 mm voxels) and two float maps, overlapping;
    // fibers partly outside them.
    using ShortImage = itk::Image<short, 3>;
    using FloatImage = itk::Image<float, 3>;
    std::mt19937 rng(11);
    auto parcellation = [&rng](double origin, short firstLabel, short numLabels) {
      auto parc = ShortImage::New();
      parc->SetRegions(ShortImage::RegionType(ShortImage::SizeType{{12, 12, 12}}));
      parc->SetOrigin(itk::Point<double, 3>(origin));
      parc->Allocate();
      std::uniform_int_distribution<int> label(0, numLabels - 1);
      for (itk::ImageRegionIteratorWithIndex<ShortImage> it(parc, parc->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
        it.Set(static_cast<short>(firstLabel + label(rng)));
      return parc;
    };
    mitk::ClusteringMetricAnatomic anatomic;
    anatomic.SetParcellations({parcellation(0, 0, 7), parcellation(4, 20, 4)});
    anatomic.SetScale(1.5f);
    const unsigned int numLabels = 11;

    std::vector<FloatImage::Pointer> maps;
    std::vector<FloatInterpolator::Pointer> interpolators;
    std::uniform_real_distribution<float> value(-5.0f, 5.0f);
    for (int m = 0; m < 2; ++m)
    {
      auto map = FloatImage::New();
      map->SetRegions(FloatImage::RegionType(FloatImage::SizeType{{14, 14, 14}}));
      map->SetOrigin(itk::Point<double, 3>(m * 2.0));
      map->Allocate();
      for (itk::ImageRegionIteratorWithIndex<FloatImage> it(map, map->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
        it.Set(value(rng));
      maps.push_back(map);
      auto interpolator = FloatInterpolator::New();
      interpolator->SetInputImage(map);
      interpolators.push_back(interpolator);
    }
    mitk::ClusteringMetricScalarMap scalarMap;
    scalarMap.SetImages(maps);

    mitk::ClusteringMetricEuclideanMean mean;
    mitk::ClusteringMetricEuclideanStd std;
    mitk::ClusteringMetricEuclideanMax max;
    mitk::ClusteringMetricLength length;
    mean.SetScale(1.25f);
    std.SetScale(0.75f);

    std::uniform_real_distribution<float> coordinate(-2.0f, 18.0f);
    for (const unsigned int numPoints : {12u, 3u, 1u, 40u})
      for (int pair = 0; pair < 50; ++pair)
      {
        vnl_matrix<float> s(3, numPoints), t(3, numPoints);
        for (unsigned int k = 0; k < numPoints; ++k)
          for (unsigned int d = 0; d < 3; ++d)
          {
            s(d, k) = coordinate(rng);
            t(d, k) = coordinate(rng);
          }
        // the second fiber close to the first one reversed in every other pair: flips both ways
        if (pair % 2 == 1)
          for (unsigned int k = 0; k < numPoints; ++k)
            for (unsigned int d = 0; d < 3; ++d)
              t(d, k) = s(d, numPoints - 1 - k) + 0.1f * coordinate(rng);

        const std::string message = std::to_string(numPoints) + " points, pair " + std::to_string(pair);
        bool flipped = false, expectedFlip = false;
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_MEAN, " + message, VnlEuclideanMean(1.25f, s, t, expectedFlip), mean.CalculateDistance(s, t, flipped));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_MEAN flip, " + message, expectedFlip, flipped);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_STD, " + message, VnlEuclideanStdOrMax(0.75f, s, t, false, expectedFlip), std.CalculateDistance(s, t, flipped));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_STD flip, " + message, expectedFlip, flipped);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_MAX, " + message, VnlEuclideanStdOrMax(1.0f, s, t, true, expectedFlip), max.CalculateDistance(s, t, flipped));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("EU_MAX flip, " + message, expectedFlip, flipped);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("LENGTH, " + message, VnlLength(1.0f, s, t), length.CalculateDistance(s, t, flipped));

        const float expectedAnatomic = VnlAnatomic(anatomic, numLabels, s, t, expectedFlip);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("ANAT, " + message, expectedAnatomic, anatomic.CalculateDistance(s, t, flipped));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("ANAT flip, " + message, expectedFlip, flipped);
        const auto sAnatomic = anatomic.PrepareFiber(s);
        const auto tAnatomic = anatomic.PrepareFiber(t);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("ANAT prepared, " + message, expectedAnatomic,
                                     anatomic.CalculatePreparedDistance(s, t, sAnatomic.get(), tAnatomic.get(), flipped));

        const float expectedMap = VnlScalarMap(interpolators, scalarMap.GetScale(), s, t, expectedFlip);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("MAP, " + message, expectedMap, scalarMap.CalculateDistance(s, t, flipped));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("MAP flip, " + message, expectedFlip, flipped);
      }
  }

  void Clusters_KeepInputWeights()
  {
    auto input = TwoGroups();
    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(input);
    clusterer.SetDistances({10});
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
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

  void Outputs_KeepRasFlagAndSetHeaderOnly()
  {
    // Clusters, centroids and unassigned fibers are bundles derived from the input: they keep its
    // coordinate convention and a TrackVis header set on it, but not the header the input derives
    // from its own bounding box, which would not contain a cluster's fibers.
    using namespace mitk::FiberProcessingTestHelper;
    auto plain = TwoGroups();
    plain = plain->AddBundle(CreateBundle({Line({0, 200, 0}, {50, 200, 0}, 20)}));
    auto withMetadata = plain->GetDeepCopy();
    SetMetadata(withMetadata);
    auto centroids = CreateBundle({Line({0, 40.5, 0}, {50, 40.5, 0}, 12), Line({50, 0.5, 0}, {0, 0.5, 0}, 12)});

    for (auto input : {plain, withMetadata})
      for (auto inCentroids : {mitk::FiberBundle::Pointer(), centroids})
      {
        mitk::TractClusteringFilter clusterer;
        clusterer.SetTractogram(input);
        clusterer.SetInCentroids(inCentroids);
        clusterer.SetDistances({10});
        clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
        clusterer.Update();

        CPPUNIT_ASSERT(!clusterer.GetOutTractograms().empty());
        for (const auto& tract : clusterer.GetOutTractograms())
          CPPUNIT_ASSERT(HasMetadataOf(tract, input));
        for (const auto& centroid : clusterer.GetOutCentroids())
          CPPUNIT_ASSERT(HasMetadataOf(centroid, input));
        if (inCentroids.IsNotNull())
        {
          CPPUNIT_ASSERT_EQUAL(1u, clusterer.GetOutUnassigned()->GetNumFibers());
          CPPUNIT_ASSERT(HasMetadataOf(clusterer.GetOutUnassigned(), input));
        }
      }
  }

  void FiberFitBundles_PartitionAsQuickBundles()
  {
    // Safety net on mitk-diffusion-test-data: the five FiberFit bundles (2157 fibers) joined and clustered as
    // MitkFiberClustering does with its defaults (10 mm, EU_MEAN, 12 points, duplicate merging
    // at half the cluster size). dipy's QuickBundles gives the same partition into 7 clusters.
    std::vector<mitk::FiberBundle::Pointer> bundles;
    for (int i = 0; i < 5; ++i)
      bundles.push_back(mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("FiberFit/Cluster_" + std::to_string(i) + ".fib")));
    auto fib = mitk::FiberBundle::New()->AddBundles(bundles);
    CPPUNIT_ASSERT_EQUAL(2157u, fib->GetNumFibers());

    std::vector<float> distances;
    for (int i = 1; distances.empty() || distances.back() < fib->GetGeometry()->GetDiagonalLength() / 2; ++i)
      distances.push_back(10.0f * i);

    mitk::TractClusteringFilter clusterer;
    clusterer.SetTractogram(fib);
    clusterer.SetDistances(distances);
    clusterer.AddMetric(std::make_unique<mitk::ClusteringMetricEuclideanMean>());
    clusterer.Update();

    std::vector<unsigned int> sizes;
    for (const auto& tract : clusterer.GetOutTractograms())
      sizes.push_back(tract->GetNumFibers());
    CPPUNIT_ASSERT(sizes == std::vector<unsigned int>({649, 508, 400, 302, 276, 14, 8}));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberClustering)
