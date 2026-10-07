/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "mitkStreamlineFeatureExtractor.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include <vtkFloatArray.h>
#include <vtkPolyData.h>

namespace mitk{

namespace
{
  /** Pools larger than this are reduced to the candidates closest to the positive labels in the first cycle. */
  constexpr std::size_t PoolReductionThreshold = 80000;
  constexpr double PoolReductionFraction = 0.2;
  /** Number of the most uncertain fibers considered by GetDiverseUncertainIds(). */
  constexpr std::size_t MaxDiverseCandidates = 250;

  double PointDistance(const float *a, const float *b)
  {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  double BinaryEntropy(double p)
  {
    if (p <= 0.0 || p >= 1.0)
      return 0.0;
    return -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
  }

  /**
  * Farthest-first traversal over 'size' elements starting at 'first': repeatedly picks the element with the
  * largest distance to all elements picked so far. distance(i, j) is evaluated on the fly.
  */
  template <typename DistanceFunction>
  std::vector<unsigned int> FarthestFirst(std::size_t size, std::size_t count, unsigned int first, DistanceFunction distance)
  {
    std::vector<unsigned int> selected;
    count = std::min(count, size);
    if (count == 0)
      return selected;

    std::vector<float> minDistance(size, std::numeric_limits<float>::max());
    unsigned int next = first;
    while (selected.size() < count)
    {
      selected.push_back(next);
      minDistance[next] = -1;

#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(size); ++i)
      {
        if (minDistance[i] >= 0)
          minDistance[i] = std::min(minDistance[i], static_cast<float>(distance(next, i)));
      }
      next = static_cast<unsigned int>(std::max_element(minDistance.begin(), minDistance.end()) - minDistance.begin());
    }
    return selected;
  }
}

StreamlineFeatureExtractor::StreamlineFeatureExtractor()
  : m_ActiveCycle(0)
{
}

StreamlineFeatureExtractor::~StreamlineFeatureExtractor()
{
}

void StreamlineFeatureExtractor::SetTractogramTest(const mitk::FiberBundle::Pointer &tractogram)
{
  m_TractogramTest = tractogram;
  m_PoolFibers = ToMatrix(tractogram, "The tractogram");
  m_PoolIds.resize(m_PoolFibers.size() / FiberStride);
  std::iota(m_PoolIds.begin(), m_PoolIds.end(), 0);
  m_PoolActive.assign(m_PoolIds.size(), 1);
  m_PrototypeColumns.clear();
  m_LocalColumns.clear();
}

void StreamlineFeatureExtractor::SetTractogramPrototypes(const mitk::FiberBundle::Pointer &prototypes)
{
  m_Prototypes = ToMatrix(prototypes, "The prototype bundle");
  m_PrototypeColumns.clear();
}

void StreamlineFeatureExtractor::SetTractogramPlus(const mitk::FiberBundle::Pointer &tractogram)
{
  m_TractogramPlus = tractogram;
}

void StreamlineFeatureExtractor::SetTractogramMinus(const mitk::FiberBundle::Pointer &tractogram)
{
  m_TractogramMinus = tractogram;
}

void StreamlineFeatureExtractor::SetActiveCycle(int activeCycle)
{
  m_ActiveCycle = activeCycle;
}

StreamlineFeatureExtractor::FiberMatrix StreamlineFeatureExtractor::ToMatrix(mitk::FiberBundle *fib, const std::vector<unsigned int> &ids, const std::string &name)
{
  FiberMatrix out(ids.size() * FiberStride);
  vtkPolyData* polyData = fib->GetFiberPolyData();
  for (std::size_t i = 0; i < ids.size(); ++i)
  {
    vtkCell* cell = polyData->GetCell(ids[i]);
    if (cell->GetNumberOfPoints() != static_cast<vtkIdType>(NumPoints))
      throw std::runtime_error(name + " is not resampled to " + std::to_string(NumPoints) + " points per streamline. Resample it in the preprocessing step and reset the classifier.");

    vtkPoints* points = cell->GetPoints();
    float* dst = &out[i * FiberStride];
    for (unsigned int j = 0; j < NumPoints; ++j)
    {
      double p[3];
      points->GetPoint(j, p);
      dst[3*j] = p[0];
      dst[3*j + 1] = p[1];
      dst[3*j + 2] = p[2];
    }
  }
  return out;
}

StreamlineFeatureExtractor::FiberMatrix StreamlineFeatureExtractor::ToMatrix(mitk::FiberBundle *fib, const std::string &name)
{
  if (fib == nullptr)
    return FiberMatrix();

  std::vector<unsigned int> ids(fib->GetFiberPolyData()->GetNumberOfCells());
  std::iota(ids.begin(), ids.end(), 0);
  return ToMatrix(fib, ids, name);
}

std::uint64_t StreamlineFeatureExtractor::Hash(const float *fiber)
{
  // FNV-1a over the coordinates. Labeled fibers are exact copies of pool fibers, so identical geometry is a stable key.
  const auto* bytes = reinterpret_cast<const unsigned char*>(fiber);
  std::uint64_t hash = 14695981039346656037ull;
  for (std::size_t i = 0; i < FiberStride * sizeof(float); ++i)
  {
    hash ^= bytes[i];
    hash *= 1099511628211ull;
  }
  return hash;
}

void StreamlineFeatureExtractor::Distances(const float *a, const float *b, float &mdf, float &endpoints)
{
  const unsigned int last = NumPoints - 1;

  double direct = 0;
  double flipped = 0;
  for (unsigned int j = 0; j < NumPoints; ++j)
  {
    direct += PointDistance(a + 3*j, b + 3*j);
    flipped += PointDistance(a + 3*j, b + 3*(last - j));
  }
  mdf = std::min(direct, flipped) / NumPoints;

  const double endDirect = (PointDistance(a, b) + PointDistance(a + 3*last, b + 3*last)) / 2;
  const double endFlipped = (PointDistance(a, b + 3*last) + PointDistance(a + 3*last, b)) / 2;
  endpoints = std::min(endDirect, endFlipped);
}

std::vector<StreamlineFeatureExtractor::FeatureColumn> StreamlineFeatureExtractor::ComputeColumns(const FiberMatrix &references) const
{
  const std::size_t numReferences = references.size() / FiberStride;
  const std::size_t numCandidates = m_PoolIds.size();

  std::vector<FeatureColumn> columns(numReferences);
  for (auto& column : columns)
  {
    column.mdf.resize(numCandidates);
    column.endpoints.resize(numCandidates);
  }

#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(numCandidates); ++i)
  {
    const float* fiber = &m_PoolFibers[i * FiberStride];
    for (std::size_t r = 0; r < numReferences; ++r)
      Distances(fiber, &references[r * FiberStride], columns[r].mdf[i], columns[r].endpoints[i]);
  }
  return columns;
}

void StreamlineFeatureExtractor::AddLocalColumns(const std::vector<std::uint64_t> &keys, const FiberMatrix &fibers)
{
  std::unordered_set<std::uint64_t> pending;
  std::vector<std::uint64_t> missingKeys;
  FiberMatrix missingFibers;
  for (std::size_t i = 0; i < keys.size(); ++i)
  {
    if (m_LocalColumns.count(keys[i]) || !pending.insert(keys[i]).second)
      continue;
    missingKeys.push_back(keys[i]);
    missingFibers.insert(missingFibers.end(), fibers.begin() + i * FiberStride, fibers.begin() + (i + 1) * FiberStride);
  }

  auto columns = ComputeColumns(missingFibers);
  for (std::size_t i = 0; i < missingKeys.size(); ++i)
    m_LocalColumns[missingKeys[i]] = std::move(columns[i]);

  MITK_INFO << "Features of " << missingKeys.size() << " newly labeled streamlines computed, " << keys.size() - missingKeys.size() << " reused from cache";
}

void StreamlineFeatureExtractor::PruneLocalColumns(const std::vector<std::uint64_t> &keys)
{
  const std::unordered_set<std::uint64_t> used(keys.begin(), keys.end());
  for (auto it = m_LocalColumns.begin(); it != m_LocalColumns.end();)
  {
    if (used.count(it->first))
      ++it;
    else
      it = m_LocalColumns.erase(it);
  }
}

void StreamlineFeatureExtractor::ReducePool(const std::vector<std::uint64_t> &positiveKeys)
{
  const std::size_t numCandidates = m_PoolIds.size();

  // distance of each candidate to the closest positive label
  std::vector<float> score(numCandidates, std::numeric_limits<float>::max());
  for (auto key : positiveKeys)
  {
    const auto& mdf = m_LocalColumns.at(key).mdf;
    for (std::size_t i = 0; i < numCandidates; ++i)
      score[i] = std::min(score[i], mdf[i]);
  }

  const std::size_t numKeep = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(PoolReductionFraction * numCandidates)));
  std::vector<float> sorted = score;
  std::nth_element(sorted.begin(), sorted.begin() + numKeep - 1, sorted.end());
  const float threshold = sorted[numKeep - 1];

  std::vector<std::size_t> keep;
  for (std::size_t i = 0; i < numCandidates; ++i)
  {
    if (m_PoolActive[i] && score[i] <= threshold)
      keep.push_back(i);
  }

  auto compact = [&keep](std::vector<float>& values) {
    for (std::size_t k = 0; k < keep.size(); ++k)
      values[k] = values[keep[k]];
    values.resize(keep.size());
  };

  FiberMatrix fibers(keep.size() * FiberStride);
  for (std::size_t k = 0; k < keep.size(); ++k)
  {
    std::copy_n(m_PoolFibers.begin() + keep[k] * FiberStride, FiberStride, fibers.begin() + k * FiberStride);
    m_PoolIds[k] = m_PoolIds[keep[k]];
  }
  m_PoolFibers.swap(fibers);
  m_PoolIds.resize(keep.size());
  m_PoolActive.assign(keep.size(), 1);

  for (auto& column : m_PrototypeColumns)
  {
    compact(column.mdf);
    compact(column.endpoints);
  }
  for (auto& entry : m_LocalColumns)
  {
    compact(entry.second.mdf);
    compact(entry.second.endpoints);
  }

  MITK_INFO << "Pool reduced from " << numCandidates << " to " << keep.size() << " streamlines close to the positive labels";
}

cv::Mat StreamlineFeatureExtractor::ComputeFeatures(const FiberMatrix &fibers, const FiberMatrix &localReferences) const
{
  const int numFibers = static_cast<int>(fibers.size() / FiberStride);
  const int numPrototypes = static_cast<int>(m_Prototypes.size() / FiberStride);
  const int numReferences = numPrototypes + static_cast<int>(localReferences.size() / FiberStride);

  // layout: [mdf to all references | endpoint distance to all references]
  cv::Mat features(numFibers, 2 * numReferences, CV_32F);

#pragma omp parallel for
  for (int i = 0; i < numFibers; ++i)
  {
    const float* fiber = &fibers[i * FiberStride];
    float* row = features.ptr<float>(i);
    for (int r = 0; r < numReferences; ++r)
    {
      const float* reference = r < numPrototypes ? &m_Prototypes[r * FiberStride] : &localReferences[(r - numPrototypes) * FiberStride];
      Distances(fiber, reference, row[r], row[numReferences + r]);
    }
  }
  return features;
}

cv::Mat StreamlineFeatureExtractor::AssembleCandidateFeatures(const std::vector<std::uint64_t> &localKeys)
{
  std::vector<const FeatureColumn*> columns;
  for (const auto& column : m_PrototypeColumns)
    columns.push_back(&column);
  for (auto key : localKeys)
    columns.push_back(&m_LocalColumns.at(key));

  std::vector<std::size_t> positions;
  m_CandidateIds.clear();
  for (std::size_t i = 0; i < m_PoolIds.size(); ++i)
  {
    if (m_PoolActive[i])
    {
      positions.push_back(i);
      m_CandidateIds.push_back(m_PoolIds[i]);
    }
  }

  const int numReferences = static_cast<int>(columns.size());
  cv::Mat features(static_cast<int>(positions.size()), 2 * numReferences, CV_32F);

#pragma omp parallel for
  for (int i = 0; i < features.rows; ++i)
  {
    float* row = features.ptr<float>(i);
    for (int r = 0; r < numReferences; ++r)
    {
      row[r] = columns[r]->mdf[positions[i]];
      row[numReferences + r] = columns[r]->endpoints[positions[i]];
    }
  }
  return features;
}

void StreamlineFeatureExtractor::Update()
{
  if (m_PoolIds.empty())
    throw std::runtime_error("The tractogram to dissect is empty.");
  if (m_Prototypes.empty())
    throw std::runtime_error("No prototypes are set.");

  const FiberMatrix plus = ToMatrix(m_TractogramPlus, "The positive bundle");
  const FiberMatrix minus = ToMatrix(m_TractogramMinus, "The negative bundle");
  if (plus.empty() || minus.empty())
    throw std::runtime_error("Label at least one positive and one negative streamline before training the classifier.");

  // local references: negative labels followed by positive labels
  FiberMatrix localReferences = minus;
  localReferences.insert(localReferences.end(), plus.begin(), plus.end());

  std::vector<std::uint64_t> plusKeys;
  for (std::size_t i = 0; i < plus.size(); i += FiberStride)
    plusKeys.push_back(Hash(&plus[i]));
  std::vector<std::uint64_t> localKeys;
  for (std::size_t i = 0; i < localReferences.size(); i += FiberStride)
    localKeys.push_back(Hash(&localReferences[i]));

  PruneLocalColumns(localKeys);

  // Large pools are reduced in the first cycle before the prototype features are computed, which saves most of the work.
  if (m_ActiveCycle == 0 && m_PoolIds.size() > PoolReductionThreshold)
  {
    AddLocalColumns(plusKeys, plus);
    ReducePool(plusKeys);
  }

  AddLocalColumns(localKeys, localReferences);
  if (m_PrototypeColumns.empty())
  {
    MITK_INFO << "Computing prototype features of " << m_PoolIds.size() << " streamlines (done once per session)";
    m_PrototypeColumns = ComputeColumns(m_Prototypes);
  }

  MITK_INFO << "Calculate features of training data";
  TrainModel(ComputeFeatures(plus, localReferences), ComputeFeatures(minus, localReferences));

  MITK_INFO << "Predict " << std::count(m_PoolActive.begin(), m_PoolActive.end(), 1) << " streamlines";
  Predict(AssembleCandidateFeatures(localKeys));
}

void StreamlineFeatureExtractor::TrainModel(const cv::Mat &plusFeatures, const cv::Mat &minusFeatures)
{
  cv::Mat samples;
  cv::Mat labels;
  samples.push_back(plusFeatures);
  labels.push_back(cv::Mat(plusFeatures.rows, 1, CV_32S, cv::Scalar(1)));
  samples.push_back(minusFeatures);
  labels.push_back(cv::Mat(minusFeatures.rows, 1, CV_32S, cv::Scalar(0)));

  // Balance the classes by repeating the samples of the smaller class
  const bool plusIsMinority = plusFeatures.rows < minusFeatures.rows;
  const cv::Mat& minority = plusIsMinority ? plusFeatures : minusFeatures;
  const int missing = std::abs(plusFeatures.rows - minusFeatures.rows);
  for (int i = 0; i < missing; ++i)
  {
    samples.push_back(minority.row(i % minority.rows));
    labels.push_back(plusIsMinority ? 1 : 0);
  }

  MITK_INFO << "Start training on " << plusFeatures.rows << " positive and " << minusFeatures.rows << " negative streamlines ("
            << missing << " " << (plusIsMinority ? "positive" : "negative") << " samples repeated for class balance)";

  m_Model = cv::ml::RTrees::create();
  m_Model->setMaxDepth(10);
  m_Model->setMinSampleCount(2);
  m_Model->setTruncatePrunedTree(false);
  m_Model->setUse1SERule(false);
  m_Model->setUseSurrogates(false);
  m_Model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 300, 0));
  m_Model->setCVFolds(1);
  m_Model->train(cv::ml::TrainData::create(samples, cv::ml::ROW_SAMPLE, labels));
}

void StreamlineFeatureExtractor::Predict(const cv::Mat &features)
{
  const int numCandidates = features.rows;
  m_CandidateFeatures = features;
  m_Entropies.assign(numCandidates, 0.0f);
  std::vector<char> positive(numCandidates, 0);

#pragma omp parallel for
  for (int i = 0; i < numCandidates; ++i)
  {
    // first row: class labels, second row: number of trees voting for the class
    cv::Mat votes;
    m_Model->getVotes(features.row(i), votes, 0);

    int positiveVotes = 0;
    int negativeVotes = 0;
    for (int c = 0; c < votes.cols; ++c)
    {
      if (votes.at<int>(0, c) == 1)
        positiveVotes += votes.at<int>(1, c);
      else
        negativeVotes += votes.at<int>(1, c);
    }

    positive[i] = positiveVotes > negativeVotes;
    const int totalVotes = positiveVotes + negativeVotes;
    if (totalVotes > 0)
      m_Entropies[i] = BinaryEntropy(static_cast<double>(positiveVotes) / totalVotes);
  }

  m_PredictedIds.clear();
  for (int i = 0; i < numCandidates; ++i)
  {
    if (positive[i])
      m_PredictedIds.push_back(m_CandidateIds[i]);
  }

  const auto numVeryUncertain = std::count_if(m_Entropies.begin(), m_Entropies.end(), [](float e) { return e >= 0.95f; });
  const std::size_t numUncertain = std::min<std::size_t>(numVeryUncertain > 500 ? 500 : 100, numCandidates);

  std::vector<unsigned int> order(numCandidates);
  std::iota(order.begin(), order.end(), 0);
  std::partial_sort(order.begin(), order.begin() + numUncertain, order.end(),
                    [this](unsigned int a, unsigned int b) { return m_Entropies[a] > m_Entropies[b]; });

  m_UncertainIds.clear();
  for (std::size_t i = 0; i < numUncertain; ++i)
    m_UncertainIds.push_back(m_CandidateIds[order[i]]);

  MITK_INFO << "Number of streamlines in prediction: " << m_PredictedIds.size();
}

std::vector<unsigned int> StreamlineFeatureExtractor::GetDiverseUncertainIds(float entropyThreshold, unsigned int count) const
{
  // candidates still in the pool, most uncertain first
  std::vector<unsigned int> order;
  for (std::size_t i = 0; i < m_CandidateIds.size(); ++i)
  {
    if (IsInPool(m_CandidateIds[i]))
      order.push_back(i);
  }
  std::stable_sort(order.begin(), order.end(), [this](unsigned int a, unsigned int b) { return m_Entropies[a] > m_Entropies[b]; });

  std::size_t numCandidates = std::count_if(order.begin(), order.end(), [&](unsigned int i) { return m_Entropies[i] >= entropyThreshold; });
  numCandidates = std::min(numCandidates, MaxDiverseCandidates);
  numCandidates = std::min(std::max<std::size_t>(numCandidates, count), order.size());
  order.resize(numCandidates);

  // spread the selection in feature space, starting with the most uncertain streamline
  const auto selected = FarthestFirst(numCandidates, count, 0, [&](unsigned int a, unsigned int b) {
    return cv::norm(m_CandidateFeatures.row(order[a]), m_CandidateFeatures.row(order[b]), cv::NORM_L1) / m_CandidateFeatures.cols;
  });

  std::vector<unsigned int> ids;
  for (auto i : selected)
    ids.push_back(m_CandidateIds[order[i]]);
  return ids;
}

bool StreamlineFeatureExtractor::IsInPool(unsigned int id) const
{
  const auto it = std::lower_bound(m_PoolIds.begin(), m_PoolIds.end(), id);
  return it != m_PoolIds.end() && *it == id && m_PoolActive[it - m_PoolIds.begin()];
}

void StreamlineFeatureExtractor::RemoveFromPool(const std::vector<unsigned int> &ids)
{
  for (auto id : ids)
  {
    const auto it = std::lower_bound(m_PoolIds.begin(), m_PoolIds.end(), id);
    if (it != m_PoolIds.end() && *it == id)
      m_PoolActive[it - m_PoolIds.begin()] = 0;
  }
}

mitk::FiberBundle::Pointer StreamlineFeatureExtractor::ExtractFibers(const std::vector<unsigned int> &ids) const
{
  return SubsetByIds(m_TractogramTest, ids);
}

mitk::FiberBundle::Pointer StreamlineFeatureExtractor::SubsetByIds(mitk::FiberBundle *fib, const std::vector<unsigned int> &ids)
{
  auto weights = vtkSmartPointer<vtkFloatArray>::New();
  auto polyData = fib->GeneratePolyDataByIds(ids, weights);
  auto subset = mitk::FiberBundle::New(polyData);
  subset->SetFiberWeights(weights);
  return subset;
}

std::vector<unsigned int> StreamlineFeatureExtractor::SelectDiverseFibers(mitk::FiberBundle *fib, unsigned int count, unsigned int subsetSize)
{
  std::vector<unsigned int> ids(fib->GetFiberPolyData()->GetNumberOfCells());
  std::iota(ids.begin(), ids.end(), 0);
  std::shuffle(ids.begin(), ids.end(), std::mt19937(std::random_device()()));
  ids.resize(std::min<std::size_t>(std::max(subsetSize, count), ids.size()));

  const FiberMatrix fibers = ToMatrix(fib, ids, "The tractogram");
  const auto selected = FarthestFirst(ids.size(), count, 0, [&](unsigned int a, unsigned int b) {
    float mdf, endpoints;
    Distances(&fibers[a * FiberStride], &fibers[b * FiberStride], mdf, endpoints);
    return mdf;
  });

  std::vector<unsigned int> result;
  for (auto i : selected)
    result.push_back(ids[i]);
  return result;
}

}
