/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkActiveLearningSession.h"
#include "mitkStreamlineFeatureExtractor.h"

#include <mitkExceptionMacro.h>

#include <opencv2/core.hpp>

#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <numeric>
#include <stdexcept>

namespace
{
  std::vector<unsigned int> Sorted(const std::set<unsigned int>& ids)
  {
    return std::vector<unsigned int>(ids.begin(), ids.end());
  }

  /** Subset size of the subset farthest first selection of count fibers: 3 n ln n, as the view chooses it. */
  unsigned int SubsetSize(unsigned int count, unsigned int minimum)
  {
    return std::max(minimum, static_cast<unsigned int>(std::ceil(3.0 * count * std::log(std::max(count, 1u)))));
  }

  /** The coordinates of a fiber as the bundle stores them (float). */
  std::vector<float> Coordinates(vtkPolyData* polyData, unsigned int fiber, vtkIdList* pointIds)
  {
    polyData->GetCellPoints(fiber, pointIds);
    std::vector<float> coordinates(3 * pointIds->GetNumberOfIds());

    for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); ++i)
    {
      double p[3];
      polyData->GetPoint(pointIds->GetId(i), p);
      std::copy(p, p + 3, coordinates.begin() + 3 * i);
    }

    return coordinates;
  }

  /** A cheap key of a fiber's geometry: its point count and first, middle and last point. */
  std::vector<float> GeometryKey(const std::vector<float>& coordinates)
  {
    std::vector<float> key{ static_cast<float>(coordinates.size()) };

    if (coordinates.empty())
      return key;

    const std::size_t numPoints = coordinates.size() / 3;

    for (const auto point : { std::size_t(0), numPoints / 2, numPoints - 1 })
      key.insert(key.end(), coordinates.begin() + 3 * point, coordinates.begin() + 3 * point + 3);

    return key;
  }
}

mitk::ActiveLearningSession::ActiveLearningSession(FiberBundle* tractogram, unsigned int numPrototypes, std::optional<std::uint32_t> seed)
  : m_Tractogram(tractogram),
    m_Seed(seed),
    m_Engine(seed ? *seed : std::random_device{}())
{
  if (tractogram == nullptr || tractogram->GetNumFibers() == 0)
    mitkThrow() << "The tractogram has no fibers.";

  if (numPrototypes == 0)
    mitkThrow() << "A session needs at least one prototype.";

  // A shuffled copy, so that the order of the fibers in the file does not matter, as the view's
  // resampling makes it.
  m_PoolToOriginal.resize(tractogram->GetNumFibers());
  std::iota(m_PoolToOriginal.begin(), m_PoolToOriginal.end(), 0);
  std::shuffle(m_PoolToOriginal.begin(), m_PoolToOriginal.end(), m_Engine);
  m_Pool = StreamlineFeatureExtractor::SubsetByIds(tractogram, m_PoolToOriginal);
  m_Pool->ResampleToNumPoints(StreamlineFeatureExtractor::NumPoints);

  std::vector<unsigned int> poolIds(m_Pool->GetNumFibers());
  std::iota(poolIds.begin(), poolIds.end(), 0);
  const auto prototypeIds = StreamlineFeatureExtractor::SelectDiverseFibers(m_Pool, poolIds, numPrototypes, SubsetSize(numPrototypes, 1), m_Engine);
  m_NumPrototypes = static_cast<unsigned int>(prototypeIds.size());

  m_Classifier = std::make_unique<StreamlineFeatureExtractor>();
  m_Classifier->SetTractogramTest(m_Pool);
  m_Classifier->SetTractogramPrototypes(StreamlineFeatureExtractor::SubsetByIds(m_Pool, prototypeIds));
}

mitk::ActiveLearningSession::~ActiveLearningSession() = default;

void mitk::ActiveLearningSession::CheckIds(const std::vector<unsigned int>& ids) const
{
  for (const auto id : ids)
  {
    if (id >= m_Pool->GetNumFibers())
      mitkThrow() << "The session has no fiber " << id << "; its ids go from 0 to " << m_Pool->GetNumFibers() - 1 << ".";
  }
}

void mitk::ActiveLearningSession::RequireTraining() const
{
  if (m_Cycle == 0)
    mitkThrow() << "Train the classifier before sampling uncertain fibers.";
}

void mitk::ActiveLearningSession::HandOut(const std::vector<unsigned int>& ids)
{
  // The classifier keeps them from sampling until the next training, which classifies the fibers
  // left unlabeled again.
  m_Classifier->ExcludeFromSampling(ids);
}

std::vector<mitk::Point3D> mitk::ActiveLearningSession::GetPoolFiberPoints(unsigned int poolId) const
{
  this->CheckIds({ poolId });
  auto pointIds = vtkSmartPointer<vtkIdList>::New();
  const auto coordinates = Coordinates(m_Pool->GetFiberPolyData(), poolId, pointIds);
  std::vector<Point3D> points(coordinates.size() / 3);

  for (std::size_t i = 0; i < points.size(); ++i)
  {
    for (unsigned int d = 0; d < 3; ++d)
      points[i][d] = coordinates[3 * i + d];
  }

  return points;
}

std::vector<unsigned int> mitk::ActiveLearningSession::GetOriginalIds(const std::vector<unsigned int>& poolIds) const
{
  this->CheckIds(poolIds);
  std::vector<unsigned int> ids;
  ids.reserve(poolIds.size());

  for (const auto id : poolIds)
    ids.push_back(m_PoolToOriginal[id]);

  return ids;
}

std::vector<unsigned int> mitk::ActiveLearningSession::GetPoolIdsOf(const FiberBundle* fibers) const
{
  std::vector<unsigned int> result;

  if (fibers == nullptr || fibers->GetNumFibers() == 0)
    return result;

  auto pointIds = vtkSmartPointer<vtkIdList>::New();
  auto polyData = fibers->GetFiberPolyData();
  std::multimap<std::vector<float>, std::vector<float>> wanted;

  for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
  {
    auto coordinates = Coordinates(polyData, i, pointIds);
    auto key = GeometryKey(coordinates);
    wanted.emplace(std::move(key), std::move(coordinates));
  }

  auto pool = m_Pool->GetFiberPolyData();

  for (unsigned int id = 0; id < m_Pool->GetNumFibers(); ++id)
  {
    const auto coordinates = Coordinates(pool, id, pointIds);
    const auto range = wanted.equal_range(GeometryKey(coordinates));

    if (std::any_of(range.first, range.second, [&](const auto& entry) { return entry.second == coordinates; }))
      result.push_back(id);
  }

  return result;
}

std::vector<unsigned int> mitk::ActiveLearningSession::SampleInitial(unsigned int count, const std::vector<unsigned int>& enrichment)
{
  this->CheckIds(enrichment);

  std::vector<unsigned int> available;
  for (unsigned int id = 0; id < m_Pool->GetNumFibers(); ++id)
  {
    if (m_Classifier->IsAvailableForSampling(id))
      available.push_back(id);
  }

  // Half from the enrichment, drawn at random, as the view draws them from the fibers between its
  // two spheres.
  std::vector<unsigned int> enriched;
  for (const auto id : std::set<unsigned int>(enrichment.begin(), enrichment.end()))
  {
    if (m_Classifier->IsAvailableForSampling(id))
      enriched.push_back(id);
  }
  std::shuffle(enriched.begin(), enriched.end(), m_Engine);
  enriched.resize(std::min<std::size_t>(count / 2, enriched.size()));

  std::vector<unsigned int> ids = enriched;
  std::sort(enriched.begin(), enriched.end());
  std::vector<unsigned int> rest;
  std::set_difference(available.begin(), available.end(), enriched.begin(), enriched.end(), std::back_inserter(rest));

  // The rest covers the whole pool.
  const auto numDiverse = count - static_cast<unsigned int>(ids.size());
  const auto diverse = StreamlineFeatureExtractor::SelectDiverseFibers(m_Pool, std::move(rest), numDiverse, SubsetSize(numDiverse, 500), m_Engine);
  ids.insert(ids.end(), diverse.begin(), diverse.end());

  this->HandOut(ids);
  return ids;
}

std::vector<unsigned int> mitk::ActiveLearningSession::SampleUncertain(unsigned int count)
{
  this->RequireTraining();

  const auto& uncertain = m_Classifier->GetUncertainIds();
  std::vector<unsigned int> ids;

  while (m_UncertainPosition < uncertain.size() && ids.size() < count)
  {
    const auto id = uncertain[m_UncertainPosition++];

    if (m_Classifier->IsAvailableForSampling(id))
      ids.push_back(id);
  }

  this->HandOut(ids);
  return ids;
}

std::vector<unsigned int> mitk::ActiveLearningSession::SampleDiverse(float entropyThreshold, unsigned int count)
{
  this->RequireTraining();

  const auto ids = count > 0 ? m_Classifier->GetDiverseUncertainIds(entropyThreshold, count) : std::vector<unsigned int>();
  this->HandOut(ids);
  return ids;
}

void mitk::ActiveLearningSession::Label(const std::vector<unsigned int>& positive, const std::vector<unsigned int>& negative, const std::vector<unsigned int>& unlabel)
{
  this->CheckIds(positive);
  this->CheckIds(negative);
  this->CheckIds(unlabel);

  const std::set<unsigned int> positiveSet(positive.begin(), positive.end());
  const std::set<unsigned int> negativeSet(negative.begin(), negative.end());
  const std::set<unsigned int> unlabelSet(unlabel.begin(), unlabel.end());

  for (const auto id : positiveSet)
  {
    if (negativeSet.count(id) != 0 || unlabelSet.count(id) != 0)
      mitkThrow() << "Fiber " << id << " is in two of the lists positive, negative and unlabel.";
  }

  for (const auto id : negativeSet)
  {
    if (unlabelSet.count(id) != 0)
      mitkThrow() << "Fiber " << id << " is in two of the lists positive, negative and unlabel.";
  }

  for (const auto id : positiveSet)
  {
    m_Negative.erase(id);
    m_Positive.insert(id);
  }

  for (const auto id : negativeSet)
  {
    m_Positive.erase(id);
    m_Negative.insert(id);
  }

  for (const auto id : unlabelSet)
  {
    m_Positive.erase(id);
    m_Negative.erase(id);
  }

  // Labeled fibers are no candidates for sampling; the training makes them training samples.
  this->HandOut(positive);
  this->HandOut(negative);
}

mitk::ActiveLearningSession::Metrics mitk::ActiveLearningSession::Train()
{
  if (m_Positive.empty() || m_Negative.empty())
    mitkThrow() << "Label at least one fiber positive and one negative before training.";

  // OpenCV's random forest draws its bootstrap samples and variable subsets from the thread's
  // cv::theRNG(), whose state depends on what the thread trained before; a seeded session sets it.
  if (m_Seed)
    cv::theRNG() = cv::RNG((static_cast<std::uint64_t>(*m_Seed) << 32) + m_Cycle + 1);

  m_Classifier->SetActiveCycle(static_cast<int>(m_Cycle));
  m_Classifier->SetTractogramPlus(StreamlineFeatureExtractor::SubsetByIds(m_Pool, Sorted(m_Positive)));
  m_Classifier->SetTractogramMinus(StreamlineFeatureExtractor::SubsetByIds(m_Pool, Sorted(m_Negative)));

  try
  {
    m_Classifier->Update();
  }
  catch (const std::runtime_error& e)
  {
    mitkThrow() << e.what();
  }

  // Labeled fibers are not classified; the positive ones belong to the tract.
  std::vector<unsigned int> prediction = m_Classifier->GetPredictedIds();
  prediction.insert(prediction.end(), m_Positive.begin(), m_Positive.end());
  std::sort(prediction.begin(), prediction.end());
  prediction.erase(std::unique(prediction.begin(), prediction.end()), prediction.end());

  Metrics metrics;
  metrics.Cycle = m_Cycle;
  metrics.Positives = m_Positive.size();
  metrics.Negatives = m_Negative.size();
  metrics.Predicted = prediction.size();

  if (m_Cycle > 0)
  {
    std::vector<unsigned int> changed;
    std::set_symmetric_difference(m_Prediction.begin(), m_Prediction.end(), prediction.begin(), prediction.end(), std::back_inserter(changed));
    const std::size_t both = (m_Prediction.size() + prediction.size() + changed.size()) / 2;
    metrics.ChangedFraction = both > 0 ? static_cast<double>(changed.size()) / both : 0.0;
  }

  const auto& entropies = m_Classifier->GetEntropies();
  metrics.Uncertain = std::count_if(entropies.begin(), entropies.end(), [](float e) { return e >= UncertainEntropy; });

  m_Prediction.swap(prediction);
  m_UncertainPosition = 0;
  ++m_Cycle;
  return metrics;
}

mitk::FiberBundle::Pointer mitk::ActiveLearningSession::GetPoolSubset(const std::vector<unsigned int>& poolIds) const
{
  this->CheckIds(poolIds);
  return StreamlineFeatureExtractor::SubsetByIds(m_Pool, poolIds);
}

mitk::FiberBundle::Pointer mitk::ActiveLearningSession::GetFinishedTract() const
{
  auto ids = this->GetOriginalIds(m_Prediction);
  std::sort(ids.begin(), ids.end());
  auto tract = StreamlineFeatureExtractor::SubsetByIds(m_Tractogram, ids);

  // The metadata that FiberBundle's own operations keep (mitk::CopyFiberBundleMetadata of
  // FiberProcessing, which this module does not depend on): a set TrackVis header and the
  // coordinate convention, so that the tract is saved as the tractogram was.
  if (m_Tractogram->HasTrackVisHeader())
    tract->SetTrackVisHeader(m_Tractogram->GetTrackVisHeader());
  tract->setIsRAS(m_Tractogram->IsRAS());
  return tract;
}
