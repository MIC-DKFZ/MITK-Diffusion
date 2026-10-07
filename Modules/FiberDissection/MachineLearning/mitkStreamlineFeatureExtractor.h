/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef StreamlineFeatureExtractor_h
#define StreamlineFeatureExtractor_h

#include "MitkFiberDissectionExports.h"

// MITK
#include <mitkFiberBundle.h>

// OpenCV
#include <opencv2/ml.hpp>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace mitk{

/**
* \brief Active-learning streamline classifier (atTRACTive).
*
* Every streamline is described by its distances to a set of reference streamlines: the global
* prototypes (fixed for a session) and the labeled streamlines (local prototypes). For each
* reference two features are used: the mean direct-flip distance (MDF) and the mean endpoint distance.
*
* The object is meant to live for a whole labeling session. Features of the test pool are cached per
* reference streamline, so a new training cycle only computes distances to newly labeled streamlines
* instead of recomputing the full feature matrix.
*
* All fiber ids taken or returned by this class refer to the fibers of the tractogram passed to SetTractogramTest().
*/
class MITKFIBERDISSECTION_EXPORT StreamlineFeatureExtractor
{
public:

  /** All streamlines (pool, prototypes, labels) must be resampled to this number of points. */
  static constexpr unsigned int NumPoints = 40;

  StreamlineFeatureExtractor();
  ~StreamlineFeatureExtractor();

  /** Sets the pool of streamlines to classify. Resets the pool state and all cached features. */
  void SetTractogramTest(const mitk::FiberBundle::Pointer &tractogram);
  /** Sets the global prototypes. Resets the cached prototype features. */
  void SetTractogramPrototypes(const mitk::FiberBundle::Pointer &prototypes);
  void SetTractogramPlus(const mitk::FiberBundle::Pointer &tractogram);
  void SetTractogramMinus(const mitk::FiberBundle::Pointer &tractogram);
  void SetActiveCycle(int activeCycle);

  /** Updates the features, trains the classifier and predicts the pool. Throws std::runtime_error on invalid input. */
  void Update();

  /** Pool fibers predicted as part of the target tract. */
  const std::vector<unsigned int>& GetPredictedIds() const { return m_PredictedIds; }
  /** Pool fibers sorted by decreasing prediction uncertainty. */
  const std::vector<unsigned int>& GetUncertainIds() const { return m_UncertainIds; }
  /** Pool fibers classified in the last Update(); aligned with GetEntropies(). */
  const std::vector<unsigned int>& GetCandidateIds() const { return m_CandidateIds; }
  const std::vector<float>& GetEntropies() const { return m_Entropies; }

  /**
  * Picks up to 'count' uncertain fibers (entropy >= entropyThreshold) that are far apart from each other in
  * feature space. Fibers removed from the pool since the last Update() are skipped.
  */
  std::vector<unsigned int> GetDiverseUncertainIds(float entropyThreshold, unsigned int count) const;

  /** Excludes fibers from future predictions, e.g. because they were handed to the user for labeling. */
  void RemoveFromPool(const std::vector<unsigned int> &ids);
  bool IsInPool(unsigned int id) const;

  /** Creates a new bundle containing the given pool fibers. */
  mitk::FiberBundle::Pointer ExtractFibers(const std::vector<unsigned int> &ids) const;

  /** Copies the given fibers of a bundle (with weights) into a new bundle. */
  static mitk::FiberBundle::Pointer SubsetByIds(mitk::FiberBundle *fib, const std::vector<unsigned int> &ids);

  /**
  * Selects 'count' fibers that cover the bundle well: a random subset of 'subsetSize' fibers is drawn and
  * farthest-first traversal on the MDF distance picks the fibers from it. All fibers need NumPoints points.
  */
  static std::vector<unsigned int> SelectDiverseFibers(mitk::FiberBundle *fib, unsigned int count, unsigned int subsetSize);

protected:

  /** Fibers as contiguous floats; fiber i occupies [i*FiberStride, (i+1)*FiberStride) as x,y,z per point. */
  typedef std::vector<float> FiberMatrix;
  static constexpr unsigned int FiberStride = 3 * NumPoints;

  /** Distances of all pool candidates to one reference streamline. */
  struct FeatureColumn
  {
    std::vector<float> mdf;
    std::vector<float> endpoints;
  };

  static FiberMatrix ToMatrix(mitk::FiberBundle *fib, const std::vector<unsigned int> &ids, const std::string &name);
  static FiberMatrix ToMatrix(mitk::FiberBundle *fib, const std::string &name);
  static std::uint64_t Hash(const float *fiber);
  static void Distances(const float *a, const float *b, float &mdf, float &endpoints);

  /** Computes the feature columns of the given reference fibers for all pool candidates. */
  std::vector<FeatureColumn> ComputeColumns(const FiberMatrix &references) const;
  /** Computes the columns of labeled fibers that are not cached yet. */
  void AddLocalColumns(const std::vector<std::uint64_t> &keys, const FiberMatrix &fibers);
  /** Drops cached columns of fibers that are no longer labeled. */
  void PruneLocalColumns(const std::vector<std::uint64_t> &keys);
  /** Restricts the pool to the candidates closest to the positive labels. */
  void ReducePool(const std::vector<std::uint64_t> &positiveKeys);
  /** Features of arbitrary fibers: global prototypes followed by the given local references. */
  cv::Mat ComputeFeatures(const FiberMatrix &fibers, const FiberMatrix &localReferences) const;
  /** Features of all active pool candidates from the cache, same layout as ComputeFeatures(). */
  cv::Mat AssembleCandidateFeatures(const std::vector<std::uint64_t> &localKeys);
  void TrainModel(const cv::Mat &plusFeatures, const cv::Mat &minusFeatures);
  void Predict(const cv::Mat &features);

  int                                         m_ActiveCycle;

  mitk::FiberBundle::Pointer                  m_TractogramTest;
  mitk::FiberBundle::Pointer                  m_TractogramPlus;
  mitk::FiberBundle::Pointer                  m_TractogramMinus;
  FiberMatrix                                 m_Prototypes;

  // pool candidates: sorted pool ids, their geometry and whether they are still available
  std::vector<unsigned int>                   m_PoolIds;
  std::vector<char>                           m_PoolActive;
  FiberMatrix                                 m_PoolFibers;

  // cached features of the pool candidates, one column per reference fiber
  std::vector<FeatureColumn>                  m_PrototypeColumns;
  std::unordered_map<std::uint64_t, FeatureColumn> m_LocalColumns;

  // results of the last Update()
  cv::Ptr<cv::ml::RTrees>                     m_Model;
  cv::Mat                                     m_CandidateFeatures;
  std::vector<unsigned int>                   m_CandidateIds;
  std::vector<float>                          m_Entropies;
  std::vector<unsigned int>                   m_PredictedIds;
  std::vector<unsigned int>                   m_UncertainIds;
};
}

#endif
