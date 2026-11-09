/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkActiveLearningSession_h
#define mitkActiveLearningSession_h

#include "MitkFiberDissectionExports.h"

#include <mitkFiberBundle.h>
#include <mitkPoint.h>

#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <vector>

namespace mitk
{
  class StreamlineFeatureExtractor;

  /**
   * \brief An atTRACTive labeling session without GUI: the fibers of a tractogram, their labels
   *        and the classifier, for callers such as the REST endpoints of the Workbench agent.
   *
   * The session works on a pool: a shuffled copy of the tractogram with every fiber resampled to
   * StreamlineFeatureExtractor::NumPoints points, as the Interactive Fiber Dissection view
   * prepares it. All fiber ids the session takes and gives are pool ids, from 0 to
   * GetNumberOfPoolFibers() - 1; GetOriginalIds() maps them to the fibers of the tractogram.
   *
   * A cycle samples fibers (SampleInitial, then SampleUncertain or SampleDiverse after a
   * training), labels them (Label) and trains the classifier (Train), which predicts the tract:
   * the classified fibers it puts into the tract and the positive labels. GetFinishedTract() gives
   * the prediction as the tractogram's own fibers.
   *
   * With a seed, pool, samples and predictions are the same at every run.
   *
   * Not thread-safe: the caller serializes the calls of a session.
   */
  class MITKFIBERDISSECTION_EXPORT ActiveLearningSession
  {
  public:
    /** Entropy from which a classified fiber counts as uncertain (Metrics::Uncertain). */
    static constexpr float UncertainEntropy = 0.95f;

    struct Metrics
    {
      /** The cycle of the training, from 0. */
      unsigned int Cycle = 0;
      std::size_t Positives = 0;
      std::size_t Negatives = 0;
      /** Fibers of the prediction, the positive labels included. */
      std::size_t Predicted = 0;
      /** Fibers that entered or left the prediction since the last training, relative to the fibers
       *  of either prediction; 1 in cycle 0. */
      double ChangedFraction = 1;
      /** Classified fibers with an entropy of at least UncertainEntropy. */
      std::size_t Uncertain = 0;
    };

    /**
     * \brief Builds the pool and selects the prototypes (subset farthest first on a subset of
     *        3 n ln n fibers, as the view does).
     * \param seed Seeds pool order, samples and classifier; without it every session differs.
     * \throw mitk::Exception for a tractogram without fibers or no prototypes.
     */
    ActiveLearningSession(FiberBundle* tractogram, unsigned int numPrototypes = 75, std::optional<std::uint32_t> seed = {});
    ~ActiveLearningSession();

    ActiveLearningSession(const ActiveLearningSession&) = delete;
    ActiveLearningSession& operator=(const ActiveLearningSession&) = delete;

    FiberBundle* GetTractogram() const { return m_Tractogram; }
    /** The resampled fibers; pool id i is fiber i of this bundle. */
    FiberBundle* GetPool() const { return m_Pool; }
    unsigned int GetNumberOfPoolFibers() const { return m_Pool->GetNumFibers(); }
    unsigned int GetNumberOfPrototypes() const { return m_NumPrototypes; }

    /** The points of a pool fiber. \throw mitk::Exception for an id out of range. */
    std::vector<Point3D> GetPoolFiberPoints(unsigned int poolId) const;
    /** The tractogram's fiber ids of pool fibers. \throw mitk::Exception for an id out of range. */
    std::vector<unsigned int> GetOriginalIds(const std::vector<unsigned int>& poolIds) const;
    /** The pool ids of fibers whose points equal those of a pool fiber, like the fibers of a
     *  bundle extracted from GetPool(); other fibers are left out. Sorted. */
    std::vector<unsigned int> GetPoolIdsOf(const FiberBundle* fibers) const;

    /**
     * \brief Fibers to label before the first training: half of count drawn at random from
     *        enrichment (pool ids, e.g. the fibers that end in two spheres at the tract's ends), the
     *        rest spread over the pool (subset farthest first on max(500, 3 n ln n) fibers).
     *
     * Only fibers that are neither labeled nor handed out since the last training are drawn.
     */
    std::vector<unsigned int> SampleInitial(unsigned int count, const std::vector<unsigned int>& enrichment = {});
    /** The most uncertain classified fibers not handed out yet. \throw mitk::Exception before a training. */
    std::vector<unsigned int> SampleUncertain(unsigned int count);
    /** Uncertain fibers (entropy at least entropyThreshold) far apart from each other in feature
     *  space. \throw mitk::Exception before a training. */
    std::vector<unsigned int> SampleDiverse(float entropyThreshold, unsigned int count);

    /**
     * \brief Labels fibers as part of the tract (positive) or not (negative), or removes their
     *        labels. A labeled fiber may be labeled again; a predicted fiber labeled negative
     *        leaves the tract at the next training.
     * \throw mitk::Exception for an id out of range or in two of the lists; then no label changes.
     */
    void Label(const std::vector<unsigned int>& positive, const std::vector<unsigned int>& negative, const std::vector<unsigned int>& unlabel = {});
    const std::set<unsigned int>& GetPositiveIds() const { return m_Positive; }
    const std::set<unsigned int>& GetNegativeIds() const { return m_Negative; }

    /**
     * \brief Trains the classifier on the labels and predicts the tract.
     * \throw mitk::Exception without at least one positive and one negative label.
     */
    Metrics Train();
    /** Trainings done so far: the cycle of the next training. */
    unsigned int GetNumberOfTrainings() const { return m_Cycle; }
    /** The pool ids of the last prediction, sorted. */
    const std::vector<unsigned int>& GetPrediction() const { return m_Prediction; }

    /** The fibers of a pool subset, resampled, with their weights. */
    FiberBundle::Pointer GetPoolSubset(const std::vector<unsigned int>& poolIds) const;
    /** The predicted fibers as the tractogram holds them, in its order, with their weights. */
    FiberBundle::Pointer GetFinishedTract() const;

  private:
    void CheckIds(const std::vector<unsigned int>& ids) const;
    void RequireTraining() const;
    void HandOut(const std::vector<unsigned int>& ids);

    FiberBundle::Pointer m_Tractogram;
    FiberBundle::Pointer m_Pool;
    std::vector<unsigned int> m_PoolToOriginal;
    std::optional<std::uint32_t> m_Seed;
    std::mt19937 m_Engine;
    unsigned int m_NumPrototypes = 0;
    std::unique_ptr<StreamlineFeatureExtractor> m_Classifier;

    std::set<unsigned int> m_Positive;
    std::set<unsigned int> m_Negative;
    unsigned int m_Cycle = 0;
    std::vector<unsigned int> m_Prediction;
    /** Position in the classifier's uncertain fibers up to which SampleUncertain handed them out. */
    std::size_t m_UncertainPosition = 0;
  };
}

#endif
