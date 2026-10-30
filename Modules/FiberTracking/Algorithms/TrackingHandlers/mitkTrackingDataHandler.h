/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _TrackingDataHandler
#define _TrackingDataHandler

#include <MitkFiberTrackingExports.h>
#include <algorithm>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <cstdint>
#include <deque>
#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkPoint.h>
#include <mitkBaseData.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkStreamlineTractographyParameters.h>

namespace mitk
{
  /**
   * \brief  Abstract class for tracking handler. A tracking handler deals with determining the next progression
   * direction of a streamline fiber. There are different handlers for tensor images, peak images, ...
   * Tensors, ODF directions and peaks are world (LPS) vectors (the project's gradient convention), so
   * the handlers use them without applying the image rotation.
   *
   * ProposeDirection does the steps all handlers share (voxel lookup, the same-voxel shortcut
   * without interpolation, the flip options) and asks the derived handler only for the direction
   * in the data frame (ProposeDataDirection). */

  class MITKFIBERTRACKING_EXPORT TrackingDataHandler
  {
  public:
    TrackingDataHandler();
    virtual ~TrackingDataHandler() {}

    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator ItkRngType;
    typedef boost::mt19937 BoostRngType;
    typedef itk::Image<unsigned char, 3> ItkUcharImgType;
    typedef itk::Image<short, 3> ItkShortImgType;
    typedef itk::Image<float, 3> ItkFloatImgType;
    typedef itk::Image<double, 3> ItkDoubleImgType;
    typedef vnl_vector_fixed<float, 3> TrackingDirectionType;
    typedef mitk::StreamlineTractographyParameters::MODE MODE;

    /** Next progression direction at pos in the tracking frame, (0,0,0) to stop. olddirs holds the
     *  previous directions (tracking frame), oldIndex the voxel of the previous step. */
    TrackingDirectionType ProposeDirection(const itk::Point<float, 3> &pos,
                                           std::deque<TrackingDirectionType> &olddirs,
                                           itk::Index<3> &oldIndex);

    virtual void InitForTracking() = 0;
    virtual itk::Vector<double, 3> GetSpacing() = 0;
    virtual itk::Point<float, 3> GetOrigin() = 0;
    virtual itk::Matrix<double, 3, 3> GetDirection() = 0;
    virtual itk::ImageRegion<3> GetLargestPossibleRegion() = 0;
    virtual bool WorldToIndex(itk::Point<float, 3> &pos, itk::Index<3> &index) = 0;

    void SetParameters(std::shared_ptr<mitk::StreamlineTractographyParameters> parameters)
    {
      m_Parameters = parameters;

      if (m_Parameters->m_FixRandomSeed)
        m_RngItk->SetSeed(0);
      else
        m_RngItk->SetSeed();
    }

    /** Uniform random number from the handler's own engine, which is seeded in SetParameters.
     *  For set-up work outside the tracking threads (seed point jitter); not thread-safe. */
    double GetRandDouble(const double &a, const double &b) { return m_RngItk->GetUniformVariate(a, b); }

    /** Creates one random engine per tracking thread. The tracker calls it before the
     *  parallel region with the number of OpenMP threads it will use. */
    void InitThreadRandomEngines(int numThreads);

    /** Seeds the random engines of the calling OpenMP thread. The tracker calls it at each seed
     *  point with a value derived from the seed point's index, so that a streamline depends only
     *  on its seed point and not on the number of threads or their scheduling. */
    void SeedThreadRandomEngines(std::uint32_t seed);

    /** Uniform random number from the calling thread's engine (see SeedThreadRandomEngines). */
    double GetThreadRandDouble(double a, double b) { return GetThreadRngItk()->GetUniformVariate(a, b); }

  protected:
    /** False where the handler's stopping value (FA, GFA) at pos is below the cutoff. Checked before
     *  the same-voxel shortcut of ProposeDirection. */
    virtual bool PassesCutoff(const itk::Point<float, 3> & /*pos*/) { return true; }

    /** Direction at pos in the data frame, (0,0,0) to stop, including the handler's angular
     *  threshold test. oldDir is the previous direction in the data frame, (0,0,0) at the seed. */
    virtual TrackingDirectionType ProposeDataDirection(const itk::Point<float, 3> &pos,
                                                       const TrackingDirectionType &oldDir) = 0;

    /**
     * Trilinear interpolation of per-voxel directions at the continuous index cIdx, with the corner
     * indices clamped to region, as ITK's interpolators do for the FA/GFA, ODF and mask images, so
     * that the outermost half voxel is tracked too. matchAt(corner, weight) gives the corner's
     * direction, already matched in sign to the previous direction; the result is the weighted sum.
     */
    template <typename MatchAt>
    static TrackingDirectionType InterpolateDirection(const itk::ContinuousIndex<float, 3> &cIdx,
                                                      const itk::ImageRegion<3> &region,
                                                      MatchAt &&matchAt)
    {
      itk::Index<3> lo, hi;
      float frac[3];
      for (int d = 0; d < 3; ++d)
      {
        auto size = static_cast<itk::IndexValueType>(region.GetSize(d));
        auto base = static_cast<itk::IndexValueType>(std::floor(cIdx[d]));
        frac[d] = 1 - (cIdx[d] - base); // weight of the lower corner
        lo[d] = std::clamp<itk::IndexValueType>(base, 0, size - 1);
        hi[d] = std::clamp<itk::IndexValueType>(base + 1, 0, size - 1);
      }
      const float frac_x = frac[0], frac_y = frac[1], frac_z = frac[2];

      vnl_vector_fixed<float, 8> interpWeights;
      interpWeights[0] = (  frac_x)*(  frac_y)*(  frac_z);
      interpWeights[1] = (1-frac_x)*(  frac_y)*(  frac_z);
      interpWeights[2] = (  frac_x)*(1-frac_y)*(  frac_z);
      interpWeights[3] = (  frac_x)*(  frac_y)*(1-frac_z);
      interpWeights[4] = (1-frac_x)*(1-frac_y)*(  frac_z);
      interpWeights[5] = (  frac_x)*(1-frac_y)*(1-frac_z);
      interpWeights[6] = (1-frac_x)*(  frac_y)*(1-frac_z);
      interpWeights[7] = (1-frac_x)*(1-frac_y)*(1-frac_z);

      // corner k uses the upper index on the axes marked true; same order as the weights
      const bool corners[8][3] = {{false, false, false}, {true, false, false}, {false, true, false}, {false, false, true},
                                  {true, true, false}, {false, true, true}, {true, false, true}, {true, true, true}};
      TrackingDirectionType dir;
      for (int k = 0; k < 8; ++k)
      {
        itk::Index<3> c = {{corners[k][0] ? hi[0] : lo[0], corners[k][1] ? hi[1] : lo[1], corners[k][2] ? hi[2] : lo[2]}};
        if (k == 0)
          dir = matchAt(c, interpWeights[k]) * interpWeights[k];
        else
          dir += matchAt(c, interpWeights[k]) * interpWeights[k];
      }
      return dir;
    }

    void CalculateMinVoxelSize()
    {
      itk::Vector<double, 3> imageSpacing = this->GetSpacing();
      float minVoxelSize = 0;
      if (imageSpacing[0] < imageSpacing[1] && imageSpacing[0] < imageSpacing[2])
        minVoxelSize = static_cast<float>(imageSpacing[0]);
      else if (imageSpacing[1] < imageSpacing[2])
        minVoxelSize = static_cast<float>(imageSpacing[1]);
      else
        minVoxelSize = static_cast<float>(imageSpacing[2]);

      if (m_Parameters == nullptr)
        mitkThrow() << "No tractography parameter opbect set!";

      m_Parameters->SetMinVoxelSizeMm(minVoxelSize);
    }

    /** Applies the flip options of the parameters. A flip is its own inverse, so the same call maps
     *  directions from the tracking frame to the data frame and back. */
    void FlipDirection(vnl_vector_fixed<float, 3> &dir) const
    {
      if (m_Parameters->m_FlipX)
        dir[0] *= -1;
      if (m_Parameters->m_FlipY)
        dir[1] *= -1;
      if (m_Parameters->m_FlipZ)
        dir[2] *= -1;
    }

    /** Engines of the calling OpenMP thread. */
    ItkRngType *GetThreadRngItk();
    BoostRngType &GetThreadRng();

    ItkRngType::Pointer m_RngItk;
    std::vector<ItkRngType::Pointer> m_ThreadRngItk;
    std::vector<BoostRngType> m_ThreadRng;
    bool m_NeedsDataInit;
    std::shared_ptr<mitk::StreamlineTractographyParameters> m_Parameters;

    void DataModified() { m_NeedsDataInit = true; }
  };

} // namespace mitk

#endif
