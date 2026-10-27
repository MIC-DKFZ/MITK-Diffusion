/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _TrackingHandlerOdf
#define _TrackingHandlerOdf

#include "mitkTrackingDataHandler.h"
#include <mitkOdfImage.h>
#include <mitkTensorImage.h>
#include <itkOrientationDistributionFunction.h>
#include <MitkFiberTrackingExports.h>

namespace mitk
{

/**
* \brief Enables streamline tracking on sampled ODF images. */

class MITKFIBERTRACKING_EXPORT TrackingHandlerOdf : public TrackingDataHandler
{

public:

  TrackingHandlerOdf();
  ~TrackingHandlerOdf() override;

  typedef OdfImage::ItkOdfImageType ItkOdfImageType;
  typedef itk::Image< vnl_vector_fixed<float,3>, 3>  ItkPDImgType;


  void InitForTracking() override;     ///< calls InputDataValidForTracking() and creates feature images
  vnl_vector_fixed<float,3> ProposeDirection(const itk::Point<float, 3>& pos, std::deque< vnl_vector_fixed<float,3> >& olddirs, itk::Index<3>& oldIndex) override;  ///< predicts next progression direction at the given position
  bool WorldToIndex(itk::Point<float, 3>& pos, itk::Index<3>& index) override;

  void SetOdfImage( ItkOdfImageType::Pointer img ){ m_OdfImage = img; DataModified(); }
  void SetGfaImage( ItkFloatImgType::Pointer img ){ m_GfaImage = img; DataModified(); }

  ItkUcharImgType::SpacingType GetSpacing() override{ return m_OdfImage->GetSpacing(); }
  itk::Point<float,3> GetOrigin() override{ return m_OdfImage->GetOrigin(); }
  ItkUcharImgType::DirectionType GetDirection() override{ return m_OdfImage->GetDirection(); }
  ItkUcharImgType::RegionType GetLargestPossibleRegion() override{ return m_OdfImage->GetLargestPossibleRegion(); }

  int OdfPower() const;

  bool GetIsOdfFromTensor() const;
  void SetIsOdfFromTensor(bool OdfFromTensor);

protected:

  int SampleOdf(vnl_vector< float >& probs, vnl_vector< float >& angles);

  /** Peak direction near the hemisphere sampling direction idx (same orientation as that direction):
   *  climbs to the local maximum of the sampled values and fits a quadratic to their logarithm in
   *  the tangent plane of the neighbouring directions, so that deterministic tracking is not
   *  limited to the sampling directions. With a previous direction, idx is the maximum of the ODF
   *  weighted with |cos| to that direction, but the climb uses the unweighted values: the result is
   *  the ODF peak whose basin holds idx (as MRtrix's SD_Stream searches the peak from the previous
   *  direction), at most 10 neighbour steps away. */
  vnl_vector_fixed<float,3> RefinePeak(int idx, const vnl_vector< float >& values) const;

  struct OdfNeighbor
  {
    int index;    ///< hemisphere index
    float sign;   ///< +1 or -1: orientation of the hemisphere direction relative to the centre
  };
  std::vector< std::vector< OdfNeighbor > > m_OdfNeighbors;  ///< hemisphere directions within 25 degrees

  ItkFloatImgType::Pointer        m_GfaImage;     ///< GFA image used to determine streamline termination.
  ItkOdfImageType::Pointer        m_OdfImage;     ///< Input odf image.
  std::vector< int >              m_OdfHemisphereIndices;
  vnl_matrix< float >             m_OdfFloatDirs;
  int                             m_NumProbSamples;     ///< samples drawn per step, the most probable one is kept
  bool                            m_OdfFromTensor;

  itk::LinearInterpolateImageFunction< itk::Image< float, 3 >, float >::Pointer   m_GfaInterpolator;
  itk::LinearInterpolateImageFunction< itk::Image< ItkOdfImageType::PixelType, 3 >, float >::Pointer   m_OdfInterpolator;
};

}

#endif
