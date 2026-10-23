/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkFiberExtractionFilter_h
#define itkFiberExtractionFilter_h

// MITK
#include <mitkFiberBundle.h>

// ITK
#include <itkProcessObject.h>
#include <itkLinearInterpolateImageFunction.h>

#include <string>
#include <vector>

namespace itk{

/**
* \brief Extracts streamlines from a tractogram with ROI images. OVERLAP mode: a fiber is positive for
* the ROI that holds the largest fraction of its length (at least OverlapFraction) in voxels at or
* above Threshold; it ignores Interpolate and label maps. ENDPOINTS mode: a fiber is positive if
* both (BothEnds) or one of its end points lie in the ROI (SCALAR_MAP, value at or above Threshold)
* or in one of the labels (LABEL_MAP). Each Update replaces the outputs.
*/
template< class PixelType >
class FiberExtractionFilter : public ProcessObject
{
public:

  enum class MODE {
    OVERLAP,
    ENDPOINTS
  };

  enum class INPUT {
    /** positive means ROI image value >= threshold */
    SCALAR_MAP,
    /** positive means ROI image value in the labels */
    LABEL_MAP
  };

  using Self = FiberExtractionFilter;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;
  using ItkInputImgType = itk::Image< PixelType , 3>;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( FiberExtractionFilter, ProcessObject )

  void Update() override{
    this->GenerateData();
  }

  itkSetConstObjectMacro( InputFiberBundle, mitk::FiberBundle )
  /** Overlap or endpoints. The setters are written out because itkSetMacro streams the value,
   *  which an enum class does not support. */
  void SetMode(MODE mode)
  {
    if (m_Mode != mode)
    {
      m_Mode = mode;
      this->Modified();
    }
  }
  /** Scalar map or label image */
  void SetInputType(INPUT inputType)
  {
    if (m_InputType != inputType)
    {
      m_InputType = inputType;
      this->Modified();
    }
  }
  /** Both streamline ends need to be inside the ROI for the streamline to be positive. */
  itkSetMacro( BothEnds, bool )
  /** Fraction of the streamline length inside the ROI needed for the streamline to be positive */
  itkSetMacro( OverlapFraction, float )
  /** No output tractogram of the negative streamlines */
  itkSetMacro( NoNegatives, bool )
  /** No output tractograms of the positive streamlines */
  itkSetMacro( NoPositives, bool )
  /** Interpolate the ROI image (ENDPOINTS mode only) */
  itkSetMacro( Interpolate, bool )
  /** Threshold on the ROI image value for positives (SCALAR_MAP) */
  itkSetMacro( Threshold, float )
  /** Ignore streamlines between two identical labels. */
  itkSetMacro( SkipSelfConnections, bool )
  /** Only keep streamlines between two identical labels. */
  itkSetMacro( OnlySelfConnections, bool )
  /** A separate tractogram for each label-to-label tract */
  itkSetMacro( SplitLabels, bool )
  /** A separate tractogram for each ROI image */
  itkSetMacro( SplitByRoi, bool )
  /** Positive tracts with fewer fibers are discarded. */
  itkSetMacro( MinFibersPerTract, unsigned int )
  /** Start and end labels are pairs: label i of the start labels belongs to label i of the end labels. */
  itkSetMacro( PairedStartEndLabels, bool )

  /** The filter keeps raw pointers: the images must live until Update has run. */
  void SetRoiImages(const std::vector< ItkInputImgType* > &rois);
  void SetRoiImageNames(const std::vector< std::string > roi_names);
  void SetLabels(const std::vector<unsigned short> &Labels);
  void SetStartLabels(const std::vector<unsigned short> &Labels);
  void SetEndLabels(const std::vector<unsigned short> &Labels);

  /** Positive tracts (selected by the ROIs) */
  std::vector<mitk::FiberBundle::Pointer> GetPositives() const;
  /** Negative tracts (not selected by any ROI) */
  std::vector<mitk::FiberBundle::Pointer> GetNegatives() const;
  /** For label extraction and split ROIs, the name of each positive tract */
  std::vector< std::string > GetPositiveLabels() const;

protected:

  void GenerateData() override;

  FiberExtractionFilter();
  ~FiberExtractionFilter() override = default;

  mitk::FiberBundle::Pointer CreateFib(const std::vector< unsigned int >& ids);
  void AddPositives(const std::vector< std::vector< unsigned int > >& positive_ids);
  void ExtractOverlap(const mitk::FiberBundle* fib);
  void ExtractEndpoints(const mitk::FiberBundle* fib);
  void ExtractLabels(const mitk::FiberBundle* fib);
  bool IsPositive(const itk::Point<float, 3>& itkP);

  mitk::FiberBundle::ConstPointer             m_InputFiberBundle;
  std::vector< mitk::FiberBundle::Pointer >   m_Positives;
  std::vector< mitk::FiberBundle::Pointer >   m_Negatives;
  std::vector< ItkInputImgType* >             m_RoiImages;
  std::vector< std::string >                  m_RoiImageNames;
  MODE                                        m_Mode;
  INPUT                                       m_InputType;
  bool                                        m_BothEnds;
  float                                       m_OverlapFraction;
  bool                                        m_NoNegatives;
  bool                                        m_NoPositives;
  bool                                        m_Interpolate;
  float                                       m_Threshold;
  std::vector< unsigned short >               m_Labels;
  bool                                        m_SkipSelfConnections;
  bool                                        m_OnlySelfConnections;
  bool                                        m_SplitByRoi;
  bool                                        m_SplitLabels;
  unsigned int                                m_MinFibersPerTract;
  std::vector< unsigned short >               m_StartLabels;
  std::vector< unsigned short >               m_EndLabels;
  bool                                        m_PairedStartEndLabels;
  std::vector< std::string >                  m_PositiveLabels;
  typename itk::LinearInterpolateImageFunction< itk::Image< PixelType, 3 >, float >::Pointer   m_Interpolator;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberExtractionFilter.cpp"
#endif

#endif
