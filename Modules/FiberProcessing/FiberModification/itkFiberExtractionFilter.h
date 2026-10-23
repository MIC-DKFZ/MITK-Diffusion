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

namespace itk{

/**
* \brief Extract streamlines from tractograms using ROI images */

template< class PixelType >
class FiberExtractionFilter : public ProcessObject
{
public:

  enum class MODE {
    OVERLAP,
    ENDPOINTS
  };

  enum class INPUT {
    SCALAR_MAP,  ///< In this case, positive means roi image vlaue > threshold
    LABEL_MAP    ///< In this case, positive means roi image value in labels vector
  };

  typedef FiberExtractionFilter Self;
  typedef ProcessObject                   Superclass;
  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;
  typedef itk::Image< PixelType , 3>      ItkInputImgType;


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
  itkSetMacro( BothEnds, bool )             ///< Both streamline ends need to be inside of the ROI for the streamline to be considered positive
  itkSetMacro( OverlapFraction, float )     ///< Necessary fraction of streamline points inside the ROI for the streamline to be considered positive
  itkSetMacro( NoNegatives, bool )          ///< Don't create output tractograms from negative streamlines (save the computation)
  itkSetMacro( NoPositives, bool )          ///< Don't create output tractograms from positive streamlines (save the computation)
  itkSetMacro( Interpolate, bool )          ///< Interpolate input ROI image (only relevant for ENDPOINTS mode
  itkSetMacro( Threshold, float )           ///< Threshold on input ROI image value to determine positives/negatives
  itkSetMacro( SkipSelfConnections, bool )  ///< Ignore streamlines between two identical labels
  itkSetMacro( OnlySelfConnections, bool )  ///< Only keep streamlines between two identical labels
  itkSetMacro( SplitLabels, bool )          ///< Output a separate tractogram for each label-->label tract
  itkSetMacro( SplitByRoi, bool )           ///< Output a separate tractogram for each ROI image
  itkSetMacro( MinFibersPerTract, unsigned int )  ///< Discard positives with less fibers
  itkSetMacro( PairedStartEndLabels, bool )

  void SetRoiImages(const std::vector< ItkInputImgType* > &rois);
  void SetRoiImageNames(const std::vector< std::string > roi_names);
  void SetLabels(const std::vector<unsigned short> &Labels);
  void SetStartLabels(const std::vector<unsigned short> &Labels);
  void SetEndLabels(const std::vector<unsigned short> &Labels);

  std::vector<mitk::FiberBundle::Pointer> GetPositives() const; ///< Get positive tracts (filtered by the input ROIs)
  std::vector<mitk::FiberBundle::Pointer> GetNegatives() const; ///< Get negative tracts (not filtered by the ROIs)
  std::vector< std::string > GetPositiveLabels() const;  ///< In case of label extraction, this vector contains the labels corresponding to the positive tracts

protected:

  void GenerateData() override;

  FiberExtractionFilter();
  ~FiberExtractionFilter() override;

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
