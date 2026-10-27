/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkStreamlineTrackingFilter_h
#define itkStreamlineTrackingFilter_h

#include <itkProcessObject.h>
#include <itkLinearInterpolateImageFunction.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <atomic>
#include <chrono>
#include <deque>
#include <TrackingHandlers/mitkTrackingDataHandler.h>
#include <MitkFiberTrackingExports.h>
#include <mitkFiberBundle.h>
#include <mitkStreamlineTractographyParameters.h>

namespace itk{

/**
* \brief Performs streamline tracking on the input image. Depending on the tracking handler this can be a tensor, ODF or peak based tracking. */

class MITKFIBERTRACKING_EXPORT StreamlineTrackingFilter : public ProcessObject
{

public:

  typedef StreamlineTrackingFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ProcessObject Superclass;

  /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(StreamlineTrackingFilter, ProcessObject)

  typedef itk::Image<unsigned char, 3>                ItkUcharImgType;
  typedef itk::Image<unsigned int, 3>                 ItkUintImgType;
  typedef itk::Image<double, 3>                       ItkDoubleImgType;
  typedef itk::Image<float, 3>                        ItkFloatImgType;
  typedef vtkSmartPointer< vtkPolyData >              PolyDataType;
  typedef mitk::StreamlineTractographyParameters::EndpointConstraints EndpointConstraints;
  typedef mitk::StreamlineTractographyParameters::MODE MODE;

  typedef std::deque< vnl_vector_fixed<float,3> > DirectionContainer;
  typedef std::deque< itk::Point<float> > FiberType;
  typedef std::vector< FiberType > BundleType;

  std::atomic<bool> m_AbortTracking;

  void SetDicomProperties(mitk::FiberBundle::Pointer fib);

  itkGetMacro( OutputProbabilityMap, ItkDoubleImgType::Pointer)    ///< Output probability map
  itkGetMacro( FiberPolyData, PolyDataType )            ///< Output fibers

  itkSetMacro( SeedImage, ItkFloatImgType::Pointer)     ///< Seeds are only placed inside of this mask.
  itkSetMacro( MaskImage, ItkFloatImgType::Pointer)     ///< Tracking is only performed inside of this mask image.
  itkSetMacro( ExclusionRegions, ItkFloatImgType::Pointer)///< Fibers passing any of the ROIs in this image are discarded.

  itkSetMacro( StoppingRegions, ItkFloatImgType::Pointer) ///< Streamlines entering a stopping region will stop immediately
  itkSetMacro( TargetRegions, ItkFloatImgType::Pointer)    ///< Only streamline starting and ending in this mask are retained
  itkSetMacro( Verbose, bool )                        ///< If true, output tracking progress (might be slower)

  /** Called from the GUI thread while the tracking threads run, hence no itkSetMacro (Modified() is not thread-safe). */
  void SetStopTracking(bool stop) { m_StopTracking = stop; }

  ///< Use manually defined points in physical space as seed points instead of seed image
  void SetSeedPoints( const std::vector< itk::Point<float> >& sP) {
    m_SeedPoints = sP;
  }

  void SetTrackingHandler( mitk::TrackingDataHandler* h )   ///<
  {
    m_TrackingHandler = h;
  }

  void Update() override{
    this->GenerateData();
  }

  std::string GetStatusText();

  void SetTrackingPriorHandler(mitk::TrackingDataHandler *TrackingPriorHandler);

  std::shared_ptr< mitk::StreamlineTractographyParameters > GetParameters() const;
  void SetParameters(std::shared_ptr< mitk::StreamlineTractographyParameters > Parameters);

protected:

  void GenerateData() override;

  StreamlineTrackingFilter();
  ~StreamlineTrackingFilter() override {}

  bool IsValidFiber(FiberType* fib);  ///< Check endpoints
  void FiberToProbmap(FiberType* fib);
  void GetSeedPointsFromSeedImage(const ItkFloatImgType* seedImage);
  void CalculateNewPosition(itk::Point<float, 3>& pos, vnl_vector_fixed<float,3>& dir);    ///< Calculate next integration step.
  float FollowStreamline(itk::Point<float, 3> start_pos, vnl_vector_fixed<float,3> dir, FiberType* fib, DirectionContainer* container, float tractLength, bool front, bool& exclude);       ///< Start streamline in one direction.
  vnl_vector_fixed<float,3> GetNewDirection(const itk::Point<float, 3>& pos, std::deque< vnl_vector_fixed<float,3> >& olddirs, itk::Index<3>& oldIndex); ///< Determine new direction by sample voting at the current position taking the last progression direction into account.

  std::vector< vnl_vector_fixed<float,3> > CreateDirections(unsigned int NPoints);

  void BeforeTracking();
  void AfterTracking();

  PolyDataType                        m_FiberPolyData;
  vtkSmartPointer<vtkPoints>          m_Points;
  vtkSmartPointer<vtkCellArray>       m_Cells;
  BundleType                          m_Tractogram;
  BundleType                          m_GmStubs;

  ItkFloatImgType::Pointer            m_StoppingRegions;
  ItkFloatImgType::Pointer            m_TargetRegions;
  ItkFloatImgType::Pointer            m_SeedImage;
  ItkFloatImgType::Pointer            m_MaskImage;
  ItkFloatImgType::Pointer            m_ExclusionRegions;
  ItkDoubleImgType::Pointer           m_OutputProbabilityMap;

  bool                                m_Verbose;
  std::vector< itk::Point<float> >    m_SeedPoints;
  std::atomic<unsigned int>           m_CurrentTracts;  ///< read by GetStatusText from the GUI thread
  std::atomic<unsigned int>           m_Progress;       ///< seed points processed, read by GetStatusText
  std::atomic<unsigned int>           m_NumSeeds;
  std::atomic<bool>                   m_StopTracking;

  void BuildFibers();
  float CheckCurvature(DirectionContainer *fib, bool front);

  mitk::TrackingDataHandler*          m_TrackingHandler;

  std::chrono::time_point<std::chrono::system_clock> m_StartTime;
  std::chrono::time_point<std::chrono::system_clock> m_EndTime;

  itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer   m_MaskInterpolator;
  itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer   m_StopInterpolator;
  itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer   m_TargetInterpolator;
  itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer   m_SeedInterpolator;
  itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer   m_ExclusionInterpolator;
  bool                                                                     m_SeedImageSet;
  bool                                                                     m_TargetImageSet;
  EndpointConstraints                                                      m_EpConstraints;  ///< parameters' constraint, or the automatic one for a set target image
  std::vector< vnl_vector_fixed<float,3> >                                 m_ProbeVectors;   ///< neighbourhood sampling directions

  std::shared_ptr< mitk::StreamlineTractographyParameters > m_Parameters;


  // Directional priors
  mitk::TrackingDataHandler*                  m_TrackingPriorHandler;

};

}

#endif

