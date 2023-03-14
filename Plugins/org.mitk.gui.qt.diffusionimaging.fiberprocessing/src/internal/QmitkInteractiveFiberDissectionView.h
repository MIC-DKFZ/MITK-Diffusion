/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/



#ifndef QmitkInteractiveFiberDissectionView_h
#define QmitkInteractiveFiberDissectionView_h
#include "ui_QmitkInteractiveFiberDissectionViewControls.h"

#include <QInputDialog> //Pointset
#include <QmitkPointListWidget.h>//Pointset

#include <QmitkAbstractView.h>//Pointset
#include <QmitkSingleNodeSelectionWidget.h>//Pointset


#include <mitkFiberBundle.h>
#include <mitkDataInteractor.h>
#include <mitkIRenderWindowPartListener.h> //Pointset
#include <mitkStreamlineInteractor.h>
#include <mitkStreamlineInteractorBrush.h>
#include <mitkStreamlineFeatureExtractor.h>

#include <itkCastImageFilter.h>
#include <itkVTKImageImport.h>
#include <itkVTKImageExport.h>
#include <itkRegionOfInterestImageFilter.h>

#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkImageImport.h>
#include <vtkImageExport.h>
#include <vtkImageStencil.h>
#include <vtkSmartPointer.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelectedThresholds.h>
#include <vtkFloatArray.h>
#include <vtkCellPicker.h>

#include <mitkSurface.h>

#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>

/*!
\brief View to process fiber bundles. Supplies methods to extract fibers from the bundle, fiber resampling, mirroring, join and subtract bundles and much more.
*/
class QmitkInteractiveFiberDissectionView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  typedef itk::Image< unsigned char, 3 >    ItkUCharImageType;
  typedef itk::Image< float, 3 >            ItkFloatImageType;

  static const std::string VIEW_ID;

  QmitkInteractiveFiberDissectionView();
  virtual ~QmitkInteractiveFiberDissectionView();

  virtual void CreateQtPartControl(QWidget *parent) override;

  ///
  /// Sets the focus to an internal widget.
  ///
  virtual void SetFocus() override;

protected slots:

  void RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart) override; //Pointset
  void RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart) override; //Pointset
  void OnAddPointSetClicked();//Pointset
  void CreateStreamline();
  void RemovefromBundle( bool checked );
  void RemovefromBundleBrush( bool checked );
  void ExtractRandomFibersFromTractogram();
  void StartAlgorithm();
  void CreatePredictionNode();
  void CreateCertainNode();
  void CreateUncertaintySampleNode();
  void CreateDistanceSampleNode();
  void RemovefromUncertainty( bool checked );
  void RemovefromUncertaintyBrush( bool checked );
  void RemovefromDistance( bool checked );
  void RemovefromPrediction( bool checked );
  void RemovefromPredictionBrush( bool checked );
  void RemovefromSelection( bool checked );
  void ResampleTractogram();
  void RandomPrototypes();
  void SFFPrototypes();
  void StartValidation();
  void AutomaticLabelling();
  void RemoveCertainData();
  void ResetClassifier();



  void UpdateGui();     ///< update button activity etc. dpending on current datamanager selection




protected:
  void OnCurrentSelectionChanged(QmitkSingleNodeSelectionWidget::NodeList nodes);//Pointset


  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  void OnEndInteraction();
  void CreateStreamlineInteractor();
  void CreateStreamlineInteractorBrush();
  void CleanTestArray();

  Ui::QmitkInteractiveFiberDissectionViewControls* m_Controls;



  int m_IterationCounter;                                      ///< used for data node naming
  int m_RandomExtractionCounter;                               ///< used for random extracton of different Streamlines
  int m_activeCycleCounter;
  int m_newpos;
  int m_newneg;
  int m_startneg;
  int m_createdStreamlineCounter;
  int m_uncCounter;
//  int m_thresh2;


  std::vector<mitk::DataNode::Pointer>  m_SelectedFB;       ///< selected fiber bundle nodes
  mitk::DataNode::Pointer m_testnode;
//  mitk::DataNode::Pointer              m_trainbundle;
  mitk::Image::Pointer                  m_SelectedImage;
  mitk::DataNode::Pointer               m_SelectedPS;
  mitk::DataNode::Pointer               m_SelectedImageNode;
  mitk::FiberBundle::Pointer            m_positiveBundle;
  mitk::FiberBundle::Pointer            m_newfibersBundle;
  mitk::FiberBundle::Pointer            m_negativeBundle;
  mitk::FiberBundle::Pointer            m_Prediction;
  mitk::FiberBundle::Pointer            m_CertainPlus;
  mitk::FiberBundle::Pointer            m_CertainMinus;
  mitk::FiberBundle::Pointer            m_CertainBetweenPlus;
  mitk::FiberBundle::Pointer            m_CertainBetweenMinus;
  mitk::FiberBundle::Pointer            m_UncertaintyLabel;
  mitk::FiberBundle::Pointer            m_DistanceLabel;


  mitk::DataNode::Pointer               m_positiveBundleNode;
  mitk::DataNode::Pointer               m_newfibersBundleNode;
  mitk::DataNode::Pointer               m_negativeBundleNode;
  mitk::DataNode::Pointer               m_PredictionNode;
  mitk::DataNode::Pointer               m_CertainPlusNode;
  mitk::DataNode::Pointer               m_CertainMinusNode;
  mitk::DataNode::Pointer               m_CertainBetweenPlusNode;
  mitk::DataNode::Pointer               m_CertainBetweenMinusNode;
  mitk::DataNode::Pointer               m_UncertaintyLabelNode;
  mitk::DataNode::Pointer               m_DistanceLabelNode;

  vtkSmartPointer<vtkPolyData>          m_positiveFibersData;
  vtkSmartPointer<vtkPolyData>          m_newfibersFibersData;

  vtkSmartPointer<vtkCellPicker>        m_picker1;
  mitk::StreamlineInteractor::Pointer   m_StreamlineInteractor;
  mitk::StreamlineInteractorBrush::Pointer   m_StreamlineInteractorBrush;

  std::shared_ptr< mitk::StreamlineFeatureExtractor > classifier;
  std::shared_ptr< mitk::StreamlineFeatureExtractor > validater;


  std::vector<vnl_vector<float>> m_metrics;

  std::vector<std::vector<unsigned int>> m_index;

};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

