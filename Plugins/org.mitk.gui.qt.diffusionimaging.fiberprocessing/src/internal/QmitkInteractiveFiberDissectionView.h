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


#include <QmitkAbstractView.h>
#include <mitkFiberBundle.h>
#include <mitkDataInteractor.h>
#include <mitkStreamlineInteractor.h>
#include <mitkSphereInteractor.h>
#include <mitkSphereSegmentationInteractor.h>
#include <mitkStreamlineInteractorBrush.h>
#include <mitkStreamlineFeatureExtractor.h>
#include <mitkStatisticModel.h> 

#include <mitkTensorImage.h>
#include <mitkDataStorage.h>
#include <mitkDataNode.h>
#include <mitkImage.h>
#include <itkImage.h>
#include <random>
#include <mitkPointSet.h>
#include <mitkPointSetShapeProperty.h>
#include <QTimer>
#include <QmitkStdMultiWidget.h>
#include <QmitkSliceNavigationListener.h>
#include <mitkILifecycleAwarePart.h>
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
#include <vtkExtractSelection.h>
#include <vtkFloatArray.h>
#include <vtkCellPicker.h>

#include <mitkSurface.h>

#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>

/*!
\brief View to dissect fiber bundles with active learning
*/

class QmitkInteractiveFiberDissectionView : public QmitkAbstractView
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



  virtual void SetFocus() override;

protected slots:

  void UpdateGui();     ///< update button activity etc. dpending on current datamanager selection
  void RemovefromBundleBrush( bool checked );
  void CreateSubset();
  void CreateReductionSubset();
  void CreateSphereSegmentation();
  void ExtractRandomFibersFromTractogram();
  void StartAlgorithm();
  void StartAlgorithmGen();
  void StartGeneralizationAlgorithm();
  void CreatePredictionNode();
  void CreateUncertaintySampleNode();
  void CreateDistanceSampleNode();
  void RemovefromUncertaintyBrush( bool checked );
  void RemovefromDistanceBrush( bool checked);
  void RemovefromPredictionBrush( bool checked );
  void ResampleTractogram();
  void SFFPrototypes();
  void CreateUncertantyMap();
  void ResetClassifier();


  void addGeneralizationBundles(bool reset);
  void addGeneralizationSavedPosBundles();
  void addGeneralizationSavedNegBundles();





protected:


  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  void OnEndInteraction();
  void CreateStreamlineInteractor();
  void CreateStreamlineInteractorBrush();
  void StartLabelingBrush(mitk::DataNode::Pointer toLabelNode, bool fromPrediction);
  void StopLabelingBrush();
  /** Trains the classifier on the given labels and shows the prediction. Starts a new session in the first cycle. */
  bool RunClassifier(mitk::FiberBundle::Pointer plus, mitk::FiberBundle::Pointer minus);
  void AddSavedGeneralizationBundle(const mitk::DataNode *selectedNode, mitk::DataNode::Pointer &genNode, const std::string &name, float r, float g, float b);

  Ui::QmitkInteractiveFiberDissectionViewControls* m_Controls;



  int m_IterationCounter;                                      ///< used for data node naming
  int m_RandomExtractionCounter;                               ///< used for random extracton of different Streamlines
  int m_activeCycleCounter;
  int m_newpos;
  int m_newneg;
  int m_startneg;
  int m_createdStreamlineCounter;
  int m_uncCounter;
  int m_prototypecounter;
  int m_posOffsetcounter;
  int m_init_subset;
  int m_generalizationMode;




  std::vector<mitk::DataNode::Pointer>  m_SelectedFB;       ///< selected fiber bundle nodes
  mitk::DataNode::Pointer               m_testnode;
  mitk::Image::Pointer                  m_SelectedImage;
  mitk::DataNode::Pointer               m_SelectedPS;
  mitk::DataNode::Pointer               m_SelectedImageNode;
  mitk::FiberBundle::Pointer            m_positiveBundle;
  mitk::FiberBundle::Pointer            m_newfibersBundle;
  mitk::FiberBundle::Pointer            m_negativeBundle;
  mitk::FiberBundle::Pointer            m_Prediction;
  mitk::FiberBundle::Pointer            m_UncertaintyLabel;
  mitk::FiberBundle::Pointer            m_DistanceLabel;
  mitk::FiberBundle::Pointer            cur_negativeBundle;
  mitk::FiberBundle::Pointer            cur_positiveBundle;


  mitk::DataNode::Pointer               m_positiveBundleNode;
  mitk::DataNode::Pointer               m_newfibersBundleNode;
  mitk::DataNode::Pointer               m_negativeBundleNode;
  mitk::DataNode::Pointer               m_PredictionNode;
  mitk::DataNode::Pointer               m_UncertaintyLabelNode;
  mitk::DataNode::Pointer               m_DistanceLabelNode;
  mitk::DataNode::Pointer               m_reducedFibersDataNode;
  mitk::DataNode::Pointer               m_reducedFibersTractogramNode;

  mitk::DataNode::Pointer               m_Prototypes;

 // This code is for atTRACTive generalization
  mitk::FiberBundle::Pointer            m_positiveGenBundle;
  mitk::FiberBundle::Pointer            m_negativeGenBundle;
  mitk::DataNode::Pointer               m_positiveGenBundleNode;
  mitk::DataNode::Pointer               m_negativeGenBundleNode;




  vtkSmartPointer<vtkPolyData>          m_positiveFibersData;
  vtkSmartPointer<vtkPolyData>          m_newfibersFibersData;

  mitk::DataNode::Pointer               circleNode;

  mitk::StreamlineInteractorBrush::Pointer   m_StreamlineInteractorBrush;

  mitk::SphereInteractor::Pointer   m_SphereInteractor;
  mitk::SphereSegmentationInteractor::Pointer   m_SphereSegmentationInteractor;

  std::shared_ptr< mitk::StreamlineFeatureExtractor > classifier;



};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

