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
  void CreateSphereSegmentation();
  void ExtractRandomFibersFromTractogram();
  void StartAlgorithm();
  void CreatePredictionNode();
  void CreateUncertaintySampleNode();
  void CreateDistanceSampleNode();
  void RemovefromUncertaintyBrush( bool checked );
  void RemovefromDistanceBrush( bool checked );
  void RemovefromPredictionBrush( bool checked );
  void ResampleTractogram();
  void SFFPrototypes();
  void CreateUncertantyMap();
  void ResetClassifier();






protected:


  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  void OnEndInteraction();
  void CreateStreamlineInteractor();
  void CreateStreamlineInteractorBrush();

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


  mitk::DataNode::Pointer               m_positiveBundleNode;
  mitk::DataNode::Pointer               m_newfibersBundleNode;
  mitk::DataNode::Pointer               m_negativeBundleNode;
  mitk::DataNode::Pointer               m_PredictionNode;
  mitk::DataNode::Pointer               m_UncertaintyLabelNode;
  mitk::DataNode::Pointer               m_DistanceLabelNode;
  mitk::DataNode::Pointer               m_reducedFibersDataNode;
  mitk::DataNode::Pointer               m_Prototypes;
  mitk::DataNode::Pointer               m_newtestnode;

  vtkSmartPointer<vtkPolyData>          m_positiveFibersData;
  vtkSmartPointer<vtkPolyData>          m_newfibersFibersData;

  mitk::DataNode::Pointer               circleNode;

  mitk::StreamlineInteractorBrush::Pointer   m_StreamlineInteractorBrush;

  mitk::SphereInteractor::Pointer   m_SphereInteractor;
  mitk::SphereSegmentationInteractor::Pointer   m_SphereSegmentationInteractor;

  std::shared_ptr< mitk::StreamlineFeatureExtractor > classifier;


  std::vector<std::vector<unsigned int>> m_index;

};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

