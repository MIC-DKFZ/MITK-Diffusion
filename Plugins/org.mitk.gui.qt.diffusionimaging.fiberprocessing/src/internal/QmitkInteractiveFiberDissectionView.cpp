/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/



// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchPart.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkInteractiveFiberDissectionView.h"
//#include "QmitkStreamlineTrackingView.h"

#include <QMessageBox>

#include <mitkImageCast.h>
#include <mitkPointSet.h>
#include <mitkImageAccessByItk.h>
#include <mitkDataNodeObject.h>
#include <mitkTensorImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkPlanarCircle.h>
#include <mitkSurface.h>
#include <mitkStatisticModel.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
//#include <mitkPlanarFigureInteractor.h>
//#include <mitkStreamlineFeatureExtractor.h>

#include <mitkInteractionConst.h>
#include "usModuleRegistry.h"
//#include <itkFiberCurvatureFilter.h>


#include <itkResampleImageFilter.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkTractsToFiberEndingsImageFilter.h>
#include <itkTractDensityImageFilter.h>
#include <itkImageRegion.h>
#include <itkTractsToRgbaImageFilter.h>
#include <itkFiberExtractionFilter.h>

#include <mitkInteractionEventObserver.h>
#include <vtkCellPicker.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <vnl/vnl_sparse_matrix.h>
#include <numeric>
#include <random>

#include <thread>

#include <QmitkStdMultiWidget.h>

#include <mitkBaseRenderer.h>
#include <mitkRenderingManager.h>
#include <mitkVtkPropRenderer.h>
#include <mitkInteractionConst.h>
#include <mitkInteractionPositionEvent.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <QFileDialog>

const std::string QmitkInteractiveFiberDissectionView::VIEW_ID = "org.mitk.views.interactivefiberdissection";
const std::string id_DataManager = "org.mitk.views.datamanager";
using namespace mitk;

QmitkInteractiveFiberDissectionView::QmitkInteractiveFiberDissectionView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_IterationCounter(0)
  , m_RandomExtractionCounter(0)
  , m_activeCycleCounter(0)
  , m_createdStreamlineCounter(0)
  , m_prototypecounter(0)
  , m_posOffsetcounter(0)
  , m_init_subset(0)
  , m_generalizationMode(0)
{

}

// Destructor
QmitkInteractiveFiberDissectionView::~QmitkInteractiveFiberDissectionView()
{

}


void QmitkInteractiveFiberDissectionView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkInteractiveFiberDissectionViewControls;
    m_Controls->setupUi( parent );

    m_Controls->m_TestBundleBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isTestBundle = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_TestBundleBox->SetPredicate( isTestBundle);

    m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::Image>::Pointer isImage = mitk::TNodePredicateDataType<mitk::Image>::New();
    m_Controls->m_ImageBox->SetPredicate( isImage);


    m_Controls->m_posGenBundleBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isposGenBundle = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_posGenBundleBox->SetPredicate( isposGenBundle);


    m_Controls->m_negGenBundleBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isnegGenBundle = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_negGenBundleBox->SetPredicate( isnegGenBundle);


    connect(m_Controls->m_BrushButton, SIGNAL(toggled(bool)), this, SLOT( RemovefromBundleBrush(bool) ) );

    connect(m_Controls->m_SubsetCreation, SIGNAL( clicked() ), this, SLOT( CreateSubset()));

    connect(m_Controls->m_RoiReduction, SIGNAL( clicked() ), this, SLOT( CreateReductionSubset()));

    connect(m_Controls->m_CreateSegmentation, SIGNAL( clicked() ), this, SLOT( CreateSphereSegmentation()));

    connect(m_Controls->m_AddRandomFibers, SIGNAL( clicked() ), this, SLOT( ExtractRandomFibersFromTractogram() ) ); //need

    connect(m_Controls->m_TrainClassifier, SIGNAL( clicked() ), this, SLOT( StartAlgorithm( )));

    connect(m_Controls->m_AddUncertainFibers, SIGNAL( clicked() ), this, SLOT( CreateUncertaintySampleNode( )));

    connect(m_Controls->m_AddDistanceFibers, SIGNAL( clicked() ), this, SLOT( CreateDistanceSampleNode( )));

    connect(m_Controls->m_AddDistanceFibersGen, SIGNAL( clicked() ), this, SLOT( CreateDistanceSampleNode( )));

    connect(m_Controls->m_unclabelingBrush, SIGNAL(toggled(bool)), this, SLOT( RemovefromUncertaintyBrush(bool) ) ); //need

    connect(m_Controls->m_distlabelingBrush, SIGNAL(toggled(bool)), this, SLOT( RemovefromDistanceBrush(bool) ) ); //need

    connect(m_Controls->m_distlabelingBrushGen, SIGNAL(toggled(bool)), this, SLOT( RemovefromDistanceBrush(bool) ) ); //need

    connect(m_Controls->m_predlabelingBrush, SIGNAL(toggled(bool)), this, SLOT( RemovefromPredictionBrush(bool) ) ); //need
    connect(m_Controls->m_predlabelingBrushGen, SIGNAL(toggled(bool)), this, SLOT( RemovefromPredictionBrush(bool) ) ); //need

    // the generalization tab offers the same sampling tools as the normal workflow
    connect(m_Controls->m_AddUncertainFibersGen, SIGNAL( clicked() ), this, SLOT( CreateUncertaintySampleNode( )));
    connect(m_Controls->m_unclabelingBrushGen, SIGNAL(toggled(bool)), this, SLOT( RemovefromUncertaintyBrush(bool) ) );
    connect(m_Controls->m_CreateUncertantyManGen, SIGNAL( clicked() ), this, SLOT( CreateUncertantyMap( ) ) );

    connect(m_Controls->m_ResampleButton, SIGNAL( clicked() ), this, SLOT( ResampleTractogram( ) ) );

    connect(m_Controls->m_CreateUncertantyMap, SIGNAL( clicked() ), this, SLOT( CreateUncertantyMap( ) ) );

    connect(m_Controls->m_resetClassifier, SIGNAL( clicked() ), this, SLOT( ResetClassifier( ) ) );

    connect(m_Controls->m_add_annotations, &QPushButton::clicked, this, [this]() {addGeneralizationBundles(true);});
    connect(m_Controls->m_add_annotations_2, &QPushButton::clicked, this, [this]() {addGeneralizationBundles(true);});


    connect(m_Controls->m_add_positiveFibers, &QPushButton::clicked, this, [this]() {addGeneralizationSavedPosBundles();});
    connect(m_Controls->m_add_negativeFibers, &QPushButton::clicked, this, [this]() {addGeneralizationSavedNegBundles();});

    connect(m_Controls->m_initGeneralization, SIGNAL( clicked() ), this, SLOT( StartGeneralizationAlgorithm( )));

    connect(m_Controls->m_TrainClassifierGen, SIGNAL( clicked() ), this, SLOT( StartAlgorithmGen( )));


  }

  UpdateGui();
}

void QmitkInteractiveFiberDissectionView::SetFocus()
{
  m_Controls->toolBoxx->setFocus();
  UpdateGui();
}

void QmitkInteractiveFiberDissectionView::UpdateGui()
{
  m_Controls->m_FibLabel->setText("<font color='red'>nothing selected</font>");
  m_Controls->m_InputData->setTitle("Please Select Input Data");

  // disable everything, then enable what is possible. Both tabs have the same tools.
  for (QPushButton* brush : {m_Controls->m_BrushButton,
                             m_Controls->m_unclabelingBrush, m_Controls->m_unclabelingBrushGen,
                             m_Controls->m_predlabelingBrush, m_Controls->m_predlabelingBrushGen,
                             m_Controls->m_distlabelingBrush, m_Controls->m_distlabelingBrushGen})
  {
    brush->setCheckable(true);
    brush->setEnabled(false);
  }
  for (QWidget* widget : std::initializer_list<QWidget*>{
         m_Controls->m_SubsetCreation,
         m_Controls->m_TrainClassifier, m_Controls->m_TrainClassifierGen,
         m_Controls->m_add_annotations, m_Controls->m_add_annotations_2,
         m_Controls->m_CreateUncertantyMap, m_Controls->m_CreateUncertantyManGen,
         m_Controls->m_Numtolabel, m_Controls->m_NumtolabelGen,
         m_Controls->m_AddDistanceFibers, m_Controls->m_AddDistanceFibersGen,
         m_Controls->m_NumRandomFibers2, m_Controls->m_NumRandomFibers2Ben,
         m_Controls->m_AddUncertainFibers, m_Controls->m_AddUncertainFibersGen})
  {
    widget->setEnabled(false);
  }


  mitk::DataNode::Pointer curtestnode =  m_Controls->m_TestBundleBox->GetSelectedNode();
  if (m_activeCycleCounter>0 && curtestnode!=m_testnode){
   MITK_INFO << "Testnode changed";
   m_init_subset = 0;
    m_activeCycleCounter=0;

    // Checking if m_positiveBundleNode is not empty
    if (m_positiveBundleNode.IsNotNull())
    {
        // Accessing the fiber polydata and updating the offset counter
        mitk::FiberBundle::Pointer fiberBundle = dynamic_cast<mitk::FiberBundle*>(m_positiveBundleNode->GetData());
        if (fiberBundle.IsNotNull())
        {
            m_posOffsetcounter += fiberBundle->GetFiberPolyData()->GetNumberOfCells();
        }
}

  }


  bool testnodeSelected = curtestnode.IsNotNull();
  bool fibSelected = !m_SelectedFB.empty();
  bool multipleFibsSelected = (m_SelectedFB.size()>1);
  bool posSelected = this->GetDataStorage()->Exists(m_positiveBundleNode);
  bool posGenSelected = this->GetDataStorage()->Exists(m_positiveGenBundleNode);
  bool negSelected = this->GetDataStorage()->Exists(m_negativeBundleNode);
  bool indexSelected = classifier && !classifier->GetCandidateIds().empty();
  bool uncertaintySelected = this->GetDataStorage()->Exists(m_UncertaintyLabelNode);
  bool distanceuncertanySelected = this->GetDataStorage()->Exists(m_DistanceLabelNode);
  bool predictionSelected = this->GetDataStorage()->Exists(m_PredictionNode);
  bool classifierexist = (classifier != nullptr);




  // are fiber bundles selected?
  if ( testnodeSelected)
  {

      m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
      m_Controls->m_SubsetCreation->setEnabled(true);
      m_Controls->m_AddRandomFibers->setEnabled(true);
      m_Controls->m_NumRandomFibers->setEnabled(true);
  }

  
  if (fibSelected)
  {
    m_Controls->m_FibLabel->setText(QString(m_SelectedFB.at(0)->GetName().c_str()));
    if (multipleFibsSelected )
    {
      m_Controls->m_FibLabel->setText("multiple bundles selected");
    }
  }



  if (posSelected && negSelected)
  {
      m_Controls->m_TrainClassifier->setEnabled(true);
  }

  if (posGenSelected)
  {
      m_Controls->m_initGeneralization->setEnabled(true);
  }

  if (posGenSelected && posSelected)
  {
       m_Controls->m_TrainClassifierGen->setEnabled(true);
  }

  // tools that work on the current classifier are the same in both tabs
  auto enable = [](bool enabled, std::initializer_list<QWidget*> widgets) {
    for (auto* widget : widgets)
      widget->setEnabled(enabled);
  };
  enable(indexSelected, {m_Controls->m_AddUncertainFibers, m_Controls->m_AddUncertainFibersGen,
                         m_Controls->m_NumRandomFibers2, m_Controls->m_NumRandomFibers2Ben,
                         m_Controls->m_AddDistanceFibers, m_Controls->m_AddDistanceFibersGen,
                         m_Controls->m_Numtolabel, m_Controls->m_NumtolabelGen});
  enable(uncertaintySelected, {m_Controls->m_unclabelingBrush, m_Controls->m_unclabelingBrushGen});
  enable(distanceuncertanySelected, {m_Controls->m_distlabelingBrush, m_Controls->m_distlabelingBrushGen});
  enable(predictionSelected, {m_Controls->m_predlabelingBrush, m_Controls->m_predlabelingBrushGen});
  enable(predictionSelected && indexSelected, {m_Controls->m_CreateUncertantyMap, m_Controls->m_CreateUncertantyManGen});
  enable(classifierexist, {m_Controls->m_add_annotations, m_Controls->m_add_annotations_2});
}

void QmitkInteractiveFiberDissectionView::OnEndInteraction()
{

}

namespace
{
  mitk::FiberBundle::Pointer GetBundle(const mitk::DataNode *node)
  {
    return node != nullptr ? dynamic_cast<mitk::FiberBundle*>(node->GetData()) : nullptr;
  }

  /** Creates a node for a bundle of labeled or to be labeled streamlines. */
  mitk::DataNode::Pointer CreateBundleNode(const std::string &name, mitk::FiberBundle::Pointer bundle)
  {
    auto node = mitk::DataNode::New();
    node->SetData(bundle);
    node->SetName(name);
    node->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(5.0));
    return node;
  }

  /** Concatenates two bundles into a new one; null bundles count as empty. */
  mitk::FiberBundle::Pointer MergeBundles(mitk::FiberBundle::Pointer a, mitk::FiberBundle::Pointer b)
  {
    if (a.IsNull())
      return b.IsNotNull() ? b->GetDeepCopy() : mitk::FiberBundle::New();
    return a->AddBundle(b);
  }
}

void QmitkInteractiveFiberDissectionView::ResampleTractogram()
{
  mitk::DataNode::Pointer node = m_Controls->m_TestBundleBox->GetSelectedNode();
  auto tractogram = GetBundle(node);
  if (tractogram.IsNull())
    return;

  // shuffled copy, so that the order of the streamlines in the file does not matter
  std::vector<unsigned int> ids(tractogram->GetFiberPolyData()->GetNumberOfCells());
  std::iota(ids.begin(), ids.end(), 0);
  std::shuffle(ids.begin(), ids.end(), std::mt19937(std::random_device()()));

  auto resampled = mitk::StreamlineFeatureExtractor::SubsetByIds(tractogram, ids);
  resampled->ResampleToNumPoints(mitk::StreamlineFeatureExtractor::NumPoints);
  MITK_INFO << "Resampling Done";

  mitk::DataNode::Pointer newnode = mitk::DataNode::New();
  newnode->SetData(resampled);
  newnode->SetName(node->GetName() + "_" + std::to_string(mitk::StreamlineFeatureExtractor::NumPoints));
  this->GetDataStorage()->Add(newnode, node);
  node->SetVisibility(false);
  m_Controls->m_TestBundleBox->SetSelectedNode(newnode);

  UpdateGui();
}

void QmitkInteractiveFiberDissectionView::SFFPrototypes()
{
  if (m_prototypecounter != 0)
  {
    MITK_INFO << "No Prototypes created as classifier continues to learn. To start from scratch, press 'Reset Classifier'";
    return;
  }

  m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
  const unsigned int numPrototypes = m_Controls->m_NumPrototypes->value();
  MITK_INFO << "Create " << numPrototypes << " prototypes based on SFF";

  // subset size of the subset farthest first (SFF) algorithm
  const auto subsetSize = static_cast<unsigned int>(std::max(1.0, std::ceil(3.0 * numPrototypes * std::log(numPrototypes))));
  auto fib = GetBundle(m_testnode);
  auto ids = mitk::StreamlineFeatureExtractor::SelectDiverseFibers(fib, numPrototypes, subsetSize);

  m_Prototypes = mitk::DataNode::New();
  m_Prototypes->SetData(mitk::StreamlineFeatureExtractor::SubsetByIds(fib, ids));
  m_Prototypes->SetName("Prototypes");
  m_Prototypes->SetVisibility(false);
  m_prototypecounter = 1;
}

void QmitkInteractiveFiberDissectionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& nodes)
{

    m_SelectedFB.clear();

    if (nodes.empty() || nodes.front().IsNull())
    {
      m_SelectedImageNode = nullptr;
    }
    else
    {
      m_SelectedImageNode = nodes.front();
    }



    for (auto node: nodes)
    {
      if (dynamic_cast<mitk::Image*>(node->GetData()))
          m_SelectedImage = dynamic_cast<mitk::Image*>(node->GetData());
      else if ( dynamic_cast<mitk::FiberBundle*>(node->GetData()) )
          m_SelectedFB.push_back(node);
    }



    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateSphereSegmentation()
{
    if (m_Controls->m_helpmessages->isChecked())
    {
        // Nothing selected. Inform the user and return
        QMessageBox::information(nullptr, "Helper", "Roughly draw a start and an end ROI of your target tract.\n"
                                                    " 1. Use the middle mouse button to set the center of the start region, then adjust the size and confirm with the middle mouse button\n"
                                                    " 2. Navigate to the end region and do the same again. \n"
                                                    " 3. The Segmentation of ROIS will be presented.\n"
                                                    " 4. Choose how many fibers of the subset you want to label.\n"
                                                    " 5. Click 'add' to annotate a subset of presented fibers.");
    }

    mitk::DataNode::Pointer imageNode = m_Controls->m_ImageBox->GetSelectedNode();

    // Check if imageNode is nullptr
    if (!imageNode)
    {
        QMessageBox::critical(nullptr, "Error", "No image selected. Please select an image and try again.");
        return;
    }

    mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(imageNode->GetData());

    // Check if the image pointer is nullptr
    if (!image)
    {
        QMessageBox::critical(nullptr, "Error", "Selected node does not contain a valid image. Please select a valid image and try again.");
        return;
    }

    m_SphereSegmentationInteractor = mitk::SphereSegmentationInteractor::New();
    m_SphereSegmentationInteractor->LoadStateMachine("SphereSegmentationInteractionsStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_SphereSegmentationInteractor->SetEventConfig("SphereSegmentationInteractionsConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));

    mitk::DataNode::Pointer startDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer startSurfaceDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer endDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer endSurfaceDataNode = mitk::DataNode::New();

    startDataNode->SetName("StartRegion");
    startDataNode->SetColor(1, 0, 1);
    startSurfaceDataNode->SetName("StartROI");
    startSurfaceDataNode->SetColor(1, 0, 1);

    endDataNode->SetName("EndRegion");
    endDataNode->SetColor(1, 1, 0);
    endSurfaceDataNode->SetName("EndROI");
    endSurfaceDataNode->SetColor(1, 1, 0);

    this->GetDataStorage()->Add(startDataNode);
    this->GetDataStorage()->Add(endDataNode);
    this->GetDataStorage()->Add(startSurfaceDataNode, startDataNode);
    this->GetDataStorage()->Add(endSurfaceDataNode, endDataNode);

    m_SphereSegmentationInteractor->SetImage(image);
    m_SphereSegmentationInteractor->StartEndNodes(startDataNode, endDataNode);
    m_SphereSegmentationInteractor->StartEndSurfaceNodes(startSurfaceDataNode, endSurfaceDataNode);

    // this->GetDataStorage()->Remove(startDataNode);
    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateReductionSubset()
{

    if (m_reducedFibersTractogramNode.IsNotNull())
    {
      this->GetDataStorage()->Remove(m_reducedFibersTractogramNode);
    }

    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
    // Create and configure the SphereInteractor
    m_SphereInteractor = mitk::SphereInteractor::New();
    m_SphereInteractor->LoadStateMachine("SphereInteractionsStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_SphereInteractor->SetEventConfig("SphereInteractionsConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));

    MITK_INFO << "Interactor created";

    // Create start and end data nodes
    mitk::DataNode::Pointer startDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer endDataNode = mitk::DataNode::New();

    startDataNode->SetName("StartRegion");
    startDataNode->SetColor(1, 0, 1);

    endDataNode->SetName("EndRegion");
    endDataNode->SetColor(1, 1, 0);

    // Check if m_testnode is valid before adding nodes to data storage
    if (!m_testnode)
    {
        MITK_ERROR << "m_testnode is null! Cannot add start and end nodes to data storage.";
        return;
    }

    // Add start and end nodes as children of m_testnode
    this->GetDataStorage()->Add(startDataNode, m_testnode);
    this->GetDataStorage()->Add(endDataNode, m_testnode);

    MITK_INFO << "Create bundle";

    // Create and configure the reduced fibers data node
    m_reducedFibersTractogramNode = mitk::DataNode::New();
    mitk::FiberBundle::Pointer reducedFib = mitk::FiberBundle::New();

    if (!reducedFib)
    {
        MITK_ERROR << "Failed to create FiberBundle.";
        return;
    }

    m_reducedFibersTractogramNode->SetData(reducedFib);

    std::string roi_name = m_testnode->GetName();
    m_reducedFibersTractogramNode->SetName(roi_name + "_red");

    // Add the reduced fibers node to the data storage
    this->GetDataStorage()->Add(m_reducedFibersTractogramNode);

    // Ensure m_SphereInteractor and the nodes are valid before proceeding
    if (!m_SphereInteractor || !m_reducedFibersTractogramNode)
    {
        MITK_ERROR << "SphereInteractor or reducedFibersDataNode is null!";
        return;
    }

    m_SphereInteractor->workingBundleNode(m_testnode, m_reducedFibersTractogramNode);

    // Start the interaction using the start and end nodes
    m_SphereInteractor->StartEndNodes(startDataNode, endDataNode);

    MITK_INFO << "Done";
    // Update the GUI
    m_init_subset=2;
    m_Controls->m_TestBundleBox->SetSelectedNode(m_reducedFibersTractogramNode);
    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateSubset()
{

    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
    if (m_Controls->m_helpmessages->isChecked())
    {
      // Nothing selected. Inform the user and return
      QMessageBox::information(nullptr, "Helper", "Roughly draw a start and an end ROI of your target tract.\n"
                                                  " 1. Use the middle mouse button two draw the center of the start region, then adjust the size and confirm with the middle mouse button\n"
                                                  " 2. Navigate to the end region and do the same again. \n"
                                                  " 3. The subset of streamlines will be presented.\n"
                                                  " 4. Choose how many fibers of the subset you want to label.\n"
                                                  " 5. Click 'add' to annotate a subset of presented fibers.");

    }
    if (m_reducedFibersDataNode.IsNotNull())
    {
      this->GetDataStorage()->Remove(m_reducedFibersDataNode);
    }

    m_SphereInteractor = mitk::SphereInteractor::New();
    m_SphereInteractor->LoadStateMachine("SphereInteractionsStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_SphereInteractor->SetEventConfig("SphereInteractionsConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));

    MITK_INFO << "Interarctor created";

    mitk::DataNode::Pointer startDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer endDataNode = mitk::DataNode::New();


    startDataNode->SetName("StartRegion");
    startDataNode->SetColor(1, 0, 1);

    endDataNode->SetName("EndRegion");
    endDataNode->SetColor(1, 1, 0);

    this->GetDataStorage()->Add(startDataNode, m_testnode);
    this->GetDataStorage()->Add(endDataNode, m_testnode);

    MITK_INFO << "Create bundle";
    m_reducedFibersDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    mitk::FiberBundle::Pointer reducedFib = mitk::FiberBundle::New();
    m_reducedFibersDataNode->SetData(reducedFib);
    m_reducedFibersDataNode->SetName("Reduced");

    this->GetDataStorage()->Add(m_reducedFibersDataNode, m_testnode);

    m_SphereInteractor->workingBundleNode(m_testnode, m_reducedFibersDataNode);
    m_SphereInteractor->StartEndNodes(startDataNode, endDataNode);


    MITK_INFO << "Done";
    m_Controls->m_AddRandomFibers->setEnabled(true);
    m_Controls->m_NumRandomFibers->setEnabled(true);
    m_init_subset = 1;

    UpdateGui();

}

void QmitkInteractiveFiberDissectionView::ExtractRandomFibersFromTractogram()
{
    if (m_Controls->m_helpmessages->isChecked())
    {
      QMessageBox::information(nullptr, "Helper", "Please annotate the displayed fibers.\n"
                                                  "- If a fiber belongs to the target tract press 'alt' and simply hover over the fiber\n"
                                                  "- If the fiber does not belong to the target tract press 'shift' and hover over the fiber\n"
                                                  "- To unlabel a fiber press the 'Cntrl' button and hover over the respective fiber \n"
                                                  "Once you are done, you can train the classifier by pressing 'Train Classifier'\n"
                                                  "Note: You dont have to label all Fibers "
                               );
    }

    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
    auto fib = GetBundle(m_testnode);
    if (fib.IsNull())
      return;

    const unsigned int numFibers = m_Controls->m_NumRandomFibers->value();
    MITK_INFO << "Number of Fibers to extract from Tractogram: " << numFibers;

    mitk::FiberBundle::Pointer toLabel;
    unsigned int numDiverse = numFibers;

    // After a start/end ROI subset was created, half of the fibers are drawn randomly from that subset
    auto reduced = GetBundle(m_reducedFibersDataNode);
    if (m_init_subset == 1 && reduced.IsNotNull())
    {
      m_reducedFibersDataNode->SetVisibility(false);

      std::vector<unsigned int> ids(reduced->GetFiberPolyData()->GetNumberOfCells());
      std::iota(ids.begin(), ids.end(), 0);
      std::shuffle(ids.begin(), ids.end(), std::mt19937(std::random_device()()));
      ids.resize(std::min<std::size_t>(numFibers / 2, ids.size()));

      toLabel = mitk::StreamlineFeatureExtractor::SubsetByIds(reduced, ids);
      numDiverse = numFibers - ids.size();
      m_init_subset = 0;
    }

    // the remaining fibers cover the whole tractogram (subset farthest first on a subset of 3 N ln N, at least 500 fibers)
    const auto subsetSize = std::max(500u, static_cast<unsigned int>(std::ceil(3.0 * numDiverse * std::log(std::max(numDiverse, 1u)))));
    try
    {
      auto ids = mitk::StreamlineFeatureExtractor::SelectDiverseFibers(fib, numDiverse, subsetSize);
      toLabel = MergeBundles(toLabel, mitk::StreamlineFeatureExtractor::SubsetByIds(fib, ids));
    }
    catch (const std::exception &e)
    {
      QMessageBox::warning(nullptr, "Fiber extraction", e.what());
      return;
    }

    m_testnode->SetVisibility(false);
    m_Controls->m_BrushButton->setChecked(false);
    this->GetDataStorage()->Remove(m_newfibersBundleNode);

    m_newfibersBundle = toLabel;
    m_newfibersBundle->SetFiberColors(255, 255, 255);
    m_newfibersBundleNode = CreateBundleNode("ToLabel", m_newfibersBundle);
    this->GetDataStorage()->Add(m_newfibersBundleNode);

    if (!m_negativeBundleNode)
    {
      m_negativeBundleNode = CreateBundleNode("-Bundle", mitk::FiberBundle::New());
      this->GetDataStorage()->Add(m_negativeBundleNode);
      m_negativeBundleNode->SetVisibility(false);
    }

    if (!m_positiveBundleNode)
    {
      m_positiveBundleNode = CreateBundleNode("+Bundle", mitk::FiberBundle::New());
      this->GetDataStorage()->Add(m_positiveBundleNode);
    }

    m_RandomExtractionCounter++;
    m_Controls->m_BrushButton->setEnabled(true);
    m_Controls->m_BrushButton->setChecked(true);

    UpdateGui();
}

bool QmitkInteractiveFiberDissectionView::HelpMessagesEnabled() const
{
    return (m_generalizationMode ? m_Controls->m_helpmessages_2 : m_Controls->m_helpmessages)->isChecked();
}

void QmitkInteractiveFiberDissectionView::UncheckLabelingBrushes()
{
    for (QPushButton* brush : {m_Controls->m_unclabelingBrush, m_Controls->m_unclabelingBrushGen,
                               m_Controls->m_predlabelingBrush, m_Controls->m_predlabelingBrushGen,
                               m_Controls->m_distlabelingBrush, m_Controls->m_distlabelingBrushGen})
      brush->setChecked(false);
}

void QmitkInteractiveFiberDissectionView::StartLabelingBrush(mitk::DataNode::Pointer toLabelNode, bool fromPrediction)
{
    this->CreateStreamlineInteractorBrush();
    m_StreamlineInteractorBrush->EnableInteraction(true);
    m_StreamlineInteractorBrush->LabelfromPrediction(fromPrediction);
    m_StreamlineInteractorBrush->SetNegativeNode(m_negativeBundleNode);
    m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
    m_StreamlineInteractorBrush->SetToLabelNode(toLabelNode);
}

void QmitkInteractiveFiberDissectionView::StopLabelingBrush()
{
    if (m_StreamlineInteractorBrush.IsNotNull())
      m_StreamlineInteractorBrush->EnableInteraction(false);
}

void QmitkInteractiveFiberDissectionView::RemovefromBundleBrush( bool checked )
{
    if (checked)
      this->StartLabelingBrush(m_newfibersBundleNode, false);
    else
      this->StopLabelingBrush();

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateStreamlineInteractorBrush()
{

    m_StreamlineInteractorBrush = mitk::StreamlineInteractorBrush::New();

    m_StreamlineInteractorBrush->LoadStateMachine("StreamlineBrush3DStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_StreamlineInteractorBrush->SetEventConfig("StreamlineBrush3DConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));


    UpdateGui();
}

bool QmitkInteractiveFiberDissectionView::RunClassifier(mitk::FiberBundle::Pointer plus, mitk::FiberBundle::Pointer minus)
{
    this->UncheckLabelingBrushes();
    this->GetDataStorage()->Remove(m_UncertaintyLabelNode);
    this->GetDataStorage()->Remove(m_DistanceLabelNode);
    this->GetDataStorage()->Remove(m_newfibersBundleNode);
    m_uncCounter = 0;

    try
    {
      // A session starts with the first cycle. The classifier keeps the pool and the cached features until the next session.
      if (m_activeCycleCounter == 0 || m_prototypecounter == 0 || !classifier)
      {
        this->SFFPrototypes();
        classifier = std::make_shared<mitk::StreamlineFeatureExtractor>();
        classifier->SetTractogramTest(GetBundle(m_testnode));
        classifier->SetTractogramPrototypes(GetBundle(m_Prototypes));
      }

      classifier->SetActiveCycle(m_activeCycleCounter);
      classifier->SetTractogramPlus(plus);
      classifier->SetTractogramMinus(minus);
      classifier->Update();
    }
    catch (const std::exception &e)
    {
      MITK_ERROR << e.what();
      QMessageBox::warning(nullptr, "Classifier", e.what());
      return false;
    }

    MITK_INFO << "Number of Cycles " << m_activeCycleCounter;
    MITK_INFO << "Algorithm run succesfully";

    this->CreatePredictionNode();
    m_activeCycleCounter += 1;
    return true;
}

void QmitkInteractiveFiberDissectionView::StartAlgorithm()
{
    if (m_Controls->m_helpmessages->isChecked())
    {
      QMessageBox::information(nullptr, "Helper", "Classifier training started...\n"
                                                  "Once it is finished, the prediction will be displayed.\n"
                                                  "If it is sufficient, you may save it. Feel free to exclude outliers, or to label more fibers to improve the prediction.\n"
                                                  "You can turn of this helper by unchecking the radie button in the lower part of the view."
                                                  );
    }
    m_generalizationMode = 0;
    m_Controls->m_NumRandomFibers->setEnabled(false);
    m_Controls->m_BrushButton->setEnabled(false);

    m_negativeBundle = GetBundle(m_negativeBundleNode);
    m_positiveBundle = GetBundle(m_positiveBundleNode);
    this->RunClassifier(m_positiveBundle, m_negativeBundle);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreatePredictionNode()
{
    MITK_INFO << "Create Prediction";

    m_Prediction = classifier->ExtractFibers(classifier->GetPredictedIds());

    // the positively labeled streamlines of the current tractogram are part of the prediction
    auto positive = GetBundle(m_positiveBundleNode);
    if (positive.IsNotNull())
    {
      std::vector<unsigned int> ids;
      const auto numPositive = static_cast<unsigned int>(positive->GetFiberPolyData()->GetNumberOfCells());
      for (unsigned int i = m_posOffsetcounter; i < numPositive; ++i)
        ids.push_back(i);
      m_Prediction = m_Prediction->AddBundle(mitk::StreamlineFeatureExtractor::SubsetByIds(positive, ids));
    }

    m_PredictionNode = mitk::DataNode::New();
    m_PredictionNode->SetData(m_Prediction);
    m_PredictionNode->SetName("Prediction" + std::to_string(m_activeCycleCounter));
    this->GetDataStorage()->Add(m_PredictionNode);

    if (m_generalizationMode)
      m_Controls->m_predlabelingBrushGen->setEnabled(true);
    else
      m_Controls->m_predlabelingBrush->setEnabled(true);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateUncertaintySampleNode()
{
    if (!classifier)
      return;

     MITK_INFO << "Create Fibers to label based on Uncertainty";

     // next batch of the most uncertain streamlines that are still in the pool
     const auto& uncertain = classifier->GetUncertainIds();
     const bool generalization = m_generalizationMode != 0;
     const auto numRequested = static_cast<std::size_t>((generalization ? m_Controls->m_NumRandomFibers2Ben : m_Controls->m_NumRandomFibers2)->value());
     std::vector<unsigned int> ids;
     while (static_cast<std::size_t>(m_uncCounter) < uncertain.size() && ids.size() < numRequested)
     {
       const auto id = uncertain[m_uncCounter++];
       if (classifier->IsInPool(id))
         ids.push_back(id);
     }
     if (ids.empty())
     {
       QMessageBox::information(nullptr, "Helper", "All uncertain fibers of this cycle were already added. Retrain the classifier to get new ones.");
       return;
     }

     this->GetDataStorage()->Remove(m_UncertaintyLabelNode);
     m_UncertaintyLabel = classifier->ExtractFibers(ids);
     classifier->RemoveFromPool(ids);

     m_UncertaintyLabelNode = CreateBundleNode("UncertaintyLabel" + std::to_string(m_activeCycleCounter), m_UncertaintyLabel);
     this->GetDataStorage()->Add(m_UncertaintyLabelNode);

     m_PredictionNode->SetVisibility(false);
     (generalization ? m_Controls->m_unclabelingBrushGen : m_Controls->m_unclabelingBrush)->setChecked(true);

     UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateDistanceSampleNode()
{
    if (!classifier)
      return;

     MITK_INFO << "Create Fibers to label based on Distance in Features-Space";
     if (this->HelpMessagesEnabled())
     {
       QMessageBox::information(nullptr, "Helper", "A new subset is generated to be annotated\n"
                                                   "Again, use alt, shift and cntrl to label positive or negative fibers or to unlabel wrong annotated fibers\n"
                                                   "Then retrain the classifier..."
                                                   );
     }

    const bool generalization = m_generalizationMode != 0;
    const float entropyThreshold = (generalization ? m_Controls->m_subsetfftGen : m_Controls->m_subsetfft)->value() * 0.01f;
    const unsigned int numRequested = (generalization ? m_Controls->m_NumtolabelGen : m_Controls->m_Numtolabel)->value();

    auto ids = classifier->GetDiverseUncertainIds(entropyThreshold, numRequested);
    if (ids.empty())
    {
      QMessageBox::information(nullptr, "Helper", "No fibers left to label. Retrain the classifier to get new ones.");
      return;
    }

    m_DistanceLabel = classifier->ExtractFibers(ids);
    classifier->RemoveFromPool(ids);

    m_DistanceLabelNode = CreateBundleNode("DistanceLabel" + std::to_string(m_activeCycleCounter), m_DistanceLabel);
    this->GetDataStorage()->Add(m_DistanceLabelNode);

    m_PredictionNode->SetVisibility(false);

    if (generalization)
      m_Controls->m_distlabelingBrushGen->setChecked(true);
    else
      m_Controls->m_distlabelingBrush->setChecked(true);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RemovefromUncertaintyBrush( bool checked )
{
    if (checked)
    {
        m_Controls->m_predlabelingBrush->setChecked(false);
        m_Controls->m_predlabelingBrushGen->setChecked(false);
        m_UncertaintyLabel->SetFiberColors(255, 255, 255);
        this->StartLabelingBrush(m_UncertaintyLabelNode, false);
    }
    else
    {
        this->StopLabelingBrush();
    }
    RenderingManager::GetInstance()->RequestUpdateAll();

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RemovefromDistanceBrush(bool checked)
{
    if (checked)
    {
        m_Controls->m_predlabelingBrush->setChecked(false);
        m_Controls->m_predlabelingBrushGen->setChecked(false);

        m_DistanceLabel->SetFiberColors(255, 255, 255);
        this->StartLabelingBrush(m_DistanceLabelNode, false);
    }
    else
    {
        this->StopLabelingBrush();
    }
    RenderingManager::GetInstance()->RequestUpdateAll();

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RemovefromPredictionBrush( bool checked )
{
    if (checked)
    {
        if (this->HelpMessagesEnabled())
        {
          QMessageBox::information(nullptr, "Helper", "You may now exclude false predicted streamlines from the prediction\n"
                                                      "Again, use alt, shift and cntrl to label positive or negative fibers or to unlabel wrong annotated fibers\n"
                                                      "Then save the prediction..."
                                                      );
        }
        m_Controls->m_unclabelingBrush->setChecked(false);
        m_Controls->m_unclabelingBrushGen->setChecked(false);
        this->StartLabelingBrush(m_PredictionNode, true);
    }
    else
    {
        this->StopLabelingBrush();
    }

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateUncertantyMap()
{
    if (!classifier)
      return;

    const auto& entropies = classifier->GetEntropies();
    auto uncertaintymap = classifier->ExtractFibers(classifier->GetCandidateIds());

    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    weights->SetNumberOfValues(entropies.size());
    for (std::size_t i = 0; i < entropies.size(); ++i)
      weights->SetValue(i, entropies[i]);

    uncertaintymap->SetFiberWeights(weights);
    uncertaintymap->ColorFibersByFiberWeights(false, mitk::LookupTable::JET_TRANSPARENT);

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(uncertaintymap);
    node->SetName("uncertaintymap");
    this->GetDataStorage()->Add(node);
}

void QmitkInteractiveFiberDissectionView::ResetClassifier()
{
    m_IterationCounter = 0;
    m_RandomExtractionCounter = 0;
    m_activeCycleCounter = 0;
    m_createdStreamlineCounter = 0;
    m_prototypecounter = 0;
    m_uncCounter = 0;

    classifier.reset();
    this->GetDataStorage()->Remove(m_positiveBundleNode);
    this->GetDataStorage()->Remove(m_negativeBundleNode);
    this->GetDataStorage()->Remove(m_newfibersBundleNode);

    m_positiveBundleNode = nullptr;
    m_negativeBundleNode = nullptr;
    m_newfibersBundleNode = nullptr;
    m_positiveBundle = nullptr;
    m_negativeBundle = nullptr;
    m_newfibersBundle = nullptr;
    m_Controls->m_TrainClassifier->setEnabled(false);
    m_Controls->m_CreateUncertantyMap->setEnabled(false);
    m_Controls->m_AddDistanceFibers->setEnabled(false);
    m_Controls->m_distlabelingBrush->setChecked(false);
    m_Controls->m_unclabelingBrush->setChecked(false);
    m_Controls->m_predlabelingBrush->setChecked(false);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::addGeneralizationBundles(bool reset)
{
    // Combines the labels of the current subject with the generalization bundles. With reset, the combined bundles become the
    // new generalization bundles and the classifier is reset for the next subject. Otherwise they are only used for training.
    auto merge = [this, reset](mitk::DataNode::Pointer &genNode, const mitk::DataNode *labelNode, const std::string &name, float r, float g, float b)
    {
      auto merged = MergeBundles(GetBundle(labelNode), GetBundle(genNode));
      if (reset)
      {
        this->GetDataStorage()->Remove(genNode);
        merged->SetFiberColors(r, g, b);
        genNode = CreateBundleNode(name, merged);
        this->GetDataStorage()->Add(genNode);
      }
      if (genNode.IsNotNull())
        genNode->SetVisibility(false);
      return merged;
    };

    cur_negativeBundle = merge(m_negativeGenBundleNode, m_negativeBundleNode, "-BundleGen", 255, 0, 0);
    cur_positiveBundle = merge(m_positiveGenBundleNode, m_positiveBundleNode, "+BundleGen", 0, 255, 0);
    m_negativeGenBundle = GetBundle(m_negativeGenBundleNode);
    m_positiveGenBundle = GetBundle(m_positiveGenBundleNode);

    if (m_PredictionNode.IsNotNull())
      m_PredictionNode->SetVisibility(false);

    if (reset)
    {
        MITK_INFO << "Classifier reset";
        this->ResetClassifier();
    }

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::StartGeneralizationAlgorithm()
{
    m_generalizationMode = 1;
    m_prototypecounter = 0;
    if (m_Controls->m_helpmessages_2->isChecked())
    {
      QMessageBox::information(nullptr, "Helper", "Classifier is started and initiliazed with previous annotations...\n"
                                                  "Once it is finished, the prediction will be displayed.\n"
                                                  "If it is sufficient, you may save it. Feel free to exclude outliers, or to label more fibers to improve the prediction.\n"
                                                  "You can turn of this helper by unchecking the radie button in the lower part of the view."
                                                  );
    }

    cur_negativeBundle = GetBundle(m_negativeGenBundleNode);
    cur_positiveBundle = GetBundle(m_positiveGenBundleNode);
    m_testnode->SetVisibility(false);

    if (!this->RunClassifier(cur_positiveBundle, cur_negativeBundle))
      return;

    // empty label bundles for the annotations of this subject
    this->GetDataStorage()->Remove(m_negativeBundleNode);
    this->GetDataStorage()->Remove(m_positiveBundleNode);

    m_negativeBundle = mitk::FiberBundle::New();
    m_negativeBundleNode = CreateBundleNode("-Bundle", m_negativeBundle);
    this->GetDataStorage()->Add(m_negativeBundleNode);
    m_negativeBundleNode->SetVisibility(false);

    m_positiveBundle = mitk::FiberBundle::New();
    m_positiveBundleNode = CreateBundleNode("+Bundle", m_positiveBundle);
    this->GetDataStorage()->Add(m_positiveBundleNode);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::StartAlgorithmGen( )
{
    m_generalizationMode = 1;
    this->addGeneralizationBundles(false);

    if (m_Controls->m_helpmessages_2->isChecked())
    {
      QMessageBox::information(nullptr, "Helper", "Classifier is started and is trained on data from this and previous subjects...\n"
                                                  "Once it is finished, the prediction will be displayed.\n"
                                                  "If it is sufficient, you may save it. Feel free to exclude outliers, or to label more fibers to improve the prediction.\n"
                                                  "You can turn of this helper by unchecking the radie button in the lower part of the view."
                                                  );
    }

    m_testnode->SetVisibility(false);
    this->RunClassifier(cur_positiveBundle, cur_negativeBundle);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::AddSavedGeneralizationBundle(const mitk::DataNode *selectedNode, mitk::DataNode::Pointer &genNode, const std::string &name, float r, float g, float b)
{
    auto tractogram = GetBundle(selectedNode);
    if (tractogram.IsNull())
    {
        MITK_ERROR << "Selected node does not contain a FiberBundle.";
        return;
    }

    auto merged = MergeBundles(tractogram, GetBundle(genNode));
    merged->SetFiberColors(r, g, b);

    this->GetDataStorage()->Remove(genNode);
    genNode = CreateBundleNode(name, merged);
    this->GetDataStorage()->Add(genNode);
    genNode->SetVisibility(false);

    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::addGeneralizationSavedPosBundles()
{
    this->AddSavedGeneralizationBundle(m_Controls->m_posGenBundleBox->GetSelectedNode(), m_positiveGenBundleNode, "+BundleGen", 0, 255, 0);
    m_positiveGenBundle = GetBundle(m_positiveGenBundleNode);
}

void QmitkInteractiveFiberDissectionView::addGeneralizationSavedNegBundles()
{
    this->AddSavedGeneralizationBundle(m_Controls->m_negGenBundleBox->GetSelectedNode(), m_negativeGenBundleNode, "-BundleGen", 255, 0, 0);
    m_negativeGenBundle = GetBundle(m_negativeGenBundleNode);
}
