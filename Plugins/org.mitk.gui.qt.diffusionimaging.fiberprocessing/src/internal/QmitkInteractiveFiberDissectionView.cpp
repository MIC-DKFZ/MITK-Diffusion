﻿/*===================================================================

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

#include <QMessageBox>

#include <mitkNodePredicateProperty.h>
#include <mitkImageCast.h>
#include <mitkPointSet.h>
#include <mitkImageAccessByItk.h>
#include <mitkDataNodeObject.h>
#include <mitkTensorImage.h>
#include "mitkNodePredicateDataType.h"
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkPlanarCircle.h>
#include <mitkSurface.h>
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
#include <random>    // for std::default_random_engine


#include <mitkBaseRenderer.h>
#include <mitkRenderingManager.h>
#include <mitkVtkPropRenderer.h>
#include <mitkInteractionConst.h>
#include <mitkInteractionPositionEvent.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>

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
  , m_StreamlineInteractor(nullptr)
{

}

// Destructor
QmitkInteractiveFiberDissectionView::~QmitkInteractiveFiberDissectionView()
{
  //disable interactor
  if (m_StreamlineInteractor != nullptr)
  {
//    m_StreamlineInteractor->SetStreamlineNode(nullptr);
    m_StreamlineInteractor->EnableInteraction(false);
  }
}


void QmitkInteractiveFiberDissectionView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkInteractiveFiberDissectionViewControls;
    m_Controls->setupUi( parent );

    m_Controls->m_BundleBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isBundle= mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_BundleBox->SetPredicate( isBundle );

    m_Controls->m_PrototypeBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isPrototype = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_PrototypeBox->SetPredicate( isPrototype );

    m_Controls->m_PredictionBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isPrediction = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_PredictionBox->SetPredicate( isPrediction );

    m_Controls->m_GroundtruthBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isGroundtruth = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_GroundtruthBox->SetPredicate( isGroundtruth );

    m_Controls->m_TestBundleBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isTestBundle = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_TestBundleBox->SetPredicate( isTestBundle);




    connect(m_Controls->m_ErazorButton, SIGNAL(toggled(bool)), this, SLOT( RemovefromBundle(bool) ) );

    connect(m_Controls->m_BrushButton, SIGNAL(toggled(bool)), this, SLOT( RemovefromBundleBrush(bool) ) );

    connect(m_Controls->m_StreamlineCreation, SIGNAL( clicked() ), this, SLOT( CreateStreamline()));

    connect(m_Controls->m_AddRandomFibers, SIGNAL( clicked() ), this, SLOT( ExtractRandomFibersFromTractogram() ) ); //need

    connect(m_Controls->m_TrainClassifier, SIGNAL( clicked() ), this, SLOT( StartAlgorithm( )));

    connect(m_Controls->m_CreatePrediction, SIGNAL( clicked() ), this, SLOT( CreatePredictionNode( )));

    connect(m_Controls->m_certainData, SIGNAL( clicked() ), this, SLOT( CreateCertainNode( )));

    connect(m_Controls->m_AddUncertainFibers, SIGNAL( clicked() ), this, SLOT( CreateUncertaintySampleNode( )));

    connect(m_Controls->m_AddDistanceFibers, SIGNAL( clicked() ), this, SLOT( CreateDistanceSampleNode( )));

    connect(m_Controls->m_unclabeling, SIGNAL(toggled(bool)), this, SLOT( RemovefromUncertainty(bool) ) ); //need

    connect(m_Controls->m_unclabelingBrush, SIGNAL(toggled(bool)), this, SLOT( RemovefromUncertaintyBrush(bool) ) ); //need

    connect(m_Controls->m_distlabeling, SIGNAL(toggled(bool)), this, SLOT( RemovefromDistance(bool) ) ); //need

    connect(m_Controls->m_predlabeling, SIGNAL(toggled(bool)), this, SLOT( RemovefromPrediction(bool) ) ); //need

    connect(m_Controls->m_predlabelingBrush, SIGNAL(toggled(bool)), this, SLOT( RemovefromPredictionBrush(bool) ) ); //need

    connect(m_Controls->m_sellabeling, SIGNAL(toggled(bool)), this, SLOT( RemovefromSelection(bool) ) ); //need

    connect(m_Controls->m_ResampleButton, SIGNAL( clicked() ), this, SLOT( ResampleTractogram( ) ) );

    connect(m_Controls->m_RandomPrototypesButton, SIGNAL( clicked() ), this, SLOT( RandomPrototypes( ) ) );

    connect(m_Controls->m_TestBundleBox, SIGNAL( currentIndexChanged(int) ), this, SLOT( SFFPrototypes( ) ) );

    connect(m_Controls->m_validate, SIGNAL( clicked() ), this, SLOT( StartValidation( ) ) );

    connect(m_Controls->m_automaticLabelling, SIGNAL( clicked() ), this, SLOT( AutomaticLabelling( ) ) );

    connect(m_Controls->m_RemoveCertainData, SIGNAL( clicked() ), this, SLOT( RemoveCertainData( ) ) );

    connect(m_Controls->m_resetClassifier, SIGNAL( clicked() ), this, SLOT( ResetClassifier( ) ) );




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

   // disable alle frames

  m_Controls->m_ErazorButton->setCheckable(true);
  m_Controls->m_ErazorButton->setEnabled(true);
  m_Controls->m_BrushButton->setCheckable(true);
  m_Controls->m_BrushButton->setEnabled(true);
  m_Controls->m_unclabeling->setCheckable(true);
  m_Controls->m_unclabeling->setEnabled(false);
  m_Controls->m_predlabeling->setCheckable(true);
  m_Controls->m_predlabeling->setEnabled(false);
  m_Controls->m_unclabelingBrush->setCheckable(true);
  m_Controls->m_unclabelingBrush->setEnabled(false);
  m_Controls->m_predlabelingBrush->setCheckable(true);
  m_Controls->m_predlabelingBrush->setEnabled(false);
  m_Controls->m_distlabeling->setCheckable(true);
  m_Controls->m_distlabeling->setEnabled(false);
  m_Controls->m_sellabeling->setCheckable(true);
  m_Controls->m_sellabeling->setEnabled(false);




  m_Controls->m_StreamlineCreation->setEnabled(true);
  m_Controls->m_TrainClassifier->setEnabled(false);
  m_Controls->m_CreatePrediction->setEnabled(false);
  m_Controls->m_CreateUncertantyMap->setEnabled(false);
  m_Controls->m_Numtolabel->setEnabled(false);
  m_Controls->m_Numtolabel2->setEnabled(false);

  m_Controls->m_AddRandomFibers->setEnabled(false);
  m_Controls->m_AddDistanceFibers->setEnabled(false);
  m_Controls->m_AddUncertainFibers->setEnabled(false);
//  m_Controls->m_PrototypeBox->setEditable(false);
//  m_Controls->m_useStandardP->

  mitk::DataNode::Pointer curtestnode =  m_Controls->m_TestBundleBox->GetSelectedNode();
  bool testnodeSelected = curtestnode.IsNotNull();
  bool prototypesGenerated = m_Prototypes.IsNotNull();
  bool fibSelected = !m_SelectedFB.empty();
  bool multipleFibsSelected = (m_SelectedFB.size()>1);
  bool nfibSelected = m_newfibersBundleNode.IsNotNull();
  bool posSelected = this->GetDataStorage()->Exists(m_positiveBundleNode);
  bool negSelected = this->GetDataStorage()->Exists(m_negativeBundleNode);
  bool indexSelected = !m_index.empty();
  bool uncertaintySelected = this->GetDataStorage()->Exists(m_UncertaintyLabelNode);
  bool predictionSelected = this->GetDataStorage()->Exists(m_PredictionNode);



  // toggle visibility of elements according to selected method


  // are fiber bundles selected?
  if ( testnodeSelected )
  {
    m_Controls->m_AddRandomFibers->setEnabled(true);
    m_Controls->m_sellabeling->setEnabled(true);

    // if (m_Controls->m_TestBundleBox->OnPropertyListChanged()){
    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
    if  (!prototypesGenerated){
        this->SFFPrototypes();
    }
    else{
    }
  }



  if (fibSelected)
  {
    m_Controls->m_FibLabel->setText(QString(m_SelectedFB.at(0)->GetName().c_str()));
    // more than two bundles needed to join/subtract
    if (multipleFibsSelected )
    {
      m_Controls->m_FibLabel->setText("multiple bundles selected");
    }



  }


  // is image selected

  if (nfibSelected && posSelected)
  {
      m_Controls->m_ErazorButton->setEnabled(true);
      m_Controls->m_BrushButton->setEnabled(true);
  }


  if (posSelected && negSelected)
  {
      m_Controls->m_TrainClassifier->setEnabled(true);
  }

  if (indexSelected)
  {
      m_Controls->m_CreatePrediction->setEnabled(true);
      m_Controls->m_AddUncertainFibers->setEnabled(true);
      m_Controls->m_Numtolabel->setEnabled(true);
      m_Controls->m_AddDistanceFibers->setEnabled(true);
      m_Controls->m_Numtolabel2->setEnabled(true);
  }

  if (uncertaintySelected)
  {
      m_Controls->m_unclabeling->setEnabled(true);
      m_Controls->m_unclabelingBrush->setEnabled(true);
  }

  if (predictionSelected)
  {
      m_Controls->m_predlabeling->setEnabled(true);
      m_Controls->m_predlabelingBrush->setEnabled(true);
  }

//  if (distanceSelected)
//  {
      m_Controls->m_distlabeling->setEnabled(true);
//  }

//  if (m_Controls->m_useStandardP->isChecked())
//  {
//      m_Controls->m_PrototypeBox->setEditable(true);
//  }






}

void QmitkInteractiveFiberDissectionView::OnEndInteraction()
{

}

void QmitkInteractiveFiberDissectionView::ResampleTractogram()
{
    mitk::DataNode::Pointer node = m_Controls->m_BundleBox->GetSelectedNode();
    auto tractogram = dynamic_cast<mitk::FiberBundle *>(node->GetData());
    mitk::FiberBundle::Pointer tempfib = tractogram->GetDeepCopy();

    std::vector<int> myvec;
    for (unsigned int k=0; k<tempfib->GetNumFibers(); k++)
    {
      myvec.push_back(k);
    }
//      auto rng = std::default_random_engine {};
    std::random_shuffle(std::begin(myvec), std::end(myvec));

    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

     /* Check wether all Streamlines of the bundles are labeled... If all are labeled Skip for Loop*/
    unsigned int counter = 0;

    for (unsigned int i=0; i<tempfib->GetNumFibers(); i++)
    {
      vtkCell* cell = tempfib->GetFiberPolyData()->GetCell(myvec.at(i));
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; j++)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      weights->InsertValue(counter, tempfib->GetFiberWeight(myvec.at(i)));
      vNewLines->InsertNextCell(container);
      counter++;

    }

    vNewPolyData->SetLines(vNewLines);
    vNewPolyData->SetPoints(vNewPoints);

    mitk::FiberBundle::Pointer ShuffledBundle = mitk::FiberBundle::New(vNewPolyData);
    ShuffledBundle->SetFiberWeights(weights);

    ShuffledBundle->ResampleToNumPoints(40);
    MITK_INFO << "Resampling Done";

    mitk::DataNode::Pointer newnode = mitk::DataNode::New();
    newnode->SetData( ShuffledBundle );
    newnode->SetName( node->GetName() + "_" + std::to_string(40) );
    this->GetDataStorage()->Add(newnode);
    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RandomPrototypes()
{
    // Get the number of fibers to use as prototypes from the user interface
    const int numPrototypes = m_Controls->m_NumPrototypes->value();
    MITK_INFO << "Number of Fibers to use as Prototypes: " << numPrototypes;

    // Get the bundle of fibers to extract prototypes from
    mitk::DataNode::Pointer selectedNode = m_Controls->m_BundleBox->GetSelectedNode();
    if (!selectedNode)
    {
        MITK_ERROR << "No fiber bundle selected";
        return;
    }
    mitk::FiberBundle::Pointer fiberBundle = dynamic_cast<mitk::FiberBundle*>(selectedNode->GetData());
    if (!fiberBundle)
    {
        MITK_ERROR << "Selected node does not contain a fiber bundle";
        return;
    }
    MITK_INFO << "Number of fibers in the selected bundle: " << fiberBundle->GetNumFibers();

    // Randomly shuffle the indices of fibers in the bundle
    std::vector<int> indices(fiberBundle->GetNumFibers());
    std::iota(std::begin(indices), std::end(indices), 0); // fill the vector with 0 to n-1
    std::random_shuffle(std::begin(indices), std::end(indices));

    // Create new PolyData to hold the prototype fibers
    vtkSmartPointer<vtkPolyData> prototypePolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> prototypeLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> prototypePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> prototypeWeights = vtkSmartPointer<vtkFloatArray>::New();

    // Extract a subset of fibers as prototypes
    for (int i = 0; i < numPrototypes; i++)
    {
        if (indices.empty()) break;
        const int index = indices.back();
        indices.pop_back();

        vtkCell* cell = fiberBundle->GetFiberPolyData()->GetCell(index);
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
        for (unsigned int j = 0; j < numPoints; j++)
        {
            double p[3];
            points->GetPoint(j, p);
            vtkIdType id = prototypePoints->InsertNextPoint(p);
            container->GetPointIds()->InsertNextId(id);
        }

        prototypeWeights->InsertNextValue(fiberBundle->GetFiberWeight(index));
        prototypeLines->InsertNextCell(container);
    }

    // Create a new FiberBundle containing the prototype fibers
    prototypePolyData->SetLines(prototypeLines);
    prototypePolyData->SetPoints(prototypePoints);
    mitk::FiberBundle::Pointer prototypesBundle = mitk::FiberBundle::New(prototypePolyData);
    prototypesBundle->SetFiberWeights(prototypeWeights);
    auto node = mitk::DataNode::New();
    node->SetData(prototypesBundle);
    node->SetName("Random_Prototypes");
    node->SetVisibility(false);

    GetDataStorage()->Add(node);
    m_Controls->m_PrototypeBox->SetAutoSelectNewItems(true);
    m_Controls->m_PrototypeBox->SetSelectedNode(node);
    m_Controls->m_PrototypeBox->SetAutoSelectNewItems(false);
    m_Controls->m_useStandardP->setChecked(false);

    MITK_INFO << "Created new node with " << numPrototypes << " random prototypes";

}

void QmitkInteractiveFiberDissectionView::SFFPrototypes()
{   
    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
    if (m_Prototypes.IsNotNull()){
        this->GetDataStorage()->Remove(m_Prototypes);
    }

     MITK_INFO << "Number of Fibers to use as Prototypes: ";
     MITK_INFO << m_Controls->m_NumPrototypes->value();
     MITK_INFO << "Start Creating Prototypes based on SFF";

      mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());

      /* Get Subset of Tractogram*/
      int size_subset = std::max(1.0, ceil(3.0 * m_Controls->m_NumPrototypes->value() * std::log(m_Controls->m_NumPrototypes->value())));

      MITK_INFO << fib->GetNumFibers();
      std::vector<int> myvec;
      for (unsigned int k=0; k<fib->GetNumFibers(); k++)
      {
        myvec.push_back(k);
      }

//      std::random_shuffle(std::begin(myvec), std::end(myvec));

      vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
      vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

      vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();


      unsigned int counter = 0;

      for (int i=0; i<size_subset; i++)
      {
        vtkCell* cell = fib->GetFiberPolyData()->GetCell(myvec.at(i));
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
        for (unsigned int j=0; j<numPoints; j++)
        {
          double p[3];
          points->GetPoint(j, p);

          vtkIdType id = vNewPoints->InsertNextPoint(p);
          container->GetPointIds()->InsertNextId(id);
        }
        weights->InsertValue(counter, fib->GetFiberWeight(myvec.at(i)));
        vNewLines->InsertNextCell(container);
        counter++;

      }


      vNewPolyData->SetLines(vNewLines);
      vNewPolyData->SetPoints(vNewPoints);

      mitk::FiberBundle::Pointer temp_fib = mitk::FiberBundle::New(vNewPolyData);
      temp_fib->SetFiberWeights(weights);

      MITK_INFO << temp_fib->GetFiberPolyData()->GetNumberOfCells();
      /* Create std::vector <vnl_matrix> of the SubsetBundle*/
      std::vector< vnl_matrix<float> > out_fib(temp_fib->GetFiberPolyData()->GetNumberOfCells());

      for (int i=0; i<temp_fib->GetFiberPolyData()->GetNumberOfCells(); i++)
      {
        vtkCell* cell = temp_fib->GetFiberPolyData()->GetCell(i);
        int numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        vnl_matrix<float> streamline;
        streamline.set_size(3, cell->GetNumberOfPoints());
        streamline.fill(0.0);

        for (int j=0; j<numPoints; j++)
        {
          double cand[3];
          points->GetPoint(j, cand);

          vnl_vector_fixed< float, 3 > candV;
          candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
          streamline.set_column(j, candV);
        }

    //        out_fib.push_back(streamline);
        out_fib.at(i)=streamline;
      }


      /* Calculate the distancematrix of Subset*/
      std::vector< vnl_matrix<float> >  dist_vec(out_fib.size());//
      cv::parallel_for_(cv::Range(0, out_fib.size()), [&](const cv::Range &range)
      {
      for (int i = range.start; i < range.end; i++)

  //    for (unsigned int i=0; i<tractogram.size(); i++)
      {

          vnl_matrix<float> distances;
          distances.set_size(1, out_fib.size());
          distances.fill(0.0);


          for (unsigned int j=0; j<out_fib.size(); j++)
          {
              vnl_matrix<float> single_distances;
              single_distances.set_size(1, out_fib.at(0).cols());
              single_distances.fill(0.0);
              vnl_matrix<float> single_distances_flip;
              single_distances_flip.set_size(1, out_fib.at(0).cols());
              single_distances_flip.fill(0.0);
              for (unsigned int ik=0; ik<out_fib.at(0).cols(); ik++)
              {
                  double cur_dist;
                  double cur_dist_flip;

                  cur_dist = sqrt(pow(out_fib.at(i).get(0,ik) - out_fib.at(j).get(0,ik), 2.0) +
                                         pow(out_fib.at(i).get(1,ik) - out_fib.at(j).get(1,ik), 2.0) +
                                         pow(out_fib.at(i).get(2,ik) - out_fib.at(j).get(2,ik), 2.0));

                  cur_dist_flip = sqrt(pow(out_fib.at(i).get(0,ik) - out_fib.at(j).get(0,out_fib.at(0).cols()-(ik+1)), 2.0) +
                                         pow(out_fib.at(i).get(1,ik) - out_fib.at(j).get(1,out_fib.at(0).cols()-(ik+1)), 2.0) +
                                         pow(out_fib.at(i).get(2,ik) - out_fib.at(j).get(2,out_fib.at(0).cols()-(ik+1)), 2.0));

                  single_distances.put(0,ik, cur_dist);
                  single_distances_flip.put(0,ik, cur_dist_flip);

              }

              if (single_distances_flip.mean()> single_distances.mean())
              {
                  distances.put(0,j, single_distances.mean());
              }
              else {
                  distances.put(0,j, single_distances_flip.mean());
              }



          }
  //        dist_vec.push_back(distances);
          dist_vec.at(i) = distances;
      }
      });


      /*Index to find values in distancematrix*/
      std::vector<unsigned int> myidx;
      /*Index to find actual streamlines using indexUnc*/
      std::vector<unsigned int> indexUncDist;
      /*Start with the Streamline of the highest entropy, which is in distance_matrix at idx 0*/
      myidx.push_back(0);

      /*Vecotr that stores minvalues of current iteration*/
      vnl_matrix<float> cur_vec;
      cur_vec.set_size(1, size_subset);
      cur_vec.fill(0.0);
      for (int i=0; i<m_Controls->m_NumPrototypes->value(); i++)
      {

  //        unsigned int cur_i = indexUnc.at(myidx.at(i));

          /*Save mean distance of all used Samples*/

          vnl_matrix<float> sum_matrix;
          sum_matrix.set_size(myidx.size(), size_subset);
          sum_matrix.fill(0);
          for (unsigned int ii=0; ii<myidx.size(); ii++)
          {

              sum_matrix.set_row(ii, dist_vec.at(myidx.at(ii)).get_row(0));
          }

          for (unsigned int k=0; k<sum_matrix.columns(); k++)
          {
              cur_vec.put(0,k, sum_matrix.get_column(k).min_value());
          }
          myidx.push_back(cur_vec.arg_max());
          sum_matrix.clear();

        }
          vtkSmartPointer<vtkPolyData> vNewPolyData2 = vtkSmartPointer<vtkPolyData>::New();
          vtkSmartPointer<vtkCellArray> vNewLines2 = vtkSmartPointer<vtkCellArray>::New();
          vtkSmartPointer<vtkPoints> vNewPoints2 = vtkSmartPointer<vtkPoints>::New();

          vtkSmartPointer<vtkFloatArray> weights2 = vtkSmartPointer<vtkFloatArray>::New();

           /* Check wether all Streamlines of the bundles are labeled... If all are labeled Skip for Loop*/
          counter = 0;

          for (int i=0; i<m_Controls->m_NumPrototypes->value(); i++)
          {

            vtkCell* cell = fib->GetFiberPolyData()->GetCell(myidx.at(i));
            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
            for (unsigned int j=0; j<numPoints; j++)
            {
              double p[3];
              points->GetPoint(j, p);

              vtkIdType id = vNewPoints2->InsertNextPoint(p);
              container->GetPointIds()->InsertNextId(id);
            }
            weights2->InsertValue(counter, fib->GetFiberWeight(myvec.at(i)));
            vNewLines2->InsertNextCell(container);
            counter++;

          }

          vNewPolyData2->SetLines(vNewLines2);
          vNewPolyData2->SetPoints(vNewPoints2);

          mitk::FiberBundle::Pointer PrototypesBundle = mitk::FiberBundle::New(vNewPolyData2);
          PrototypesBundle->SetFiberWeights(weights2);
          MITK_INFO << "Number of cells: " << PrototypesBundle->GetFiberPolyData()->GetNumberOfCells();



      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(PrototypesBundle);
      node->SetName("Prototypes");

      m_Prototypes = node;
////      MITK_INFO << "Number of Streamlines in first function";
////      MITK_INFO << m_newfibersBundleNode->GetData()->GetFiberPolyData()->GetNumberOfCells();
      m_Controls->m_PrototypeBox->SetAutoSelectNewItems (true);
      this->GetDataStorage()->Add(m_Prototypes, m_testnode);
      m_Prototypes->SetVisibility(false);
      m_Controls->m_PrototypeBox->SetAutoSelectNewItems (false);
      m_Controls->m_useStandardP->setChecked(false);



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




void QmitkInteractiveFiberDissectionView::CreateStreamline()
{
    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();

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


    m_reducedFibersDataNode = mitk::DataNode::New();
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    mitk::FiberBundle::Pointer reducedFib = mitk::FiberBundle::New();
    m_reducedFibersDataNode->SetData(reducedFib);
    m_reducedFibersDataNode->SetName("Reduced");
//    m_reducedFibersDataNode = node;

    this->GetDataStorage()->Add(m_reducedFibersDataNode, m_testnode);

    MITK_INFO << "Get Bundle";
//    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());

    MITK_INFO <<"Sphere to interactor";
    m_SphereInteractor->workingBundleNode(m_testnode, m_reducedFibersDataNode);
    m_SphereInteractor->StartEndNodes(startDataNode, endDataNode);




    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::ExtractRandomFibersFromTractogram()
{
    m_reducedFibersDataNode->SetVisibility(false);
    // Hide the selected node.
    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
//    m_testnode = m_SphereInteractor->m_reducedFibersBundle;
    m_testnode->SetVisibility(false);

    // Uncheck the Brush and Erazor buttons.
    m_Controls->m_ErazorButton->setChecked(false);
    m_Controls->m_BrushButton->setChecked(false);


    // Print the number of fibers to extract.
    MITK_INFO << "Number of Fibers to extract from Tractogram: ";
    MITK_INFO << m_Controls->m_NumRandomFibers->value();


    // If a newfibersBundleNode exists, remove it from the data storage.
//    if (!m_newfibersBundleNode)
//    {
//        mitk::FiberBundle::Pointer m_newfibersBundle = mitk::FiberBundle::New();
//        mitk::DataNode::Pointer node = mitk::DataNode::New();
//        node->SetName("ToLabel");
//        node->SetData(m_newfibersBundle);
//        m_newfibersBundleNode = node;
//        this->GetDataStorage()->Add(m_newfibersBundleNode);
//    }
    if (m_newfibersBundleNode)
    {
        MITK_INFO << "To Label Bundle Exists";
        this->GetDataStorage()->Remove(m_newfibersBundleNode);
//        mitk::FiberBundle::Pointer Stack = dynamic_cast<mitk::FiberBundle *>(m_newfibersBundleNode->GetData());

        // Create a new data node and set it to an empty FiberBundle.
        mitk::DataNode::Pointer node = mitk::DataNode::New();
        m_newfibersFibersData = vtkSmartPointer<vtkPolyData>::New();
        m_newfibersFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
        m_newfibersBundle = mitk::FiberBundle:: New(m_newfibersFibersData);
        m_newfibersFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());

        // Print a message to indicate the creation of the FiberBundle.
        MITK_INFO << "Create Bundle";
    }
    // Get the selected node, which is assumed to be a FiberBundle.
    MITK_INFO << "Get Test Data";
//    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());

    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_reducedFibersDataNode->GetData());


    MITK_INFO << "Loaded";
    // Create a new PolyData object, and its associated points and cell arrays.
    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

    // Initialize a counter variable.
    unsigned int counter = 0;

    if (m_Controls->m_NumRandomFibers->value() != fib->GetFiberPolyData()->GetNumberOfCells()) {

        // Create a vector of indices
        std::vector<int> indices(fib->GetFiberPolyData()->GetNumberOfCells());
        std::iota(indices.begin(), indices.end(), 0); // fill the vector with 0, 1, 2, ..., numFibers-1
        std::shuffle(indices.begin(), indices.end(), std::default_random_engine(std::random_device()())); // shuffle the indices
        MITK_INFO << "Shuffling";

        // Iterate through the fibers in a random order
        for (int i = 0; i < m_Controls->m_NumRandomFibers->value() / 2; i++) {
            int index = indices[i];

            // Get the index-th fiber in the input dataset
            vtkCell* cell = fib->GetFiberPolyData()->GetCell(index);

            // Get the number of points in the fiber
            auto numPoints = cell->GetNumberOfPoints();

            // Get the points in the fiber
            vtkPoints* points = cell->GetPoints();

            // Create a new vtkPolyLine container to store the new fiber
            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

            // Loop through each point in the fiber
            for (unsigned int j = 0; j < numPoints; j++) {
                double p[3];
                points->GetPoint(j, p);

                // Insert the new point into the vtkPoints object and get its ID
                vtkIdType id = vNewPoints->InsertNextPoint(p);
                // Add the ID to the container object to create a new fiber
                container->GetPointIds()->InsertNextId(id);
            }

            // Insert the fiber weight into the vtkDoubleArray object
            //weights->InsertValue(counter, fib->GetFiberWeight(index));

            vNewLines->InsertNextCell(container);

            // Increment the counter variable
            counter++;
        }
    }

        fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());
        if (m_Controls->m_NumRandomFibers->value()/2 != fib->GetFiberPolyData()->GetNumberOfCells()) {

            // Create a vector of indices
            std::vector<int> indices(fib->GetFiberPolyData()->GetNumberOfCells());
            std::iota(indices.begin(), indices.end(), 0); // fill the vector with 0, 1, 2, ..., numFibers-1
            std::shuffle(indices.begin(), indices.end(), std::default_random_engine(std::random_device()())); // shuffle the indices

            // Iterate through the fibers in a random order
            for (int i = 0; i < m_Controls->m_NumRandomFibers->value() / 2; i++) {
                int index = indices[i];

                // Get the index-th fiber in the input dataset
                vtkCell* cell = fib->GetFiberPolyData()->GetCell(index);

                // Get the number of points in the fiber
                auto numPoints = cell->GetNumberOfPoints();

                // Get the points in the fiber
                vtkPoints* points = cell->GetPoints();

                // Create a new vtkPolyLine container to store the new fiber
                vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

                // Loop through each point in the fiber
                for (unsigned int j = 0; j < numPoints; j++) {
                    double p[3];
                    points->GetPoint(j, p);

                    // Insert the new point into the vtkPoints object and get its ID
                    vtkIdType id = vNewPoints->InsertNextPoint(p);
                    // Add the ID to the container object to create a new fiber
                    container->GetPointIds()->InsertNextId(id);
                }

                // Insert the fiber weight into the vtkDoubleArray object
                //weights->InsertValue(counter, fib->GetFiberWeight(index));

                vNewLines->InsertNextCell(container);

                // Increment the counter variable
                counter++;
            }
        }




    // Delete the old fibers from the input dataset
//    for (int i = 0; i < m_Controls->m_NumRandomFibers->value(); i++) {
//        fib->GetFiberPolyData()->DeleteCell(i);
//    }
    // Remove any deleted cells from the input dataset
//    fib->GetFiberPolyData()->RemoveDeletedCells();

    // Output the number of fibers in the new dataset to the console
//    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();

    // Set the vtkPolyData object to the new fibers
    vNewPolyData->SetLines(vNewLines);
    vNewPolyData->SetPoints(vNewPoints);

    // Create a new vtkPolyData object and set it to the new fibers
    m_newfibersFibersData = vtkSmartPointer<vtkPolyData>::New();
    m_newfibersFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
    m_newfibersFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
    m_newfibersFibersData->SetPoints(vNewPoints);
    m_newfibersFibersData->SetLines(vNewLines);

    // Create a new mitk::FiberBundle object and set it to the new fibers
    m_newfibersBundle = mitk::FiberBundle::New(vNewPolyData);
    m_newfibersBundle->SetFiberColors(255, 255, 255);
//    m_newfibersBundle->SetFiberWeights(weights);

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_newfibersBundle);
    node->SetName("ToLabel");
    m_newfibersBundleNode = node;
//      m_newfibersBundleNode->SetData(m_newfibersBundle);
    this->GetDataStorage()->Add(m_newfibersBundleNode);

    // Create a new data node for the negativeBundle object
    if (!m_negativeBundleNode)
    {
        mitk::FiberBundle::Pointer m_negativeBundle = mitk::FiberBundle::New();
        mitk::DataNode::Pointer node2 = mitk::DataNode::New();
        node2->SetName("-Bundle");
        node2->SetData(m_negativeBundle);
        m_negativeBundleNode = node2;
        this->GetDataStorage()->Add(m_negativeBundleNode);
    }

    if (!m_positiveBundleNode)
    {
        mitk::FiberBundle::Pointer m_positiveBundle = mitk::FiberBundle::New();
        mitk::DataNode::Pointer node3 = mitk::DataNode::New();
        node3->SetName("+Bundle");
        node3->SetData(m_positiveBundle);
        m_positiveBundleNode = node3;
        this->GetDataStorage()->Add(m_positiveBundleNode);
    }


    // Increment the RandomExtractionCounter
    m_RandomExtractionCounter++;

    m_Controls->m_BrushButton->setChecked(true);




    UpdateGui();


}

void QmitkInteractiveFiberDissectionView::RemovefromBundle( bool checked )
{
    if (checked)
    {
        m_Controls->m_BrushButton->setChecked(false);

//        if (m_StreamlineInteractor.IsNull())
//        {
            this->CreateStreamlineInteractorBrush();
            this->CreateStreamlineInteractor();
//        }
        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->LabelfromPrediction(false);
        m_StreamlineInteractor->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractor->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractor->SetToLabelNode(m_newfibersBundleNode);
        m_StreamlineInteractor->EnableInteraction(true);

//            m_StreamlineInteractor->EnableInteraction(true);
//            m_StreamlineInteractor->LabelfromPrediction(false);
//            m_StreamlineInteractor->SetNegativeNode(m_negativeBundleNode);
//            m_StreamlineInteractor->SetPositiveNode(m_positiveBundleNode);
//            m_StreamlineInteractor->SetToLabelNode(m_newfibersBundleNode);
//        }
//        else
//        {
//            m_StreamlineInteractor->EnableInteraction(true);
//            m_StreamlineInteractor->LabelfromPrediction(false);
//            m_StreamlineInteractor->SetPositiveNode(m_positiveBundleNode);
//            m_StreamlineInteractor->SetToLabelNode(m_newfibersBundleNode);
//        }
    }
    else
    {
      m_StreamlineInteractor->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
    }


   UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RemovefromBundleBrush( bool checked )
{
    if (checked)
    {
        m_Controls->m_ErazorButton->setChecked(false);

//        if (m_StreamlineInteractorBrush.IsNull())
//        {
            this->CreateStreamlineInteractor();
            this->CreateStreamlineInteractorBrush();
//        }


        m_StreamlineInteractorBrush->EnableInteraction(true);
        m_StreamlineInteractorBrush->LabelfromPrediction(false);
        m_StreamlineInteractorBrush->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractorBrush->SetToLabelNode(m_newfibersBundleNode);
        m_StreamlineInteractorBrush->EnableInteraction(true);
//            m_StreamlineInteractorBrush->EnableInteraction(true);
//            m_StreamlineInteractorBrush->LabelfromPrediction(false);
//            m_StreamlineInteractorBrush->SetNegativeNode(m_negativeBundleNode);
//            m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
//            m_StreamlineInteractorBrush->SetToLabelNode(m_newfibersBundleNode);
//        }
//        else
//        {
//            m_StreamlineInteractorBrush->EnableInteraction(true);
//            m_StreamlineInteractorBrush->LabelfromPrediction(false);
//            m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
//            m_StreamlineInteractorBrush->SetToLabelNode(m_newfibersBundleNode);
//        }
    }
    else
    {
      m_StreamlineInteractorBrush->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
    }


   UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateStreamlineInteractor()
{

    m_StreamlineInteractor = mitk::StreamlineInteractor::New();

    m_StreamlineInteractor->LoadStateMachine("Streamline3DStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_StreamlineInteractor->SetEventConfig("Streamline3DConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));

//  m_StreamlineInteractor->SetRotationEnabled(rotationEnabled);
}

void QmitkInteractiveFiberDissectionView::CreateStreamlineInteractorBrush()
{

    m_StreamlineInteractorBrush = mitk::StreamlineInteractorBrush::New();

    m_StreamlineInteractorBrush->LoadStateMachine("StreamlineBrush3DStates.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));
    m_StreamlineInteractorBrush->SetEventConfig("StreamlineBrush3DConfig.xml", us::ModuleRegistry::GetModule("MitkFiberDissection"));

//  m_StreamlineInteractor->SetRotationEnabled(rotationEnabled);
}

void QmitkInteractiveFiberDissectionView::StartAlgorithm()
{

    m_negativeBundle = dynamic_cast<mitk::FiberBundle*>(m_negativeBundleNode->GetData());
    m_positiveBundle = dynamic_cast<mitk::FiberBundle*>(m_positiveBundleNode->GetData());


    this->GetDataStorage()->Remove(m_UncertaintyLabelNode);
    this->GetDataStorage()->Remove(m_DistanceLabelNode);


    MITK_INFO << "Clean Test Data";

    m_testnode = m_Controls->m_TestBundleBox->GetSelectedNode();
//     mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());

    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();
    fib->GetFiberPolyData()->RemoveDeletedCells();
    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();


    m_Controls->m_unclabeling->setChecked(false);
    m_Controls->m_distlabeling->setChecked(false);
    m_Controls->m_unclabelingBrush->setChecked(false);
    m_Controls->m_predlabeling->setChecked(false);
    m_Controls->m_predlabelingBrush->setChecked(false);
    m_uncCounter = 0;

//    classifier.reset();
    MITK_INFO << "Extract Features";

    classifier = std::make_shared<mitk::StreamlineFeatureExtractor>();
    classifier->SetActiveCycle(m_activeCycleCounter);
    classifier->SetTractogramPlus(m_positiveBundle);
    classifier->SetTractogramMinus(m_negativeBundle);
    classifier->SetTractogramPrototypes(dynamic_cast<mitk::FiberBundle*>(m_Prototypes->GetData()), m_Controls->m_useStandardP->isChecked());
//    classifier->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData()), m_SelectedFB.at(0)->GetName());
    classifier->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData()), m_testnode->GetName());
//    classifier->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_trainbundle->GetData()), m_trainbundle->GetName());



    classifier->Update();

    m_index = classifier->m_index;
    MITK_INFO << "Number of Cycles";
    MITK_INFO << m_activeCycleCounter;
    m_activeCycleCounter += 1;


    MITK_INFO << "Algorithm run succesfully";
//    mitk::DataNode::Pointer node = mitk::DataNode::New();
//    this->GetDataStorage()->Add(node);
    m_Controls->m_CreatePrediction->setEnabled(true);

    this->CreatePredictionNode();
    UpdateGui();

}

void QmitkInteractiveFiberDissectionView::CreatePredictionNode()
{
    MITK_INFO << "Create Prediction";

    m_Prediction = classifier->CreatePrediction(m_index.at(0), false);



//    m_Prediction->ColorFibersByFiberWeights(false, mitk::LookupTable::JET);
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_Prediction);
    auto s = std::to_string(m_activeCycleCounter);
    node->SetName("Prediction"+s);
    m_PredictionNode = node;
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_PredictionNode->GetData());
    mitk::FiberBundle::Pointer fib_true = dynamic_cast<mitk::FiberBundle*>(m_positiveBundleNode->GetData());
    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

     /* Check wether all Streamlines of the bundles are labeled... If all are labeled Skip for Loop*/
    unsigned int counter = 0;

    for (int i=0; i<fib->GetFiberPolyData()->GetNumberOfCells(); i++)
    {
      vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; j++)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      weights->InsertValue(counter, fib->GetFiberWeight(i));
      vNewLines->InsertNextCell(container);
      counter++;

    }
    for (int i=0; i<fib_true->GetFiberPolyData()->GetNumberOfCells(); i++)
    {
      vtkCell* cell = fib_true->GetFiberPolyData()->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; j++)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
//      weights->InsertValue(counter, fib_true->GetFiberWeight(i));
      vNewLines->InsertNextCell(container);
      counter++;

    }


    vNewPolyData->SetLines(vNewLines);
    vNewPolyData->SetPoints(vNewPoints);

    mitk::FiberBundle::Pointer m_Prediction = mitk::FiberBundle::New(vNewPolyData);
//    m_Prediction->SetFiberWeights(weights);

//    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_Prediction);
    node->SetName("Prediction"+s);
//    this->GetDataStorage()->Add(node);
    m_PredictionNode = node;
    this->GetDataStorage()->Add(m_PredictionNode);
    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::CreateCertainNode()
{
    MITK_INFO << "Create Certain Data";


    m_CertainMinus = classifier->CreatePrediction(m_index.at(2), false);
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_CertainMinus);
    auto s = std::to_string(m_activeCycleCounter);
    node->SetName("m_CertainMinus"+s);
    m_CertainMinusNode = node;
    this->GetDataStorage()->Add(m_CertainMinusNode);

//    m_CertainPlus = classifier->CreatePrediction(m_index.at(4));
//    mitk::DataNode::Pointer node2= mitk::DataNode::New();
//    node2->SetData(m_CertainPlus);
//    node2->SetName("m_CertainPlus"+s);
//    m_CertainPlusNode = node2;
//    this->GetDataStorage()->Add(m_CertainPlusNode);


//    m_CertainBetweenMinus = classifier->CreatePrediction(m_index.at(4));
//    mitk::DataNode::Pointer node3 = mitk::DataNode::New();
//    node3->SetData(m_CertainBetweenMinus);
//    node3->SetName("m_CertainBetweenMinus"+s);
//    m_CertainBetweenMinusNode = node3;
//    this->GetDataStorage()->Add(m_CertainBetweenMinusNode);

//    m_CertainBetweenPlus = classifier->CreatePrediction(m_index.at(6));
//    mitk::DataNode::Pointer node4= mitk::DataNode::New();
//    node4->SetData(m_CertainBetweenPlus);
//    node4->SetName("m_CertainBetweenPlus"+s);
//    m_CertainBetweenPlusNode = node4;
//    this->GetDataStorage()->Add(m_CertainBetweenPlusNode);

}

void QmitkInteractiveFiberDissectionView::RemoveCertainData()
{
//    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());
    MITK_INFO << "Length of certain negativ data:";
    MITK_INFO << m_index.at(2).size();
    MITK_INFO << "Length of Testdata";
    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();
    for (unsigned int i=0; i< m_index.at(2).size(); i++)
    {
        fib->GetFiberPolyData()->DeleteCell(m_index.at(2).at(i));
    }
//    fib->GetFiberPolyData()->RemoveDeletedCells();

    MITK_INFO << "Length of Testdata";
    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();
}

void QmitkInteractiveFiberDissectionView::CreateUncertaintySampleNode()
{
     this->GetDataStorage()->Remove(m_UncertaintyLabelNode);
     MITK_INFO << "Create Fibers to label based on Uncertainty";

     std::vector<unsigned int> vec = m_index.at(1);
     std::vector<unsigned int> myvec = {vec.begin() + m_uncCounter, vec.begin() + m_uncCounter + m_Controls->m_Numtolabel->value()};
     m_uncCounter = m_uncCounter + m_Controls->m_Numtolabel->value();
     MITK_INFO << m_index.at(1).size();
     MITK_INFO << myvec.size();


    m_UncertaintyLabel = classifier->CreatePrediction(myvec, true);
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_UncertaintyLabel);
    float d = 5.0;

    auto s = std::to_string(m_activeCycleCounter);
    node->SetName("UncertaintyLabel"+s);
    m_UncertaintyLabelNode = node;
    m_UncertaintyLabelNode->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(d));
    this->GetDataStorage()->Add(m_UncertaintyLabelNode);
    UpdateGui();
    m_PredictionNode->SetVisibility(false);
    m_Controls->m_unclabelingBrush->setChecked(true);


}

void QmitkInteractiveFiberDissectionView::CreateDistanceSampleNode()
{
//     MITK_INFO << "Create Fibers to label based on Distance in Features-Space";
//     std::vector<unsigned int> myvec = m_index.at(2);
//     myvec.resize(m_Controls->m_Numtolabel2->value());
//     MITK_INFO << m_index.at(2).size();
//     MITK_INFO << myvec.size();

    float myval = m_Controls->m_subsetfft->value() * 0.01;
    MITK_INFO << myval;
    std::vector<std::vector<unsigned int>> curidx;
    curidx =  classifier->GetDistanceData(myval);
    std::vector<unsigned int> myvec = curidx.at(0);
    myvec.resize(m_Controls->m_Numtolabel2->value());
    MITK_INFO << m_index.at(2).size();
    MITK_INFO << myvec.size();
    m_DistanceLabel = classifier->CreatePrediction(myvec, true);
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(m_DistanceLabel);
    auto s = std::to_string(m_activeCycleCounter);
    node->SetName("DistanceLabel"+s);
    m_DistanceLabelNode = node;
    this->GetDataStorage()->Add(m_DistanceLabelNode);
    UpdateGui();
}

void QmitkInteractiveFiberDissectionView::RemovefromUncertainty( bool checked )
{
//    m_PredictionNode->SetVisibility(false);
    if (checked)
    {
        m_Controls->m_unclabelingBrush->setChecked(false);
        m_Controls->m_predlabelingBrush->setChecked(false);
        m_Controls->m_predlabeling->setChecked(false);


        m_UncertaintyLabel->SetFiberColors(255, 255, 255);
        this->CreateStreamlineInteractor();
        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->LabelfromPrediction(false);
        m_StreamlineInteractor->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractor->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractor->SetToLabelNode(m_UncertaintyLabelNode);
    }
    else
    {
      m_StreamlineInteractor->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkInteractiveFiberDissectionView::RemovefromUncertaintyBrush( bool checked )
{
    if (checked)
    {
        m_Controls->m_unclabeling->setChecked(false);
        m_Controls->m_predlabelingBrush->setChecked(false);
        m_Controls->m_predlabeling->setChecked(false);

        m_UncertaintyLabel->SetFiberColors(255, 255, 255);
        this->CreateStreamlineInteractorBrush();
        m_StreamlineInteractorBrush->EnableInteraction(true);
        m_StreamlineInteractorBrush->LabelfromPrediction(false);
        m_StreamlineInteractorBrush->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractorBrush->SetToLabelNode(m_UncertaintyLabelNode);
    }
    else
    {
      m_StreamlineInteractorBrush->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkInteractiveFiberDissectionView::RemovefromDistance( bool checked )
{
    if (checked)
    {

        m_DistanceLabel->SetFiberColors(255, 255, 255);
        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->LabelfromPrediction(false);
        m_StreamlineInteractor->SetToLabelNode(m_DistanceLabelNode);
    }
    else
    {
      m_StreamlineInteractor->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkInteractiveFiberDissectionView::RemovefromPrediction( bool checked )
{
    if (checked)
    {
        m_Controls->m_predlabelingBrush->setChecked(false);
        m_Controls->m_unclabelingBrush->setChecked(false);
        m_Controls->m_unclabeling->setChecked(false);

//        m_Prediction->SetFiberColors(255, 255, 255);
        this->CreateStreamlineInteractor();
        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->LabelfromPrediction(true);
        m_StreamlineInteractor->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractor->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractor->SetToLabelNode(m_PredictionNode);
    }
    else
    {
      m_StreamlineInteractor->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
//      this->CleanTestArray();
    }

}

void QmitkInteractiveFiberDissectionView::RemovefromPredictionBrush( bool checked )
{
    if (checked)
    {
        m_Controls->m_predlabeling->setChecked(false);
        m_Controls->m_unclabelingBrush->setChecked(false);
        m_Controls->m_unclabeling->setChecked(false);

//        m_Prediction->SetFiberColors(255, 255, 255);
        this->CreateStreamlineInteractorBrush();
        m_StreamlineInteractorBrush->EnableInteraction(true);
        m_StreamlineInteractorBrush->LabelfromPrediction(true);
        m_StreamlineInteractorBrush->SetNegativeNode(m_negativeBundleNode);
        m_StreamlineInteractorBrush->SetPositiveNode(m_positiveBundleNode);
        m_StreamlineInteractorBrush->SetToLabelNode(m_PredictionNode);
    }
    else
    {
      m_StreamlineInteractorBrush->EnableInteraction(false);
//      m_StreamlineInteractor = nullptr;
//      this->CleanTestArray();
    }

}

void QmitkInteractiveFiberDissectionView::RemovefromSelection( bool checked )
{
    if (checked)
    {

    //        m_Prediction->SetFiberColors(255, 255, 255);
        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->LabelfromPrediction(true);
//        m_StreamlineInteractor->SetToLabelNode(m_SelectedFB.at(0));
        m_StreamlineInteractor->SetToLabelNode(m_testnode);
    }
    else
    {
      m_StreamlineInteractor->EnableInteraction(false);
    //      m_StreamlineInteractor = nullptr;
    }
}

void QmitkInteractiveFiberDissectionView::StartValidation()
{




    validater.reset();
//    mitk::DataNode::Pointer prednode = m_Controls->m_BundleBox->GetSelectedNode();
//    mitk::FiberBundle::Pointer pred = dynamic_cast<mitk::FiberBundle *>(prednode->GetData());
//    mitk::DataNode::Pointer gtnode = m_Controls->m_BundleBox->GetSelectedNode();
//    mitk::FiberBundle::Pointer gt = dynamic_cast<mitk::FiberBundle *>(gtnode->GetData());
//    mitk::FiberBundle::Pointer pred = dynamic_cast<mitk::FiberBundle*>(m_Controls->m_PredictionBox->GetSelectedNode()->GetData());
//    mitk::FiberBundle::Pointer gt = dynamic_cast<mitk::FiberBundle*>(m_Controls->m_GroundtruthBox->GetSelectedNode()->GetData());

    validater= std::make_shared<mitk::StreamlineFeatureExtractor>();
    validater->SetTractogramPrototypes(dynamic_cast<mitk::FiberBundle*>(m_Prototypes->GetData()), m_Controls->m_useStandardP->isChecked());
    MITK_INFO << "Prototypes loaded";
    validater->SetTractogramPrediction(dynamic_cast<mitk::FiberBundle*>(m_Controls->m_PredictionBox->GetSelectedNode()->GetData()));
    MITK_INFO << "Prediction loaded";
    validater->SetTractogramGroundtruth(dynamic_cast<mitk::FiberBundle*>(m_Controls->m_GroundtruthBox->GetSelectedNode()->GetData()));
    MITK_INFO << "Groundtruth loaded";

    validater->SetActiveCycle(m_activeCycleCounter);
//    validater->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData()), m_SelectedFB.at(0)->GetName());
    validater->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData()), m_testnode->GetName());
//    classifier->SetTractogramTest(dynamic_cast<mitk::FiberBundle*>(m_trainbundle->GetData()), m_trainbundle->GetName());



    MITK_INFO << "Testdata loaded";
    vnl_vector<float> metrics;
    metrics = validater->ValidationPipe();

   m_metrics.push_back(metrics);


    MITK_INFO << "Validation run succesfully";

    UpdateGui();

}

void QmitkInteractiveFiberDissectionView::AutomaticLabelling()
{

    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());
    
    
    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

     /* Check wether all Streamlines of the bundles are labeled... If all are labeled Skip for Loop*/
    unsigned int counter = 0;

    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();

    for (int i=0; i<fib->GetFiberPolyData()->GetNumberOfCells(); i++)
    {
      vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; j++)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      weights->InsertValue(counter, classifier->entropy_vector.at(i));
      vNewLines->InsertNextCell(container);
      counter++;

    }



    vNewPolyData->SetLines(vNewLines);
    vNewPolyData->SetPoints(vNewPoints);

    mitk::FiberBundle::Pointer uncertaintymap = mitk::FiberBundle::New(vNewPolyData);
    uncertaintymap->SetFiberWeights(weights);
    uncertaintymap->ColorFibersByFiberWeights(false, mitk::LookupTable::JET_TRANSPARENT);

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(uncertaintymap);
    node->SetName("uncertaintymap");
    this->GetDataStorage()->Add(node);
    MITK_INFO << 3;


}

void QmitkInteractiveFiberDissectionView::CleanTestArray()
{

    mitk::FiberBundle::Pointer fib_pos = dynamic_cast<mitk::FiberBundle*>(m_positiveBundleNode->GetData());
    mitk::FiberBundle::Pointer fib_neg = dynamic_cast<mitk::FiberBundle*>(m_negativeBundleNode->GetData());
//    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_testnode->GetData());
    vtkCell* cur_cell;
    vtkCell* cur_cell2;
    MITK_INFO << fib_neg->GetFiberPolyData()->GetNumberOfCells();
    MITK_INFO << fib_pos->GetFiberPolyData()->GetNumberOfCells();
    MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();


    std::vector<float> array;
    std::vector<float> array2;


    //    cv::parallel_for_(cv::Range(0, fib->GetFiberPolyData()->GetNumberOfCells()), [&](const cv::Range &range)
    //    {
    //    for (int i = range.start; i < range.end; i++)
        for (int i=0; i<fib->GetFiberPolyData()->GetNumberOfCells(); i++)
            {
                cur_cell = fib->GetFiberPolyData()->GetCell(i);
                auto numPoints = cur_cell->GetNumberOfPoints();
                vtkPoints* points = cur_cell->GetPoints();
                for (unsigned int j=0; j<numPoints; j++)
                {
                  double p[3];
                  points->GetPoint(j, p);
                  array.push_back(*p);
                }
                for (int ik=0; ik<fib_neg->GetFiberPolyData()->GetNumberOfCells(); ik++)
                {
                    cur_cell2 = fib_neg->GetFiberPolyData()->GetCell(ik);
                    auto numPoints2 = cur_cell2->GetNumberOfPoints();
                    vtkPoints* points2 = cur_cell2->GetPoints();
                    for (unsigned int j2=0; j2<numPoints2; j2++)
                    {
                      double p2[3];
                      points2->GetPoint(j2, p2);
                      array2.push_back(*p2);

                    }
                    if (array==array2)
                    {
                        fib->GetFiberPolyData()->DeleteCell(i);
                        MITK_INFO << "Delete Cell";
                    }
                    array2.clear();
                }
                for (int ik=0; ik<fib_pos->GetFiberPolyData()->GetNumberOfCells(); ik++)
                {

                    cur_cell2 = fib_pos->GetFiberPolyData()->GetCell(ik);
                    auto numPoints2 = cur_cell2->GetNumberOfPoints();
                    vtkPoints* points2 = cur_cell2->GetPoints();
                    for (unsigned int j2=0; j2<numPoints2; j2++)
                    {
                      double p2[3];
                      points2->GetPoint(j2, p2);
                      array2.push_back(*p2);
                    }
                    if (array==array2)
                    {
                        fib->GetFiberPolyData()->DeleteCell(i);
                        MITK_INFO << "Delete Cell";
                    }
                    array2.clear();
                }
                array.clear();
            }
    //    });

//        fib->GetFiberPolyData()->RemoveDeletedCells();

        MITK_INFO << fib->GetFiberPolyData()->GetNumberOfCells();



}

void QmitkInteractiveFiberDissectionView::ResetClassifier()
{
    m_IterationCounter = 0;
    m_RandomExtractionCounter = 0;
    m_activeCycleCounter = 0;
    m_createdStreamlineCounter = 0;

    classifier.reset();
    this->GetDataStorage()->Remove(m_positiveBundleNode);
    this->GetDataStorage()->Remove(m_negativeBundleNode);
    this->GetDataStorage()->Remove(m_negativeBundleNode);

    m_positiveBundleNode = NULL;
    m_negativeBundleNode = NULL;
    m_negativeBundleNode = NULL;
    m_positiveBundle = NULL;
    m_negativeBundle = NULL;
    m_negativeBundle = NULL;
    m_Controls->m_TrainClassifier->setEnabled(false);
    m_Controls->m_CreatePrediction->setEnabled(false);
    m_Controls->m_CreateUncertantyMap->setEnabled(false);



}
