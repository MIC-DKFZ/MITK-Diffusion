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
#include <QmitkRenderWindow.h> //Pointset
#include <QmitkPointListWidget.h> //Pointset

#include <QMessageBox>

#include <mitkNodePredicateProperty.h>
#include <mitkImageCast.h>
#include <mitkPointSet.h>
#include <mitkImageAccessByItk.h>
#include <mitkDataNodeObject.h>
#include <mitkTensorImage.h>

#include <mitkInteractionConst.h>
#include "usModuleRegistry.h"
//#include <itkFiberCurvatureFilter.h>

#include "mitkNodePredicateDataType.h"
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>


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

const std::string QmitkInteractiveFiberDissectionView::VIEW_ID = "org.mitk.views.interactivefiberdissection";
const std::string id_DataManager = "org.mitk.views.datamanager";
using namespace mitk;

QmitkInteractiveFiberDissectionView::QmitkInteractiveFiberDissectionView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_IterationCounter(0)
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
    m_Controls->m_selectedPointSetWidget->SetDataStorage(GetDataStorage());//pointset
    m_Controls->m_selectedPointSetWidget->SetNodePredicate(mitk::NodePredicateAnd::New(//pointset
      mitk::TNodePredicateDataType<mitk::PointSet>::New(),//pointset
      mitk::NodePredicateNot::New(mitk::NodePredicateOr::New(//pointset
        mitk::NodePredicateProperty::New("helper object"),//pointset
        mitk::NodePredicateProperty::New("hidden object")))));//pointset

    m_Controls->m_selectedPointSetWidget->SetSelectionIsOptional(true);//pointset
    m_Controls->m_selectedPointSetWidget->SetAutoSelectNewNodes(true);//pointset
    m_Controls->m_selectedPointSetWidget->SetEmptyInfo(QString("Please select a point set"));//pointset
    m_Controls->m_selectedPointSetWidget->SetPopUpTitel(QString("Select point set"));//pointsett

    connect(m_Controls->m_ErazorButton, SIGNAL( clicked() ), this, SLOT( RemovefromBundle() ) ); //need

    connect(m_Controls->m_StreamlineCreation, SIGNAL( clicked() ), this, SLOT( CreateStreamline()));

    connect(m_Controls->m_AddRandomFibers, SIGNAL( clicked() ), this, SLOT( ExtractRandomFibersFromTractogram() ) ); //need

    connect(m_Controls->m_TrainClassifier, SIGNAL( clicked() ), this, SLOT( StartAlgorithm( )));




    connect(m_Controls->m_addPointSetPushButton, &QPushButton::clicked,//pointset
      this, &QmitkInteractiveFiberDissectionView::OnAddPointSetClicked);//pointset
    connect(m_Controls->m_selectedPointSetWidget, &QmitkSingleNodeSelectionWidget::CurrentSelectionChanged,//pointset
      this, &QmitkInteractiveFiberDissectionView::OnCurrentSelectionChanged);//pointset




    auto renderWindowPart = this->GetRenderWindowPart();//pointset
    if (nullptr != renderWindowPart)//pointset
      this->RenderWindowPartActivated(renderWindowPart);//pointset
    this->OnCurrentSelectionChanged(m_Controls->m_selectedPointSetWidget->GetSelectedNodes());//pointset

  }

  UpdateGui();
}


void QmitkInteractiveFiberDissectionView::SetFocus()
{
  m_Controls->toolBoxx->setFocus();
  //m_Controls->m_addPointSetPushButton->setFocus();//pointset
}

void QmitkInteractiveFiberDissectionView::UpdateGui()
{
  m_Controls->m_FibLabel->setText("<font color='red'>mandatory</font>");
  m_Controls->m_InputData->setTitle("Please Select Input Data");

   // disable alle frames

  m_Controls->m_ErazorButton->setEnabled(false);



  m_Controls->m_addPointSetPushButton->setEnabled(false);
  m_Controls->m_StreamlineCreation->setEnabled(false);
  m_Controls->m_TrainClassifier->setEnabled(false);
  m_Controls->m_CreatePrediction->setEnabled(false);
  m_Controls->m_CreateUncertantyMap->setEnabled(false);
  m_Controls->m_Numtolabel->setEnabled(false);
  m_Controls->m_addPointSetPushButton->setEnabled(false);
  m_Controls->m_AddRandomFibers->setEnabled(false);
  m_Controls->m_AddUncertainFibers->setEnabled(false);

  bool fibSelected = !m_SelectedFB.empty();
  bool multipleFibsSelected = (m_SelectedFB.size()>1);
  bool sthSelected = m_SelectedImageNode.IsNotNull();
  bool psSelected = m_SelectedPS.IsNotNull();
  bool nfibSelected = !m_newfibersSelectedBundles.empty();
  bool posSelected = !m_positivSelectedBundles.empty();
  bool negSelected = !m_negativeSelectedBundles.IsNotNull();

  // toggle visibility of elements according to selected method


  // are fiber bundles selected?
  if ( fibSelected )
  {
    m_Controls->m_FibLabel->setText(QString(m_SelectedFB.at(0)->GetName().c_str()));
    m_Controls->m_addPointSetPushButton->setEnabled(true);
    m_Controls->m_AddRandomFibers->setEnabled(true);



    // more than two bundles needed to join/subtract
    if (multipleFibsSelected)
    {
      m_Controls->m_FibLabel->setText("multiple bundles selected");
    }



  }


  // is image selected
  if (sthSelected)
  {
    m_Controls->m_addPointSetPushButton->setEnabled(true);
  }

  if (psSelected)
  {
      m_Controls->m_StreamlineCreation->setEnabled(true);
  }

  if (nfibSelected)
  {

      m_Controls->m_ErazorButton->setEnabled(true);
  }


  if (posSelected && negSelected)
  {
      m_Controls->m_TrainClassifier->setEnabled(true);
  }


}

void QmitkInteractiveFiberDissectionView::OnEndInteraction()
{

}

void QmitkInteractiveFiberDissectionView::OnAddPointSetClicked()//pointset
{
  // ask for the name of the point set
  bool ok = false;
  QString name = QInputDialog::getText(QApplication::activeWindow(),
    tr("Add point set..."), tr("Enter name for the new point set"), QLineEdit::Normal, tr("PointSet").arg(++m_IterationCounter), &ok);
//  QString name = "PointSet";
  if (!ok || name.isEmpty())
  {
    return;
  }
  mitk::PointSet::Pointer pointSet = mitk::PointSet::New();
  mitk::DataNode::Pointer pointSetNode = mitk::DataNode::New();
  pointSetNode->SetData(pointSet);
  pointSetNode->SetProperty("name", mitk::StringProperty::New(name.toStdString()));
  pointSetNode->SetProperty("opacity", mitk::FloatProperty::New(1));
  pointSetNode->SetColor(1.0, 1.0, 0.0);
  this->GetDataStorage()->Add(pointSetNode, m_SelectedImageNode);


  m_Controls->m_selectedPointSetWidget->SetCurrentSelectedNode(pointSetNode);
}

void QmitkInteractiveFiberDissectionView::OnCurrentSelectionChanged(QmitkSingleNodeSelectionWidget::NodeList /*nodes*/)//pointset
{

  m_Controls->m_poinSetListWidget->SetPointSetNode(m_Controls->m_selectedPointSetWidget->GetSelectedNode());
  m_SelectedPS = m_Controls->m_selectedPointSetWidget->GetSelectedNode();

  UpdateGui();

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

void QmitkInteractiveFiberDissectionView::RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart)//pointset
{
  if (nullptr != m_Controls)
  {
    m_Controls->m_poinSetListWidget->AddSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("axial")->GetSliceNavigationController());
    m_Controls->m_poinSetListWidget->AddSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("sagittal")->GetSliceNavigationController());
    m_Controls->m_poinSetListWidget->AddSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("coronal")->GetSliceNavigationController());
  }
}

void QmitkInteractiveFiberDissectionView::RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart)//pointset
{
  if (nullptr != m_Controls)
  {
    m_Controls->m_poinSetListWidget->RemoveSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("axial")->GetSliceNavigationController());
    m_Controls->m_poinSetListWidget->RemoveSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("sagittal")->GetSliceNavigationController());
    m_Controls->m_poinSetListWidget->RemoveSliceNavigationController(renderWindowPart->GetQmitkRenderWindow("coronal")->GetSliceNavigationController());
  }
}

void QmitkInteractiveFiberDissectionView::CreateStreamline()
{

    if (m_positivSelectedBundles.empty())
    {
        mitk::DataNode::Pointer node = mitk::DataNode::New();

        m_positiveFibersData = vtkSmartPointer<vtkPolyData>::New();
        m_positiveFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
        m_positiveFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
        m_positiveBundle = mitk::FiberBundle:: New(m_positiveFibersData);

        node->SetData( m_positiveBundle );
        m_positivSelectedBundles.push_back(node);
        MITK_INFO << m_positivSelectedBundles.size();
        this->GetDataStorage()->Add(m_positivSelectedBundles.at(0));
      MITK_INFO << "Create Bundle";
    }

    if (!m_positivSelectedBundles.empty())
    {
        this->GetDataStorage()->Remove(m_positivSelectedBundles.at(m_positivSelectedBundles.size()-1));
        MITK_INFO << "Adding fibers";
    }


    mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet *>(m_SelectedPS->GetData());

    vnl_matrix<float> streamline;
    streamline.set_size(3, pointSet->GetSize());
    streamline.fill(0.0);


    mitk::PointSet::PointsIterator begin = pointSet->Begin();
    mitk::PointSet::PointsIterator end = pointSet->End();
    unsigned int i;
    mitk::PointSet::PointsContainer::Iterator it;

    for (it = begin, i = 0; it != end; ++it, ++i)
    {
        PointSet::PointType pt = pointSet->GetPoint(it->Index());
        vnl_vector_fixed< float, 3 > candV;
        candV[0]=pt[0]; candV[1]=pt[1]; candV[2]=pt[2];
        streamline.set_column(i, candV);
    }




    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();



    unsigned int counter = 0;
    for (unsigned int i=0; i<m_positiveFibersData->GetNumberOfCells(); ++i)
    {
      vtkCell* cell = m_positiveFibersData->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; ++j)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }

      vNewLines->InsertNextCell(container);
      counter++;
    }

      // build Fiber
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();


      for (unsigned int j=0; j<streamline.cols(); j++)
      {
        double p[3];
        p[0] = streamline.get(0,j);
        p[1] = streamline.get(1,j);
        p[2] = streamline.get(2,j);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      vNewLines->InsertNextCell(container);

      vNewPolyData->SetPoints(vNewPoints);
      vNewPolyData->SetLines(vNewLines);

      m_positiveFibersData = vtkSmartPointer<vtkPolyData>::New();
      m_positiveFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
      m_positiveFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
      m_positiveFibersData->SetPoints(vNewPoints);
      m_positiveFibersData->SetLines(vNewLines);

      m_positiveBundle = mitk::FiberBundle::New(vNewPolyData);
      m_positiveBundle->SetTrackVisHeader(dynamic_cast<mitk::Image*>(m_SelectedImageNode->GetData())->GetGeometry());
      m_positiveBundle->SetFiberColors(0, 255, 0);




      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(m_positiveBundle);
      node->SetName("+Bundle");
      m_positivSelectedBundles.push_back(node);
      MITK_INFO << m_positivSelectedBundles.size();



      MITK_INFO << "The + Bundle has Streamlines:";
      auto m_NegStreamline= dynamic_cast<mitk::FiberBundle *>(m_positivSelectedBundles.at(m_positivSelectedBundles.size()-1)->GetData());
      MITK_INFO << m_NegStreamline->GetFiberPolyData()->GetNumberOfCells();

      this->GetDataStorage()->Add(m_positivSelectedBundles.at(m_positivSelectedBundles.size()-1));

      UpdateGui();

}

void QmitkInteractiveFiberDissectionView::ExtractRandomFibersFromTractogram()
{
     MITK_INFO << "Number of Fibers to extract from Tractogram: ";
     MITK_INFO << m_Controls->m_NumRandomFibers->value();
     if (m_newfibersSelectedBundles.empty())
     {
         mitk::DataNode::Pointer node = mitk::DataNode::New();

         m_newfibersFibersData = vtkSmartPointer<vtkPolyData>::New();
         m_newfibersFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
         m_newfibersFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
         m_newfibersBundle = mitk::FiberBundle:: New(m_newfibersFibersData);

         node->SetData( m_newfibersBundle );
         m_newfibersSelectedBundles.push_back(node);
         MITK_INFO << m_newfibersSelectedBundles.size();
    //     this->GetDataStorage()->Add(m_newfibersSelectedBundles.at(0));
    //     UpdateGui();
       MITK_INFO << "Create Bundle";
     }

      mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());

    //   mitk::DataNode::Pointer node = mitk::DataNode::New();
    //   node ->SetData(pi);
    //   this->GetDataStorage()->Add(cur_bundle);


      vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
      vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();


      unsigned int counter = 0;
      for ( int i=0; i<m_Controls->m_NumRandomFibers->value(); i++)
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
    //    weights->InsertValue(counter, fib->GetFiberWeight(i));
        vNewLines->InsertNextCell(container);
        counter++;
        MITK_INFO << counter;
        MITK_INFO << vNewLines;

      }


      vNewPolyData->SetLines(vNewLines);
      vNewPolyData->SetPoints(vNewPoints);

      m_newfibersFibersData = vtkSmartPointer<vtkPolyData>::New();
      m_newfibersFibersData->SetPoints(vtkSmartPointer<vtkPoints>::New());
      m_newfibersFibersData->SetLines(vtkSmartPointer<vtkCellArray>::New());
      m_newfibersFibersData->SetPoints(vNewPoints);
      m_newfibersFibersData->SetLines(vNewLines);

      m_newfibersBundle = mitk::FiberBundle::New(vNewPolyData);
      m_newfibersBundle->SetFiberColors(255, 255, 255);

      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(m_newfibersBundle);
      node->SetName("ToLabel");
      m_newfibersSelectedBundles.push_back(node);
      MITK_INFO << m_newfibersSelectedBundles.size();

      this->GetDataStorage()->Add(m_newfibersSelectedBundles.at(m_newfibersSelectedBundles.size()-1));

    UpdateGui();


}

void QmitkInteractiveFiberDissectionView::RemovefromBundle()
{
    if (m_StreamlineInteractor.IsNull())
    {
        this->CreateStreamlineInteractor();
        mitk::FiberBundle::Pointer m_negativeBundle = mitk::FiberBundle::New();
        mitk::DataNode::Pointer m_negativeSelectedBundles = mitk::DataNode::New();
        m_negativeSelectedBundles->SetName("-Bundle");
        m_negativeSelectedBundles->SetData(m_negativeBundle);
        this->GetDataStorage()->Add(m_negativeSelectedBundles);

//        if (m_positivSelectedBundles.IsNull())
//        {

//            mitk::FiberBundle::Pointer m_positiveBundle = mitk::FiberBundle::New();
//            mitk::DataNode::Pointer m_positiveSelectedBundles = mitk::DataNode::New();
//            m_positiveSelectedBundles->SetName("-Bundle");
//            m_positiveSelectedBundles->SetData(m_positiveBundle);
//            this->GetDataStorage()->Add(m_positiveSelectedBundles);)
//        }

        m_StreamlineInteractor->EnableInteraction(true);
        m_StreamlineInteractor->SetNegativeNode(m_negativeSelectedBundles);
        m_StreamlineInteractor->SetPositiveNode(m_positivSelectedBundles.at(m_positivSelectedBundles.size()-1));
        m_StreamlineInteractor->SetToLabelNode(m_newfibersSelectedBundles.at(m_newfibersSelectedBundles.size()-1));
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


void QmitkInteractiveFiberDissectionView::StartAlgorithm()
{

}
