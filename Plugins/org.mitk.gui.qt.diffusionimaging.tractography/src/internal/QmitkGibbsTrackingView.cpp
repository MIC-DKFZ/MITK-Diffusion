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
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkGibbsTrackingView.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QDir>
#include <QTimer>

// MITK
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkImageAccessByItk.h>
#include <mitkProgressBar.h>
#include <mitkIOUtil.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateProperty.h>

// ITK
#include <itkGibbsTrackingFilter.h>
#include <itkResampleImageFilter.h>
#include <itksys/SystemTools.hxx>

// MISC
#include <tinyxml2.h>




QmitkTrackingWorker::QmitkTrackingWorker(QmitkGibbsTrackingView* view)
  : m_View(view)
{
}

void QmitkTrackingWorker::run()
{
  m_View->m_GlobalTracker = QmitkGibbsTrackingView::GibbsTrackingFilterType::New();

  m_View->m_GlobalTracker->SetOdfImage(m_View->m_ItkOdfImage);
  m_View->m_GlobalTracker->SetTensorImage(m_View->m_ItkTensorImage);
  m_View->m_GlobalTracker->SetMaskImage(m_View->m_MaskImage);
  m_View->m_GlobalTracker->SetStartTemperature((float)m_View->m_Controls->m_StartTempSlider->value()/100);
  m_View->m_GlobalTracker->SetEndTemperature((float)m_View->m_Controls->m_EndTempSlider->value()/10000);
  m_View->m_GlobalTracker->SetIterations(m_View->m_Controls->m_IterationsBox->text().toDouble());
  m_View->m_GlobalTracker->SetParticleWeight((float)m_View->m_Controls->m_ParticleWeightSlider->value()/10000);
  m_View->m_GlobalTracker->SetParticleWidth((float)(m_View->m_Controls->m_ParticleWidthSlider->value())/10);
  m_View->m_GlobalTracker->SetParticleLength((float)(m_View->m_Controls->m_ParticleLengthSlider->value())/10);
  m_View->m_GlobalTracker->SetInexBalance((float)m_View->m_Controls->m_InExBalanceSlider->value()/10);
  m_View->m_GlobalTracker->SetMinFiberLength(m_View->m_Controls->m_FiberLengthSlider->value());
  m_View->m_GlobalTracker->SetCurvatureThreshold(cos((float)m_View->m_Controls->m_CurvatureThresholdSlider->value()*itk::Math::pi/180));
  m_View->m_GlobalTracker->SetRandomSeed(m_View->m_Controls->m_RandomSeedSlider->value());
  try{
    m_View->m_GlobalTracker->Update();
  }
  catch( mitk::Exception& e )
  {
    MITK_ERROR << "Internal error occured: " <<  e.what() << "\nAborting";
  }
  m_View->m_TrackingThread.quit();
}

const std::string QmitkGibbsTrackingView::VIEW_ID = "org.mitk.views.gibbstracking";

QmitkGibbsTrackingView::QmitkGibbsTrackingView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_TrackingNode(nullptr)
  , m_FiberBundle(nullptr)
  , m_MaskImage(nullptr)
  , m_TensorImage(nullptr)
  , m_OdfImage(nullptr)
  , m_ItkOdfImage(nullptr)
  , m_ItkTensorImage(nullptr)
  , m_ImageNode(nullptr)
  , m_MaskImageNode(nullptr)
  , m_FiberBundleNode(nullptr)
  , m_ThreadIsRunning(false)
  , m_ElapsedTime(0)
  , m_GlobalTracker(nullptr)
  , m_TrackingWorker(this)
{
  m_TrackingWorker.moveToThread(&m_TrackingThread);
  connect(&m_TrackingThread, SIGNAL(started()), this, SLOT(BeforeThread()));
  connect(&m_TrackingThread, SIGNAL(started()), &m_TrackingWorker, SLOT(run()));
  connect(&m_TrackingThread, SIGNAL(finished()), this, SLOT(AfterThread()));
}

QmitkGibbsTrackingView::~QmitkGibbsTrackingView()
{
  if (m_GlobalTracker.IsNull())
    return;

  m_GlobalTracker->SetAbortTracking(true);
  m_TrackingThread.wait();
}

// update tracking status and generate fiber bundle
void QmitkGibbsTrackingView::TimerUpdate()
{
  UpdateTrackingStatus();
  GenerateFiberBundle();
}

// tell global tractography filter to stop after current step
void QmitkGibbsTrackingView::StopGibbsTracking()
{
  if (m_GlobalTracker.IsNull())
    return;

  m_GlobalTracker->SetAbortTracking(true);
  m_Controls->m_TrackingStop->setEnabled(false);
  m_Controls->m_TrackingStop->setText("Stopping Tractography ...");
  m_TrackingNode = nullptr;
}

// update gui elements and generate fiber bundle after tracking is finished
void QmitkGibbsTrackingView::AfterThread()
{
  m_ThreadIsRunning = false;

  UpdateGUI();

  if( !m_GlobalTracker->GetIsInValidState() )
  {
    QMessageBox::critical( nullptr, "Gibbs Tracking", "An internal error occured. Tracking aborted.\n Please check the log for details." );
    m_FiberBundleNode = nullptr;
    return;
  }
  UpdateTrackingStatus();

  if(m_Controls->m_ParticleWeightSlider->value()==0)
  {
    m_Controls->m_ParticleWeightLabel->setText(QString::number(m_GlobalTracker->GetParticleWeight()));
    m_Controls->m_ParticleWeightSlider->setValue(m_GlobalTracker->GetParticleWeight()*10000);
  }
  if(m_Controls->m_ParticleWidthSlider->value()==0)
  {
    m_Controls->m_ParticleWidthLabel->setText(QString::number(m_GlobalTracker->GetParticleWidth()));
    m_Controls->m_ParticleWidthSlider->setValue(m_GlobalTracker->GetParticleWidth()*10);
  }
  if(m_Controls->m_ParticleLengthSlider->value()==0)
  {
    m_Controls->m_ParticleLengthLabel->setText(QString::number(m_GlobalTracker->GetParticleLength()));
    m_Controls->m_ParticleLengthSlider->setValue(m_GlobalTracker->GetParticleLength()*10);
  }

  GenerateFiberBundle();
  m_FiberBundleNode = 0;
  m_GlobalTracker = 0;

  // images not needed anymore ( relevant only for computation )
  // we need to release them to remove the memory access block created through CastToItk<> calls
  this->m_ItkOdfImage = 0;
  this->m_ItkTensorImage = 0;
}

// start tracking timer and update gui elements before tracking is started
void QmitkGibbsTrackingView::BeforeThread()
{
  m_ThreadIsRunning = true;
  m_TrackingTime = QTime::currentTime();
  m_ElapsedTime = 0;
  m_TrackingTimer.restart();

  UpdateGUI();
}

// setup gui elements and signal/slot connections
void QmitkGibbsTrackingView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkGibbsTrackingViewControls;
    m_Controls->setupUi( parent );

    AdvancedSettings();

    // connect( m_TrackingTimer, SIGNAL(elapsed()), this, SLOT(TimerUpdate()) );
    connect( m_Controls->m_TrackingStop, SIGNAL(clicked()), this, SLOT(StopGibbsTracking()) );
    connect( m_Controls->m_TrackingStart, SIGNAL(clicked()), this, SLOT(StartGibbsTracking()) );
    connect( m_Controls->m_AdvancedSettingsCheckbox, SIGNAL(clicked()), this, SLOT(AdvancedSettings()) );
    connect( m_Controls->m_SaveTrackingParameters, SIGNAL(clicked()), this, SLOT(SaveTrackingParameters()) );
    connect( m_Controls->m_LoadTrackingParameters, SIGNAL(clicked()), this, SLOT(LoadTrackingParameters()) );
    connect( m_Controls->m_ParticleWidthSlider, SIGNAL(valueChanged(int)), this, SLOT(SetParticleWidth(int)) );
    connect( m_Controls->m_ParticleLengthSlider, SIGNAL(valueChanged(int)), this, SLOT(SetParticleLength(int)) );
    connect( m_Controls->m_InExBalanceSlider, SIGNAL(valueChanged(int)), this, SLOT(SetInExBalance(int)) );
    connect( m_Controls->m_FiberLengthSlider, SIGNAL(valueChanged(int)), this, SLOT(SetFiberLength(int)) );
    connect( m_Controls->m_ParticleWeightSlider, SIGNAL(valueChanged(int)), this, SLOT(SetParticleWeight(int)) );
    connect( m_Controls->m_StartTempSlider, SIGNAL(valueChanged(int)), this, SLOT(SetStartTemp(int)) );
    connect( m_Controls->m_EndTempSlider, SIGNAL(valueChanged(int)), this, SLOT(SetEndTemp(int)) );
    connect( m_Controls->m_CurvatureThresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(SetCurvatureThreshold(int)) );
    connect( m_Controls->m_RandomSeedSlider, SIGNAL(valueChanged(int)), this, SLOT(SetRandomSeed(int)) );
    connect( m_Controls->m_OutputFileButton, SIGNAL(clicked()), this, SLOT(SetOutputFile()) );

    m_Controls->m_InputImageBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::OdfImage>::Pointer isOdfImagePredicate = mitk::TNodePredicateDataType<mitk::OdfImage>::New();
    mitk::TNodePredicateDataType<mitk::TensorImage>::Pointer isTensorImagePredicate = mitk::TNodePredicateDataType<mitk::TensorImage>::New();
    m_Controls->m_InputImageBox->SetPredicate( mitk::NodePredicateOr::New(isOdfImagePredicate, isTensorImagePredicate) );

    m_Controls->m_MaskImageBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_MaskImageBox->SetZeroEntryText("--");
    mitk::TNodePredicateDataType<mitk::Image>::Pointer isImagePredicate = mitk::TNodePredicateDataType<mitk::Image>::New();
    mitk::NodePredicateProperty::Pointer isBinaryPredicate = mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
    mitk::NodePredicateDimension::Pointer is3D = mitk::NodePredicateDimension::New(3);
    m_Controls->m_MaskImageBox->SetPredicate( mitk::NodePredicateAnd::New(isBinaryPredicate, mitk::NodePredicateAnd::New(isImagePredicate, is3D)) );

    connect( (QObject*)(m_Controls->m_MaskImageBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGUI()));
    connect( (QObject*)(m_Controls->m_InputImageBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGUI()));
  }
}

void QmitkGibbsTrackingView::SetFocus()
{
  m_Controls->m_TrackingStart->setFocus();
}

void QmitkGibbsTrackingView::SetInExBalance(int value)
{
  m_Controls->m_InExBalanceLabel->setText(QString::number((float)value/10));
}

void QmitkGibbsTrackingView::SetFiberLength(int value)
{
  m_Controls->m_FiberLengthLabel->setText(QString::number(value)+"mm");
}

void QmitkGibbsTrackingView::SetRandomSeed(int value)
{
  if (value>=0)
    m_Controls->m_RandomSeedLabel->setText(QString::number(value));
  else
    m_Controls->m_RandomSeedLabel->setText("auto");
}

void QmitkGibbsTrackingView::SetParticleWeight(int value)
{
  if (value>0)
    m_Controls->m_ParticleWeightLabel->setText(QString::number((float)value/10000));
  else
    m_Controls->m_ParticleWeightLabel->setText("auto");
}

void QmitkGibbsTrackingView::SetStartTemp(int value)
{
  m_Controls->m_StartTempLabel->setText(QString::number((float)value/100));
}

void QmitkGibbsTrackingView::SetEndTemp(int value)
{
  m_Controls->m_EndTempLabel->setText(QString::number((float)value/10000));
}

void QmitkGibbsTrackingView::SetParticleWidth(int value)
{
  if (value>0)
    m_Controls->m_ParticleWidthLabel->setText(QString::number((float)value/10)+" mm");
  else
    m_Controls->m_ParticleWidthLabel->setText("auto");
}

void QmitkGibbsTrackingView::SetParticleLength(int value)
{
  if (value>0)
    m_Controls->m_ParticleLengthLabel->setText(QString::number((float)value/10)+" mm");
  else
    m_Controls->m_ParticleLengthLabel->setText("auto");
}

void QmitkGibbsTrackingView::SetCurvatureThreshold(int value)
{
  m_Controls->m_CurvatureThresholdLabel->setText(QString::number(value)+"°");
}

// called if datamanager selection changes
void QmitkGibbsTrackingView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& )
{
  UpdateGUI();
}


void QmitkGibbsTrackingView::NodeRemoved(const mitk::DataNode * node)
{
  if (m_ThreadIsRunning)
  {
    if (node==m_TrackingNode.GetPointer())
    {
      StopGibbsTracking();
    }
  }
}

// update gui elements displaying trackings status
void QmitkGibbsTrackingView::UpdateTrackingStatus()
{
  if (m_GlobalTracker.IsNull())
    return;

  m_ElapsedTime += m_TrackingTimer.elapsed()/1000;
  unsigned long hours = m_ElapsedTime/3600;
  unsigned long minutes = (m_ElapsedTime%3600)/60;
  unsigned long seconds = m_ElapsedTime%60;

  m_Controls->m_ProposalAcceptance->setText(QString::number(m_GlobalTracker->GetProposalAcceptance()*100)+"%");

  m_Controls->m_TrackingTimeLabel->setText( QString::number(hours)+QString("h ")+QString::number(minutes)+QString("m ")+QString::number(seconds)+QString("s") );
  m_Controls->m_NumConnectionsLabel->setText( QString::number(m_GlobalTracker->GetNumConnections()) );
  m_Controls->m_NumParticlesLabel->setText( QString::number(m_GlobalTracker->GetNumParticles()) );
  m_Controls->m_CurrentStepLabel->setText( QString::number(100*m_GlobalTracker->GetCurrentIteration()/m_GlobalTracker->GetIterations())+"%" );
  m_Controls->m_AcceptedFibersLabel->setText( QString::number(m_GlobalTracker->GetNumAcceptedFibers()) );
}

// update gui elements (enable/disable elements and set tooltips)
void QmitkGibbsTrackingView::UpdateGUI()
{
  if (!m_ThreadIsRunning && m_Controls->m_InputImageBox->GetSelectedNode().IsNotNull())
  {
    m_Controls->m_TrackingStop->setEnabled(false);
    m_Controls->m_TrackingStart->setEnabled(true);
    m_Controls->m_LoadTrackingParameters->setEnabled(true);
    m_Controls->m_IterationsBox->setEnabled(true);
    m_Controls->m_AdvancedFrame->setEnabled(true);
    m_Controls->m_TrackingStop->setText("Stop Tractography");
    m_Controls->m_TrackingStart->setToolTip("Start tractography. No further change of parameters possible.");
    m_Controls->m_TrackingStop->setToolTip("");
    m_Controls->m_MaskImageBox->setEnabled(true);
    m_Controls->m_InputImageBox->setEnabled(true);
  }
  else if (!m_ThreadIsRunning)
  {
    m_Controls->m_TrackingStop->setEnabled(false);
    m_Controls->m_TrackingStart->setEnabled(false);
    m_Controls->m_LoadTrackingParameters->setEnabled(true);
    m_Controls->m_IterationsBox->setEnabled(true);
    m_Controls->m_AdvancedFrame->setEnabled(true);
    m_Controls->m_TrackingStop->setText("Stop Tractography");
    m_Controls->m_TrackingStart->setToolTip("No ODF image selected.");
    m_Controls->m_TrackingStop->setToolTip("");
    m_Controls->m_MaskImageBox->setEnabled(true);
    m_Controls->m_InputImageBox->setEnabled(true);
  }
  else
  {
    m_Controls->m_TrackingStop->setEnabled(true);
    m_Controls->m_TrackingStart->setEnabled(false);
    m_Controls->m_LoadTrackingParameters->setEnabled(false);
    m_Controls->m_IterationsBox->setEnabled(false);
    m_Controls->m_AdvancedFrame->setEnabled(false);
    m_Controls->m_AdvancedFrame->setVisible(false);
    m_Controls->m_AdvancedSettingsCheckbox->setChecked(false);
    m_Controls->m_TrackingStart->setToolTip("Tracking in progress.");
    m_Controls->m_TrackingStop->setToolTip("Stop tracking and display results.");
    m_Controls->m_MaskImageBox->setEnabled(false);
    m_Controls->m_InputImageBox->setEnabled(false);
  }
}

// show/hide advanced settings frame
void QmitkGibbsTrackingView::AdvancedSettings()
{
  m_Controls->m_AdvancedFrame->setVisible(m_Controls->m_AdvancedSettingsCheckbox->isChecked());
}

// check for mask and odf and start tracking thread
void QmitkGibbsTrackingView::StartGibbsTracking()
{
  if(m_ThreadIsRunning)
  {
    MITK_WARN("QmitkGibbsTrackingView")<<"Thread already running!";
    return;
  }
  m_GlobalTracker = nullptr;


  if (m_Controls->m_InputImageBox->GetSelectedNode().IsNull())
  {
    QMessageBox::information( nullptr, "Warning", "Please load and select a Odf image before starting image processing.");
    return;
  }
  m_ImageNode = m_Controls->m_InputImageBox->GetSelectedNode();

  if (dynamic_cast<mitk::OdfImage*>(m_ImageNode->GetData()))
    m_OdfImage = dynamic_cast<mitk::OdfImage*>(m_ImageNode->GetData());
  else if (dynamic_cast<mitk::TensorImage*>(m_ImageNode->GetData()))
    m_TensorImage = dynamic_cast<mitk::TensorImage*>(m_ImageNode->GetData());

  if (m_OdfImage.IsNull() && m_TensorImage.IsNull())
    return;

  // cast odf to itk
  m_TrackingNode = m_ImageNode;
  m_ItkTensorImage = nullptr;
  m_ItkOdfImage = nullptr;
  m_MaskImage = nullptr;

  if (m_OdfImage.IsNotNull())
  {
    m_ItkOdfImage = ItkOdfImgType::New();
    mitk::CastToItkImage(m_OdfImage, m_ItkOdfImage);
  }
  else
  {
    m_ItkTensorImage = ItkTensorImage::New();
    mitk::CastToItkImage(m_TensorImage, m_ItkTensorImage);
  }

  // mask image found?
  // catch exceptions thrown by the itkAccess macros
  try{
    if(m_Controls->m_MaskImageBox->GetSelectedNode().IsNotNull())
    {
      m_MaskImageNode = m_Controls->m_MaskImageBox->GetSelectedNode();
      if (dynamic_cast<mitk::Image*>(m_MaskImageNode->GetData()))
        mitk::CastToItkImage(dynamic_cast<mitk::Image*>(m_MaskImageNode->GetData()), m_MaskImage);
    }
  }
  catch(...){};

  // start worker thread
  m_TrackingThread.start(QThread::LowestPriority);
}

// generate mitkFiberBundle from tracking filter output
void QmitkGibbsTrackingView::GenerateFiberBundle()
{
  if (m_GlobalTracker.IsNull() || (!(m_Controls->m_VisualizationCheckbox->isChecked() || m_Controls->m_VisualizeOnceButton->isChecked()) && m_ThreadIsRunning))
    return;

  if (m_Controls->m_VisualizeOnceButton->isChecked())
    m_Controls->m_VisualizeOnceButton->setChecked(false);

  vtkSmartPointer<vtkPolyData> fiberBundle = m_GlobalTracker->GetFiberBundle();
  if ( m_GlobalTracker->GetNumAcceptedFibers()==0 )
    return;
  m_FiberBundle = mitk::FiberBundle::New(fiberBundle);
  m_FiberBundle->SetTrackVisHeader(dynamic_cast<mitk::Image*>(m_ImageNode->GetData())->GetGeometry());
  mitk::DiffusionPropertyHelper::CopyDICOMProperties(m_ImageNode->GetData(), m_FiberBundle);

  if (m_FiberBundleNode.IsNotNull()){
    GetDataStorage()->Remove(m_FiberBundleNode);
    m_FiberBundleNode = 0;
  }
  m_GlobalTracker->SetDicomProperties(m_FiberBundle);

  m_FiberBundleNode = mitk::DataNode::New();
  m_FiberBundleNode->SetData(m_FiberBundle);

  QString name("FiberBundle_");
  name += m_ImageNode->GetName().c_str();
  name += "_Gibbs";
  m_FiberBundleNode->SetName(name.toStdString());
  m_FiberBundleNode->SetVisibility(true);


  if (!m_OutputFileName.isEmpty() && !m_ThreadIsRunning)
  {
    try
    {
      mitk::IOUtil::Save(m_FiberBundle.GetPointer(),m_OutputFileName.toStdString());
      QMessageBox::information(nullptr, "Fiber bundle saved to", m_OutputFileName);
    }
    catch (itk::ExceptionObject &ex)
    {
      QMessageBox::information(nullptr, "Fiber bundle could not be saved", QString("%1\n%2\n%3\n%4\n%5\n%6").arg(ex.GetNameOfClass()).arg(ex.GetFile()).arg(ex.GetLine()).arg(ex.GetLocation()).arg(ex.what()).arg(ex.GetDescription()));
    }
  }
  GetDataStorage()->Add(m_FiberBundleNode, m_ImageNode);
}

void QmitkGibbsTrackingView::SetOutputFile()
{
  // SELECT FOLDER DIALOG
  m_OutputFileName = QFileDialog::getSaveFileName(0,
                                                  tr("Set file name"),
                                                  QDir::currentPath()+"/FiberBundle.fib",
                                                  tr("Fiber Bundle (*.fib)") );
  if (m_OutputFileName.isEmpty())
    m_Controls->m_OutputFileLabel->setText("N/A");
  else
    m_Controls->m_OutputFileLabel->setText(m_OutputFileName);
}

// save current tracking paramters as xml file (.gtp)
void QmitkGibbsTrackingView::SaveTrackingParameters()
{
  tinyxml2::XMLDocument documentXML;
  //TiXmlDeclaration* declXML = new TiXmlDeclaration( "1.0", "", "" );
  //documentXML.LinkEndChild( declXML );

  auto mainXML = documentXML.NewElement("global_tracking_parameter_file");
  mainXML->SetAttribute("file_version",  0.1);
  documentXML.InsertFirstChild(mainXML);

  auto paramXML = documentXML.NewElement("parameter_set");
  paramXML->SetAttribute("iterations", m_Controls->m_IterationsBox->text().toStdString().c_str());
  paramXML->SetAttribute("particle_length", (float)m_Controls->m_ParticleLengthSlider->value()/10);
  paramXML->SetAttribute("particle_width", (float)m_Controls->m_ParticleWidthSlider->value()/10);
  paramXML->SetAttribute("particle_weight", (float)m_Controls->m_ParticleWeightSlider->value()/10000);
  paramXML->SetAttribute("temp_start", (float)m_Controls->m_StartTempSlider->value()/100);
  paramXML->SetAttribute("temp_end", (float)m_Controls->m_EndTempSlider->value()/10000);
  paramXML->SetAttribute("inexbalance", (float)m_Controls->m_InExBalanceSlider->value()/10);
  paramXML->SetAttribute("fiber_length", m_Controls->m_FiberLengthSlider->value());
  paramXML->SetAttribute("curvature_threshold", m_Controls->m_CurvatureThresholdSlider->value());
  mainXML->InsertEndChild(paramXML);
  QString filename = QFileDialog::getSaveFileName(
                       0,
                       tr("Save Parameters"),
                       QDir::currentPath()+"/param.gtp",
                       tr("Global Tracking Parameters (*.gtp)") );

  if(filename.isEmpty() || filename.isNull())
    return;
  if(!filename.endsWith(".gtp"))
    filename += ".gtp";
  documentXML.SaveFile( filename.toStdString().c_str() );
}

// load current tracking paramters from xml file (.gtp)
void QmitkGibbsTrackingView::LoadTrackingParameters()
{
  QString filename = QFileDialog::getOpenFileName(0, tr("Load Parameters"), QDir::currentPath(), tr("Global Tracking Parameters (*.gtp)") );
  if(filename.isEmpty() || filename.isNull())
    return;

  tinyxml2::XMLDocument doc;
  doc.LoadFile(filename.toStdString().c_str());

  tinyxml2::XMLNode * hRoot = doc.FirstChild();
  auto pElem = hRoot->FirstChildElement("parameter_set");

  QString iterations(pElem->Attribute("iterations"));
  m_Controls->m_IterationsBox->setText(iterations);

  QString particleLength(pElem->Attribute("particle_length"));
  float pLength = particleLength.toFloat();
  QString particleWidth(pElem->Attribute("particle_width"));
  float pWidth = particleWidth.toFloat();

  if (pLength==0)
    m_Controls->m_ParticleLengthLabel->setText("auto");
  else
    m_Controls->m_ParticleLengthLabel->setText(particleLength+" mm");
  if (pWidth==0)
    m_Controls->m_ParticleWidthLabel->setText("auto");
  else
    m_Controls->m_ParticleWidthLabel->setText(particleWidth+" mm");

  m_Controls->m_ParticleWidthSlider->setValue(pWidth*10);
  m_Controls->m_ParticleLengthSlider->setValue(pLength*10);

  QString partWeight(pElem->Attribute("particle_weight"));
  m_Controls->m_ParticleWeightSlider->setValue(partWeight.toFloat()*10000);
  m_Controls->m_ParticleWeightLabel->setText(partWeight);

  QString startTemp(pElem->Attribute("temp_start"));
  m_Controls->m_StartTempSlider->setValue(startTemp.toFloat()*100);
  m_Controls->m_StartTempLabel->setText(startTemp);

  QString endTemp(pElem->Attribute("temp_end"));
  m_Controls->m_EndTempSlider->setValue(endTemp.toFloat()*10000);
  m_Controls->m_EndTempLabel->setText(endTemp);

  QString inExBalance(pElem->Attribute("inexbalance"));
  m_Controls->m_InExBalanceSlider->setValue(inExBalance.toFloat()*10);
  m_Controls->m_InExBalanceLabel->setText(inExBalance);

  QString fiberLength(pElem->Attribute("fiber_length"));
  m_Controls->m_FiberLengthSlider->setValue(fiberLength.toInt());
  m_Controls->m_FiberLengthLabel->setText(fiberLength+"mm");

  QString curvThres(pElem->Attribute("curvature_threshold"));
  m_Controls->m_CurvatureThresholdSlider->setValue(curvThres.toInt());
  m_Controls->m_CurvatureThresholdLabel->setText(curvThres+"°");
}
