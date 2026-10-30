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
#include <mitkDiffusionPropertyHelper.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkGibbsTrackingView.h"

#include <QmitkDiffusionViewControlsGuard.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QDir>
#include <QTimer>

// MITK
#include <mitkDiffusionMasks.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkImageAccessByItk.h>
#include <mitkIOUtil.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateFunction.h>

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
  // run() is called in the tracking thread, where an escaping exception would terminate the
  // application; the message is shown by AfterThread in the GUI thread.
  m_View->m_TrackingError.clear();
  try
  {
    m_View->m_GlobalTracker->Update();
  }
  catch (const itk::ExceptionObject &e)
  {
    m_View->m_TrackingError = e.GetDescription();
  }
  catch (const std::exception &e)
  {
    m_View->m_TrackingError = e.what();
  }
  catch (...)
  {
    m_View->m_TrackingError = "Unknown error.";
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
  , m_GlobalTracker(nullptr)
  , m_TrackingWorker(this)
{
  m_TrackingWorker.moveToThread(&m_TrackingThread);
  connect(&m_TrackingThread, SIGNAL(started()), this, SLOT(BeforeThread()));
  connect(&m_TrackingThread, SIGNAL(started()), &m_TrackingWorker, SLOT(run()));
  connect(&m_TrackingThread, SIGNAL(finished()), this, SLOT(AfterThread()));
  m_MonitorTimer.setInterval(1000);
  connect(&m_MonitorTimer, &QTimer::timeout, this, &QmitkGibbsTrackingView::TimerUpdate);
}

QmitkGibbsTrackingView::~QmitkGibbsTrackingView()
{
  // m_GlobalTracker is set before the thread starts and reset only after it finished
  if (!m_TrackingThread.isRunning())
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
  m_MonitorTimer.stop();

  UpdateGUI();

  // After an exception there is no result to show.
  if (!m_TrackingError.empty())
  {
    MITK_ERROR << "Gibbs tracking failed: " << m_TrackingError;
    QMessageBox::critical(nullptr, "Gibbs Tracking", QString("Tracking aborted:\n") + QString::fromStdString(m_TrackingError));
    m_FiberBundleNode = nullptr;
    m_GlobalTracker = nullptr;
    m_ItkOdfImage = nullptr;
    m_ItkTensorImage = nullptr;
    return;
  }

  if( !m_GlobalTracker->GetIsInValidState() )
  {
    QMessageBox::critical( nullptr, "Gibbs Tracking", "An internal error occured. Tracking aborted.\n Please check the log for details." );
    m_FiberBundleNode = nullptr;
    m_GlobalTracker = nullptr;
    m_ItkOdfImage = nullptr;
    m_ItkTensorImage = nullptr;
    return;
  }
  UpdateTrackingStatus();

  // Automatic values stay automatic (slider 0) for the next run, possibly on other data; the
  // estimate is only shown.
  if(m_Controls->m_ParticleWeightSlider->value()==0)
    m_Controls->m_ParticleWeightLabel->setText("auto (" + QString::number(m_GlobalTracker->GetParticleWeight()) + ")");
  if(m_Controls->m_ParticleWidthSlider->value()==0)
    m_Controls->m_ParticleWidthLabel->setText("auto (" + QString::number(m_GlobalTracker->GetParticleWidth()) + " mm)");
  if(m_Controls->m_ParticleLengthSlider->value()==0)
    m_Controls->m_ParticleLengthLabel->setText("auto (" + QString::number(m_GlobalTracker->GetParticleLength()) + " mm)");

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
  m_TrackingTimer.restart();
  m_MonitorTimer.start();

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
    QmitkDiffusionViewControlsGuard::Install(parent, this);
    m_Controls->setupUi( parent );

    AdvancedSettings();

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
    // any 3D scalar image: the mask values weight the particle density (e.g. a white-matter
    // probability map), voxels above 0.5 are seeded; a segmentation counts as 1 inside its labels
    auto isScalar3dImage = mitk::NodePredicateFunction::New([](const mitk::DataNode *node) {
      auto image = dynamic_cast<const mitk::Image*>(node->GetData());
      return image != nullptr && image->GetDimension() == 3 && image->GetPixelType().GetNumberOfComponents() == 1;
    });
    m_Controls->m_MaskImageBox->SetPredicate(mitk::NodePredicateOr::New(isScalar3dImage, mitk::GetSegmentationMaskPredicate()));

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

  const qint64 elapsedTime = m_TrackingTimer.elapsed()/1000;
  const qint64 hours = elapsedTime/3600;
  const qint64 minutes = (elapsedTime%3600)/60;
  const qint64 seconds = elapsedTime%60;

  m_Controls->m_ProposalAcceptance->setText(QString::number(m_GlobalTracker->GetProposalAcceptance()*100)+"%");

  m_Controls->m_TrackingTimeLabel->setText( QString::number(hours)+QString("h ")+QString::number(minutes)+QString("m ")+QString::number(seconds)+QString("s") );
  m_Controls->m_NumConnectionsLabel->setText( QString::number(m_GlobalTracker->GetNumConnections()) );
  m_Controls->m_NumParticlesLabel->setText( QString::number(m_GlobalTracker->GetNumParticles()) );
  if (m_GlobalTracker->GetIterations() > 0)
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

  // Both are reset: the input of an earlier run must not take precedence over the current one.
  m_OdfImage = nullptr;
  m_TensorImage = nullptr;
  if (dynamic_cast<mitk::OdfImage*>(m_ImageNode->GetData()))
    m_OdfImage = dynamic_cast<mitk::OdfImage*>(m_ImageNode->GetData());
  else if (dynamic_cast<mitk::TensorImage*>(m_ImageNode->GetData()))
    m_TensorImage = dynamic_cast<mitk::TensorImage*>(m_ImageNode->GetData());

  if (m_OdfImage.IsNull() && m_TensorImage.IsNull())
    return;

  // the line edit has no validator; invalid text would give an empty run
  bool iterationsOk = false;
  const double iterations = m_Controls->m_IterationsBox->text().toDouble(&iterationsOk);
  if (!iterationsOk || iterations < 1)
  {
    QMessageBox::information(nullptr, "Gibbs Tracking", "The number of iterations must be a positive number.");
    return;
  }

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

  // a mask that cannot be used must not silently be replaced by the whole image
  try{
    if(m_Controls->m_MaskImageBox->GetSelectedNode().IsNotNull())
    {
      m_MaskImageNode = m_Controls->m_MaskImageBox->GetSelectedNode();
      if (auto mask = mitk::GetMaskImage(m_MaskImageNode->GetData()); mask.IsNotNull())
        mitk::CastToItkImage(mask, m_MaskImage);
    }
  }
  catch (const std::exception &e)
  {
    QMessageBox::warning(nullptr, "Gibbs Tracking", QString("The mask image cannot be used:\n") + e.what());
    m_ItkOdfImage = nullptr;
    m_ItkTensorImage = nullptr;
    m_MaskImage = nullptr;
    m_TrackingNode = nullptr;
    return;
  }

  // The filter is configured here, in the GUI thread: the worker thread must not read widgets, and
  // the view must not see a half-assigned m_GlobalTracker.
  m_GlobalTracker = GibbsTrackingFilterType::New();
  m_GlobalTracker->SetOdfImage(m_ItkOdfImage);
  m_GlobalTracker->SetTensorImage(m_ItkTensorImage);
  m_GlobalTracker->SetMaskImage(m_MaskImage);
  m_GlobalTracker->SetStartTemperature((float)m_Controls->m_StartTempSlider->value()/100);
  m_GlobalTracker->SetEndTemperature((float)m_Controls->m_EndTempSlider->value()/10000);
  m_GlobalTracker->SetIterations(iterations);
  m_GlobalTracker->SetParticleWeight((float)m_Controls->m_ParticleWeightSlider->value()/10000);
  m_GlobalTracker->SetParticleWidth((float)(m_Controls->m_ParticleWidthSlider->value())/10);
  m_GlobalTracker->SetParticleLength((float)(m_Controls->m_ParticleLengthSlider->value())/10);
  m_GlobalTracker->SetInexBalance((float)m_Controls->m_InExBalanceSlider->value()/10);
  m_GlobalTracker->SetMinFiberLength(m_Controls->m_FiberLengthSlider->value());
  m_GlobalTracker->SetCurvatureThreshold(cos((float)m_Controls->m_CurvatureThresholdSlider->value()*itk::Math::pi/180));
  m_GlobalTracker->SetRandomSeed(m_Controls->m_RandomSeedSlider->value());

  // set before start(): BeforeThread is queued, and a second click must not restart the run
  m_ThreadIsRunning = true;
  UpdateGUI();

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

  // nullptr while the filter still prepares the data (no fibers built yet)
  vtkSmartPointer<vtkPolyData> fiberBundle = m_GlobalTracker->GetFiberBundle();
  if (fiberBundle == nullptr || fiberBundle->GetNumberOfLines() == 0)
  {
    if (!m_ThreadIsRunning)
      QMessageBox::information(nullptr, "Gibbs Tracking", "No fiber was accepted, so no tractogram was added or saved. More iterations, a lower particle weight or a shorter minimum fiber length may help.");
    return;
  }
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
  if (documentXML.SaveFile( filename.toStdString().c_str() ) != tinyxml2::XML_SUCCESS)
    QMessageBox::warning(nullptr, "Gibbs Tracking", "Could not save the parameters to\n" + filename);
}

// load current tracking paramters from xml file (.gtp)
void QmitkGibbsTrackingView::LoadTrackingParameters()
{
  QString filename = QFileDialog::getOpenFileName(0, tr("Load Parameters"), QDir::currentPath(), tr("Global Tracking Parameters (*.gtp)") );
  if(filename.isEmpty() || filename.isNull())
    return;

  tinyxml2::XMLDocument doc;
  tinyxml2::XMLElement* pElem = nullptr;
  // RootElement(), not FirstChild(): with an <?xml?> declaration, that is the first child.
  if (doc.LoadFile(filename.toStdString().c_str()) == tinyxml2::XML_SUCCESS && doc.RootElement() != nullptr)
    pElem = doc.RootElement()->FirstChildElement("parameter_set");
  if (pElem == nullptr)
  {
    QMessageBox::warning(nullptr, "Gibbs Tracking", "Not a Gibbs tracking parameter file (no <parameter_set>):\n" + filename);
    return;
  }
  // checked before any widget changes, so that an incomplete file changes nothing
  for (const char* name : {"iterations", "particle_length", "particle_width", "particle_weight", "temp_start",
                           "temp_end", "inexbalance", "fiber_length", "curvature_threshold"})
  {
    if (pElem->Attribute(name) == nullptr)
    {
      QMessageBox::warning(nullptr, "Gibbs Tracking", QString("The parameter file has no attribute %1:\n%2").arg(name).arg(filename));
      return;
    }
  }

  m_Controls->m_IterationsBox->setText(pElem->Attribute("iterations"));

  // Rounded to the sliders' steps (truncation turned 0.29 into 0.28); the labels are then set from
  // the sliders, so they show the values that are used, also when a slider clamps.
  auto toSlider = [pElem](const char* name, float scale) { return qRound(QString(pElem->Attribute(name)).toFloat()*scale); };
  m_Controls->m_ParticleLengthSlider->setValue(toSlider("particle_length", 10));
  m_Controls->m_ParticleWidthSlider->setValue(toSlider("particle_width", 10));
  m_Controls->m_ParticleWeightSlider->setValue(toSlider("particle_weight", 10000));
  m_Controls->m_StartTempSlider->setValue(toSlider("temp_start", 100));
  m_Controls->m_EndTempSlider->setValue(toSlider("temp_end", 10000));
  m_Controls->m_InExBalanceSlider->setValue(toSlider("inexbalance", 10));
  m_Controls->m_FiberLengthSlider->setValue(toSlider("fiber_length", 1));
  m_Controls->m_CurvatureThresholdSlider->setValue(toSlider("curvature_threshold", 1));

  SetParticleLength(m_Controls->m_ParticleLengthSlider->value());
  SetParticleWidth(m_Controls->m_ParticleWidthSlider->value());
  SetParticleWeight(m_Controls->m_ParticleWeightSlider->value());
  SetStartTemp(m_Controls->m_StartTempSlider->value());
  SetEndTemp(m_Controls->m_EndTempSlider->value());
  SetInExBalance(m_Controls->m_InExBalanceSlider->value());
  SetFiberLength(m_Controls->m_FiberLengthSlider->value());
  SetCurvatureThreshold(m_Controls->m_CurvatureThresholdSlider->value());
}
