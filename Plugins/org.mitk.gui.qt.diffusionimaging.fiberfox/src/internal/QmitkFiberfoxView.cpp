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
#include "QmitkFiberfoxView.h"
#include <QmitkDiffusionViewControlsGuard.h>

// MITK
#include <mitkDiffusionMasks.h>
#include <mitkImage.h>
#include <mitkImageToItk.h>
#include <mitkImageCast.h>
#include <mitkProperties.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkDataStorage.h>
#include <itkTractsToDWIImageFilter.h>
#include <mitkTensorImage.h>
#include <mitkILinkedRenderWindowPart.h>
#include <mitkImageToItk.h>
#include <mitkImageCast.h>
#include <mitkImageGenerator.h>
#include <mitkNodePredicateDataType.h>
#include <itkScalableAffineTransform.h>
#include <mitkLevelWindowProperty.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateIsDWI.h>
#include <boost/property_tree/ptree.hpp>
#define RAPIDXML_NO_EXCEPTIONS
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <QFileDialog>
#include <QMessageBox>
#include <QPushButton>
#include <mitkMemoryUtilities.h>
#include "usModuleRegistry.h"
#include <mitkNoncentralChiNoiseModel.h>
#include <itksys/SystemTools.hxx>
#include <mitkIOUtil.h>
#include <QScrollBar>
#include <itkInvertIntensityImageFilter.h>
#include <QDialogButtonBox>
#include <itkAdcImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <mitkITKImageImport.h>

#include "mitkNodePredicateDataType.h"
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <algorithm>
#include <cmath>
#include <QFontDatabase>
#include <QMouseEvent>
#include <QPainter>
#include <QStyleOptionSizeGrip>
#include <QTextCharFormat>
#include <QTextCursor>

namespace
{
  /**
   * Handle in the lower right corner of a text field, below its vertical scroll bar: dragging it makes
   * the field taller or lower, like the corner of a text area in a web browser. The view lies in a
   * scroll area, which takes a taller page; the width follows the view.
   */
  class HeightGrip : public QWidget
  {
  public:
    explicit HeightGrip(QAbstractScrollArea* field)
      : QWidget(field), m_Field(field)
    {
      this->setCursor(Qt::SizeVerCursor);
      this->setToolTip("Drag to change the height of the field.");
    }

  protected:
    void paintEvent(QPaintEvent*) override
    {
      QPainter painter(this);
      QStyleOptionSizeGrip option;
      option.initFrom(this);
      option.corner = Qt::BottomRightCorner;
      this->style()->drawControl(QStyle::CE_SizeGrip, &option, &painter, this);
    }

    void mousePressEvent(QMouseEvent* event) override
    {
      m_StartY = event->globalPosition().y();
      m_StartHeight = m_Field->height();
    }

    void mouseMoveEvent(QMouseEvent* event) override
    {
      if (!(event->buttons() & Qt::LeftButton))
        return;
      // at least three lines of text
      const int minimum = 3*m_Field->fontMetrics().lineSpacing() + 2*m_Field->frameWidth();
      m_Field->setFixedHeight(std::max(minimum, m_StartHeight + static_cast<int>(std::lround(event->globalPosition().y() - m_StartY))));
    }

  private:
    QAbstractScrollArea* m_Field;
    double m_StartY = 0;
    int m_StartHeight = 0;
  };

  /**
   * The model names of a signal model string as the view writes it, one per model in list order
   * ("StickZeppelinBallBall"); empty if the string is not of that form.
   */
  std::vector<std::string> ModelNames(std::string modelString)
  {
    static const std::vector<std::string> names = {"Stick", "Zeppelin", "Tensor", "Ball", "Astrosticks", "Dot", "Prototype"};
    std::vector<std::string> result;
    while (!modelString.empty())
    {
      auto name = std::find_if(names.begin(), names.end(), [&](const std::string& n) { return modelString.compare(0, n.size(), n)==0; });
      if (name==names.end())
        return {};
      result.push_back(*name);
      modelString.erase(0, name->size());
    }
    return result;
  }

  /** The parameters of a started simulation as saved next to its image (ListMovedVolumes). */
  mitk::FiberfoxParameters SimulatedParameters(const mitk::FiberfoxParameters& parameters)
  {
    mitk::FiberfoxParameters simulated = parameters;
    simulated.ListMovedVolumes();
    return simulated;
  }
}

QmitkFiberfoxWorker::QmitkFiberfoxWorker(QmitkFiberfoxView* view)
  : m_View(view)
{

}

void QmitkFiberfoxWorker::run()
{
  // Written only here, in the worker thread; AfterThread reads it in the GUI thread after
  // finished(), i.e. after this function has returned.
  m_View->m_SimulationError.clear();
  try
  {
    m_View->m_TractsToDwiFilter->Update();
  }
  catch (const itk::ExceptionObject& e)
  {
    m_View->m_SimulationError = e.GetDescription();
  }
  catch (const std::exception& e)
  {
    m_View->m_SimulationError = e.what();
  }
  catch (...)
  {
    m_View->m_SimulationError = "Unknown error.";
  }
  m_View->m_Thread.quit();
}

const std::string QmitkFiberfoxView::VIEW_ID = "org.mitk.views.fiberfoxview";

QmitkFiberfoxView::QmitkFiberfoxView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_SelectedImageNode( nullptr )
  , m_Worker(this)
  , m_ThreadIsRunning(false)
{
  m_Worker.moveToThread(&m_Thread);
  connect(&m_Thread, SIGNAL(started()), this, SLOT(BeforeThread()));
  connect(&m_Thread, SIGNAL(started()), &m_Worker, SLOT(run()));
  connect(&m_Thread, SIGNAL(finished()), this, SLOT(AfterThread()));
  //    connect(&m_Thread, SIGNAL(terminated()), this, SLOT(AfterThread()));
  m_SimulationTimer = new QTimer(this);
}

void QmitkFiberfoxView::KillThread()
{
  MITK_INFO << "Aborting DWI simulation.";
  m_TractsToDwiFilter->SetAbortGenerateData(true);
  m_Controls->m_AbortSimulationButton->setEnabled(false);
  m_Controls->m_AbortSimulationButton->setText("Aborting simulation ...");
}

void QmitkFiberfoxView::BeforeThread()
{
  m_SimulationTime = QTime::currentTime();
  m_SimulationTimer->start(100);
  m_Controls->m_AbortSimulationButton->setVisible(true);
  m_Controls->m_GenerateImageButton->setVisible(false);
  m_Controls->m_SimulationStatusText->setVisible(true);
}

void QmitkFiberfoxView::AfterThread()
{
  UpdateSimulationStatus();
  m_SimulationTimer->stop();
  m_Controls->m_AbortSimulationButton->setVisible(false);
  m_Controls->m_AbortSimulationButton->setEnabled(true);
  m_Controls->m_AbortSimulationButton->setText("Abort simulation");
  m_Controls->m_GenerateImageButton->setVisible(true);
  m_ThreadIsRunning = false;

  QString statusText;
  mitk::FiberfoxParameters parameters;
  mitk::Image::Pointer mitkImage = mitk::Image::New();

  statusText = QString(m_TractsToDwiFilter->GetStatusText().c_str());
  if (m_TractsToDwiFilter->GetAbortGenerateData())
  {
    MITK_INFO << "Simulation aborted.";
    return;
  }

  // the filter attaches its output only at the very end, so a failed run has no image to add
  if (!m_SimulationError.empty())
  {
    MITK_ERROR << "Simulation failed: " << m_SimulationError;
    m_Controls->m_SimulationStatusText->append(QString("Simulation failed: ") + m_SimulationError.c_str());
    QMessageBox::warning(nullptr, "Simulation failed", m_SimulationError.c_str());
    return;
  }

  parameters = m_TractsToDwiFilter->GetParameters();

  mitkImage = mitk::GrabItkImageMemory( m_TractsToDwiFilter->GetOutput() );
  if (parameters.m_SignalGen.GetNumWeightedVolumes() > 0)
  {
    mitk::DiffusionPropertyHelper::SetGradientContainer(mitkImage, parameters.m_SignalGen.GetItkGradientContainer());
    mitk::DiffusionPropertyHelper::SetReferenceBValue(mitkImage, parameters.m_SignalGen.GetBvalue());
    mitk::DiffusionPropertyHelper::InitializeImage( mitkImage );
  }
  parameters.m_Misc.m_ResultNode->SetData( mitkImage );

  GetDataStorage()->Add(parameters.m_Misc.m_ResultNode, parameters.m_Misc.m_ParentNode);

  if (m_Controls->m_VolumeFractionsBox->isChecked())
  {
    if (m_TractsToDwiFilter->GetTickImage().IsNotNull())
    {
      mitk::Image::Pointer mitkImage = mitk::Image::New();
      itk::TractsToDWIImageFilter< short >::Float2DImageType::Pointer itkImage = m_TractsToDwiFilter->GetTickImage();
      mitkImage = mitk::GrabItkImageMemory( itkImage.GetPointer() );
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( mitkImage );
      node->SetName("Tick Image");
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
    }

    if (m_TractsToDwiFilter->GetRfImage().IsNotNull())
    {
      mitk::Image::Pointer mitkImage = mitk::Image::New();
      itk::TractsToDWIImageFilter< short >::Float2DImageType::Pointer itkImage = m_TractsToDwiFilter->GetRfImage();
      mitkImage = mitk::GrabItkImageMemory( itkImage.GetPointer() );
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( mitkImage );
      node->SetName("RF Image");
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
    }

    if (m_TractsToDwiFilter->GetPhaseImage().IsNotNull())
    {
      mitk::Image::Pointer phaseImage = mitk::Image::New();
      itk::TractsToDWIImageFilter< short >::DoubleDwiType::Pointer itkPhase = m_TractsToDwiFilter->GetPhaseImage();
      phaseImage = mitk::GrabItkImageMemory( itkPhase.GetPointer() );
      mitk::DataNode::Pointer phaseNode = mitk::DataNode::New();
      phaseNode->SetData( phaseImage );
      phaseNode->SetName("Phase Image");
      GetDataStorage()->Add(phaseNode, parameters.m_Misc.m_ResultNode);
    }

    if (m_TractsToDwiFilter->GetKspaceImage().IsNotNull())
    {
      mitk::Image::Pointer image = mitk::Image::New();
      itk::TractsToDWIImageFilter< short >::DoubleDwiType::Pointer itkImage = m_TractsToDwiFilter->GetKspaceImage();
      image = mitk::GrabItkImageMemory( itkImage.GetPointer() );
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( image );
      node->SetName("k-Space");
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
    }

    {
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(m_TractsToDwiFilter->GetCoilPointset());
      node->SetName("Coil Positions");
      node->SetProperty("pointsize", mitk::FloatProperty::New(parameters.m_SignalGen.m_ImageSpacing[0]/4));
      node->SetProperty("color", mitk::ColorProperty::New(0, 1, 0));
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
    }

    int c = 1;
    std::vector< itk::TractsToDWIImageFilter< short >::DoubleDwiType::Pointer > output_real = m_TractsToDwiFilter->GetOutputImagesReal();
    for (auto real : output_real)
    {
      mitk::Image::Pointer image = mitk::Image::New();
      image->InitializeByItk(real.GetPointer());
      image->SetVolume(real->GetBufferPointer());

      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( image );
      node->SetName("Coil-"+QString::number(c).toStdString()+"-real");
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
      ++c;
    }

    c = 1;
    std::vector< itk::TractsToDWIImageFilter< short >::DoubleDwiType::Pointer > output_imag = m_TractsToDwiFilter->GetOutputImagesImag();
    for (auto imag : output_imag)
    {
      mitk::Image::Pointer image = mitk::Image::New();
      image->InitializeByItk(imag.GetPointer());
      image->SetVolume(imag->GetBufferPointer());

      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( image );
      node->SetName("Coil-"+QString::number(c).toStdString()+"-imag");
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
      ++c;
    }

    std::vector< itk::TractsToDWIImageFilter< short >::ItkDoubleImgType::Pointer > volumeFractions = m_TractsToDwiFilter->GetVolumeFractions();
    for (unsigned int k=0; k<volumeFractions.size(); k++)
    {
      mitk::Image::Pointer image = mitk::Image::New();
      image->InitializeByItk(volumeFractions.at(k).GetPointer());
      image->SetVolume(volumeFractions.at(k)->GetBufferPointer());

      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData( image );
      node->SetName("CompartmentVolume-"+QString::number(k).toStdString());
      GetDataStorage()->Add(node, parameters.m_Misc.m_ResultNode);
    }
  }
  m_TractsToDwiFilter = nullptr;

  if (parameters.m_Misc.m_AfterSimulationMessage.size()>0)
    QMessageBox::information( nullptr, "Warning", parameters.m_Misc.m_AfterSimulationMessage.c_str());

  mitk::BaseData::Pointer basedata = parameters.m_Misc.m_ResultNode->GetData();
  if (basedata.IsNotNull())
  {
    mitk::RenderingManager::GetInstance()->InitializeViews(
          basedata->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true );
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }

  if (!parameters.m_Misc.m_OutputPath.empty())
  {
    try{
      // dots only in the node name (b-value, noise variance), the folder may contain some
      QString name(parameters.m_Misc.m_ResultNode->GetName().c_str());
      name.replace(QString("."), QString("_"));
      QString outputFileName = QString(parameters.m_Misc.m_OutputPath.c_str()) + name;
      // the parameters of this run (FSL bvecs, sampled prototypes), not the GUI's current state
      m_SimulatedParameters.SaveParameters((outputFileName+".ffp").toStdString());
      outputFileName += ".dwi";
      QString status("Saving output image to ");
      status += outputFileName;
      m_Controls->m_SimulationStatusText->append(status);
      mitk::IOUtil::Save(mitkImage, outputFileName.toStdString());
      m_Controls->m_SimulationStatusText->append("File saved successfully.");

    }
    catch (itk::ExceptionObject &e)
    {
      QString status("Exception during DWI writing: ");
      status += e.GetDescription();
      m_Controls->m_SimulationStatusText->append(status);
    }
    catch (...)
    {
      m_Controls->m_SimulationStatusText->append("Unknown exception during DWI writing!");
    }
  }
}

void QmitkFiberfoxView::UpdateSimulationStatus()
{
  QString statusText = QString(m_TractsToDwiFilter->GetStatusText().c_str());

  // rewritten only when the log changed, so that the user can scroll back during a run
  if (QString::compare(m_SimulationStatusText,statusText)!=0)
  {
    m_SimulationStatusText = statusText;
    m_Controls->m_SimulationStatusText->clear();
    m_Controls->m_SimulationStatusText->setPlainText(statusText);
    // The percentages above a progress bar line up with its stars only in a fixed-pitch font. The
    // theme's font rule replaces the font of the widget, not that of the text's characters.
    QTextCharFormat fixedPitch;
    fixedPitch.setFontFamilies({QFontDatabase::systemFont(QFontDatabase::FixedFont).family()});
    fixedPitch.setFontFixedPitch(true);
    QTextCursor cursor(m_Controls->m_SimulationStatusText->document());
    cursor.select(QTextCursor::Document);
    cursor.mergeCharFormat(fixedPitch);
    QScrollBar *vScrollBar = m_Controls->m_SimulationStatusText->verticalScrollBar();
    vScrollBar->triggerAction(QScrollBar::SliderToMaximum);
  }
}

// Destructor
QmitkFiberfoxView::~QmitkFiberfoxView()
{
  // Destroying a running QThread aborts the Workbench, and the worker and the filter use this view,
  // so a running simulation is aborted and joined. AfterThread must not run any more: the widgets
  // are gone. The abort flag is set repeatedly because the filter resets it when it starts.
  if (m_Thread.isRunning())
  {
    disconnect(&m_Thread, SIGNAL(finished()), this, SLOT(AfterThread()));
    do
      m_TractsToDwiFilter->SetAbortGenerateData(true);
    while (!m_Thread.wait(100));
  }
  delete m_SimulationTimer;
  delete m_Controls;
}

void QmitkFiberfoxView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkFiberfoxViewControls;
    // the selection slots of its boxes update other widgets (UpdateGui)
    QmitkDiffusionViewControlsGuard::Install(parent, this);
    m_Controls->setupUi( parent );

    m_Controls->m_StickWidget1->setVisible(true);
    m_Controls->m_StickWidget2->setVisible(false);
    m_Controls->m_ZeppelinWidget1->setVisible(false);
    m_Controls->m_ZeppelinWidget2->setVisible(false);
    m_Controls->m_TensorWidget1->setVisible(false);
    m_Controls->m_TensorWidget2->setVisible(false);

    m_Controls->m_BallWidget1->setVisible(true);
    m_Controls->m_BallWidget2->setVisible(false);
    m_Controls->m_BallWidget2->SetT1(4658);
    m_Controls->m_BallWidget2->SetT2(2200);
    m_Controls->m_AstrosticksWidget1->setVisible(false);
    m_Controls->m_AstrosticksWidget2->setVisible(false);
    m_Controls->m_AstrosticksWidget2->SetT1(4658);
    m_Controls->m_AstrosticksWidget2->SetT2(2200);

    // the handle that changes the height of the simulation log is the corner widget, which takes the
    // place below a vertical scroll bar that is always shown
    m_Controls->m_SimulationStatusText->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    m_Controls->m_SimulationStatusText->setCornerWidget(new HeightGrip(m_Controls->m_SimulationStatusText));
    m_Controls->m_DotWidget1->setVisible(false);
    m_Controls->m_DotWidget2->setVisible(false);
    m_Controls->m_DotWidget2->SetT1(4658);
    m_Controls->m_DotWidget2->SetT2(2200);

    m_Controls->m_PrototypeWidget1->setVisible(false);
    m_Controls->m_PrototypeWidget3->setVisible(false);
    m_Controls->m_PrototypeWidget4->setVisible(false);

    m_Controls->m_PrototypeWidget3->SetMinFa(0.0);
    m_Controls->m_PrototypeWidget3->SetMaxFa(0.15);
    m_Controls->m_PrototypeWidget4->SetMinFa(0.0);
    m_Controls->m_PrototypeWidget4->SetMaxFa(0.15);
    m_Controls->m_PrototypeWidget3->SetMinAdc(0.0);
    m_Controls->m_PrototypeWidget3->SetMaxAdc(0.001);
    m_Controls->m_PrototypeWidget4->SetMinAdc(0.003);
    m_Controls->m_PrototypeWidget4->SetMaxAdc(0.004);

    m_Controls->m_Comp2FractionFrame->setVisible(false);
    m_Controls->m_Comp4FractionFrame->setVisible(false);
    m_Controls->m_DiffusionPropsMessage->setVisible(false);
    m_Controls->m_GeometryMessage->setVisible(false);
    m_Controls->m_AdvancedSignalOptionsFrame->setVisible(false);
    m_Controls->m_NoiseFrame->setVisible(false);
    m_Controls->m_NoiseCoilsLabel->setVisible(false);
    m_Controls->m_NoiseCoilsBox->setVisible(false);
    m_Controls->m_ZeroRinging->setVisible(false);
    m_Controls->m_GhostFrame->setVisible(false);
    m_Controls->m_DistortionsFrame->setVisible(false);
    m_Controls->m_EddyFrame->setVisible(false);
    m_Controls->m_SpikeFrame->setVisible(false);
    m_Controls->m_AliasingFrame->setVisible(false);
    m_Controls->m_MotionArtifactFrame->setVisible(false);
    m_Controls->m_DriftFrame->setVisible(false);
    ShowKspaceTimeSegments(m_Parameters.m_SignalGen.m_KspaceTimeSegments);   // the parameters' default
    m_ParameterFile = QDir::currentPath()+"/param.ffp";

    m_Controls->m_AbortSimulationButton->setVisible(false);
    m_Controls->m_SimulationStatusText->setVisible(false);

    m_Controls->m_FrequencyMapBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_Comp1VolumeFraction->SetDataStorage(this->GetDataStorage());
    m_Controls->m_Comp2VolumeFraction->SetDataStorage(this->GetDataStorage());
    m_Controls->m_Comp3VolumeFraction->SetDataStorage(this->GetDataStorage());
    m_Controls->m_Comp4VolumeFraction->SetDataStorage(this->GetDataStorage());
    m_Controls->m_MaskComboBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_TemplateComboBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_FiberBundleComboBox->SetDataStorage(this->GetDataStorage());

    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isFiberBundle = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    mitk::TNodePredicateDataType<mitk::Image>::Pointer isMitkImage = mitk::TNodePredicateDataType<mitk::Image>::New();
    mitk::NodePredicateIsDWI::Pointer isDwi = mitk::NodePredicateIsDWI::New( );
    mitk::NodePredicateDataType::Pointer isDti = mitk::NodePredicateDataType::New("TensorImage");
    mitk::NodePredicateDataType::Pointer isOdf = mitk::NodePredicateDataType::New("Odfmage");
    mitk::NodePredicateOr::Pointer isDiffusionImage = mitk::NodePredicateOr::New(isDwi, isDti);
    isDiffusionImage = mitk::NodePredicateOr::New(isDiffusionImage, isOdf);
    mitk::NodePredicateNot::Pointer noDiffusionImage = mitk::NodePredicateNot::New(isDiffusionImage);
    mitk::NodePredicateAnd::Pointer isNonDiffMitkImage = mitk::NodePredicateAnd::New(isMitkImage, noDiffusionImage);

    m_Controls->m_FrequencyMapBox->SetPredicate(isNonDiffMitkImage);
    m_Controls->m_Comp1VolumeFraction->SetPredicate(isNonDiffMitkImage);
    m_Controls->m_Comp1VolumeFraction->SetZeroEntryText("--");
    m_Controls->m_Comp2VolumeFraction->SetPredicate(isNonDiffMitkImage);
    m_Controls->m_Comp2VolumeFraction->SetZeroEntryText("--");
    m_Controls->m_Comp3VolumeFraction->SetPredicate(isNonDiffMitkImage);
    m_Controls->m_Comp3VolumeFraction->SetZeroEntryText("--");
    m_Controls->m_Comp4VolumeFraction->SetPredicate(isNonDiffMitkImage);
    m_Controls->m_Comp4VolumeFraction->SetZeroEntryText("--");
    m_Controls->m_MaskComboBox->SetPredicate(mitk::GetBinaryMaskPredicate());
    m_Controls->m_MaskComboBox->SetZeroEntryText("--");
    m_Controls->m_TemplateComboBox->SetPredicate(isMitkImage);
    m_Controls->m_TemplateComboBox->SetZeroEntryText("--");
    m_Controls->m_FiberBundleComboBox->SetPredicate(isFiberBundle);
    m_Controls->m_FiberBundleComboBox->SetZeroEntryText("--");

    QFont font;
    font.setFamily("Courier");
    font.setStyleHint(QFont::Monospace);
    font.setFixedPitch(true);
    font.setPointSize(7);
    m_Controls->m_SimulationStatusText->setFont(font);

    connect( m_SimulationTimer, SIGNAL(timeout()), this, SLOT(UpdateSimulationStatus()) );
    connect((QObject*) m_Controls->m_AbortSimulationButton, SIGNAL(clicked()), (QObject*) this, SLOT(KillThread()));
    connect((QObject*) m_Controls->m_GenerateImageButton, SIGNAL(clicked()), (QObject*) this, SLOT(GenerateImage()));
    connect((QObject*) m_Controls->m_AddNoise, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddNoise(int)));
    // the noise channel count is a parameter of the noncentral chi distribution (index 2) only
    connect(m_Controls->m_NoiseDistributionBox, &QComboBox::currentIndexChanged, this, [this](int index)
    {
      m_Controls->m_NoiseCoilsLabel->setVisible(index==2);
      m_Controls->m_NoiseCoilsBox->setVisible(index==2);
    });
    connect((QObject*) m_Controls->m_AddGhosts, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddGhosts(int)));
    connect((QObject*) m_Controls->m_AddDistortions, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddDistortions(int)));
    connect((QObject*) m_Controls->m_AddEddy, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddEddy(int)));
    // the segment count is a setting of the fixed time segments (index 2) only
    connect(m_Controls->m_KspacePhaseBox, &QComboBox::currentIndexChanged, this, [this](int index)
    {
      m_Controls->m_KspaceTimeSegmentsLabel->setVisible(index==2);
      m_Controls->m_KspaceTimeSegmentsBox->setVisible(index==2);
    });
    connect((QObject*) m_Controls->m_AddSpikes, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddSpikes(int)));
    connect((QObject*) m_Controls->m_AddAliasing, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddAliasing(int)));
    connect((QObject*) m_Controls->m_AddMotion, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddMotion(int)));
    connect((QObject*) m_Controls->m_AddDrift, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddDrift(int)));
    connect((QObject*) m_Controls->m_AddGibbsRinging, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddRinging(int)));

    connect((QObject*) m_Controls->m_Compartment1Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp1ModelFrameVisibility(int)));
    connect((QObject*) m_Controls->m_Compartment2Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp2ModelFrameVisibility(int)));
    connect((QObject*) m_Controls->m_Compartment3Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp3ModelFrameVisibility(int)));
    connect((QObject*) m_Controls->m_Compartment4Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp4ModelFrameVisibility(int)));

    connect((QObject*) m_Controls->m_AdvancedOptionsBox_2, SIGNAL( stateChanged(int)), (QObject*) this, SLOT(ShowAdvancedOptions(int)));
    connect((QObject*) m_Controls->m_UseBvalsBvecsBox, SIGNAL( stateChanged(int)), (QObject*) this, SLOT(OnBvalsBvecsCheck(int)));

    connect((QObject*) m_Controls->m_SaveParametersButton, SIGNAL(clicked()), (QObject*) this, SLOT(SaveParameters()));
    connect((QObject*) m_Controls->m_LoadParametersButton, SIGNAL(clicked()), (QObject*) this, SLOT(LoadParameters()));
    connect((QObject*) m_Controls->m_OutputPathButton, SIGNAL(clicked()), (QObject*) this, SLOT(SetOutputPath()));
    connect((QObject*) m_Controls->m_LoadBvalsButton, SIGNAL(clicked()), (QObject*) this, SLOT(SetBvalsEdit()));
    connect((QObject*) m_Controls->m_LoadBvecsButton, SIGNAL(clicked()), (QObject*) this, SLOT(SetBvecsEdit()));

    connect((QObject*) m_Controls->m_MaskComboBox, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(OnMaskSelected(int)));
    connect((QObject*) m_Controls->m_TemplateComboBox, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(OnTemplateSelected(int)));
    connect((QObject*) m_Controls->m_FiberBundleComboBox, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(OnFibSelected(int)));
    connect((QObject*) m_Controls->m_LineReadoutTimeBox, SIGNAL( valueChanged(double)), (QObject*) this, SLOT(OnTlineChanged()));
    connect((QObject*) m_Controls->m_SizeX, SIGNAL( valueChanged(int)), (QObject*) this, SLOT(OnTlineChanged()));
    connect((QObject*) m_Controls->m_NumGradientsBox, SIGNAL( valueChanged(int)), (QObject*) this, SLOT(OnGradientSettingsChanged()));
    connect((QObject*) m_Controls->m_BvalueBox, SIGNAL( valueChanged(int)), (QObject*) this, SLOT(OnGradientSettingsChanged()));
    connect((QObject*) m_Controls->m_LoadBvalsEdit, SIGNAL( textChanged(QString)), (QObject*) this, SLOT(OnGradientSettingsChanged()));
    connect((QObject*) m_Controls->m_LoadBvecsEdit, SIGNAL( textChanged(QString)), (QObject*) this, SLOT(OnGradientSettingsChanged()));
  }
  OnTlineChanged();
  UpdateGui();
}

void QmitkFiberfoxView::OnTlineChanged()
{
  double num_pix_line = 0;
  if (m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull())   // use geometry of selected image
  {
    mitk::Image::Pointer img = dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData());
    num_pix_line = img->GetDimension(0);
  }
  else if (m_Controls->m_MaskComboBox->GetSelectedNode().IsNotNull())   // use geometry of mask image
  {
    // the mask may be a segmentation, which is not an image: its index extent is the size
    num_pix_line = m_Controls->m_MaskComboBox->GetSelectedNode()->GetData()->GetGeometry()->GetExtent(0);
  }
  else
  {
     num_pix_line = m_Controls->m_SizeX->value();
  }

  double value = static_cast<double>(m_Controls->m_LineReadoutTimeBox->value())/1000.0;  // line readout time in seconds
  double dweel_pp = value/num_pix_line; // pixel readout time in seconds
  double bw = 1/dweel_pp;
  double bw_pp = bw/num_pix_line;
  std::string tt = "Bandwidth:\n" + boost::lexical_cast<std::string>(itk::Math::Round<int, double>(bw)) + "Hz\n" + boost::lexical_cast<std::string>(itk::Math::Round<int, double>(bw_pp)) + "Hz/Px";
  m_Controls->m_LineReadoutTimeBox->setToolTip(tt.c_str());
}

void QmitkFiberfoxView::OnMaskSelected(int )
{
  UpdateGui();
}

void QmitkFiberfoxView::OnTemplateSelected(int )
{
  UpdateGui();
}

void QmitkFiberfoxView::OnFibSelected(int )
{
  UpdateGui();
}

void QmitkFiberfoxView::OnBvalsBvecsCheck(int )
{
  OnGradientSettingsChanged();
  UpdateGui();
}

namespace
{
  std::vector<double> PrototypeSettings(QmitkPrototypeSignalParametersWidget* widget)
  {
    return {widget->GetNumberOfSamples(), widget->GetMinFa(), widget->GetMaxFa(), widget->GetMinAdc(), widget->GetMaxAdc()};
  }

  /** The non-zero b-values of the gradients, sorted. */
  std::vector<int> Shells(mitk::SignalGenerationParameters& signalGen)
  {
    std::vector<int> shells = signalGen.GetBvalues();
    shells.erase(std::remove(shells.begin(), shells.end(), 0), shells.end());
    std::sort(shells.begin(), shells.end());
    return shells;
  }
}

QmitkPrototypeSignalParametersWidget* QmitkFiberfoxView::GetPrototypeWidget(int compartmentId)
{
  if (compartmentId<1 || compartmentId>4)
    return nullptr;
  // nullptr for compartment 2, which offers no prototype model
  return static_cast<QmitkPrototypeSignalParametersWidget*>(ModelWidget(Compartments().at(compartmentId-1), CompartmentModel::Prototype));
}

mitk::RawShModel<>* QmitkFiberfoxView::NewPrototypeModel(QmitkPrototypeSignalParametersWidget* widget, int compartmentId)
{
  mitk::RawShModel<>* model = nullptr;
  auto loaded = m_LoadedPrototypes.find(compartmentId);
  // The kernels of a loaded .ffp, as long as the widget shows the file's settings and the simulated
  // b-values (from boxes, a template DWI or bvals/bvecs files) are the file's: kernels hold the
  // signal of their shell, they are evaluated for the gradient direction only.
  if (loaded!=m_LoadedPrototypes.end() && loaded->second.m_WidgetSettings==PrototypeSettings(widget)
      && loaded->second.m_Shells==Shells(m_Parameters.m_SignalGen))
  {
    model = static_cast< mitk::RawShModel<>* >(loaded->second.m_Model->Clone());
  }
  else
  {
    model = new mitk::RawShModel<>();
    model->SetMaxNumKernels(widget->GetNumberOfSamples());
    model->SetFaRange(widget->GetMinFa(), widget->GetMaxFa());
    model->SetAdcRange(widget->GetMinAdc(), widget->GetMaxAdc());
  }
  model->SetGradientList(m_Parameters.m_SignalGen.GetGradientDirections());
  model->m_CompartmentId = compartmentId;
  return model;
}

void QmitkFiberfoxView::SetGradientsFromGui()
{
  if (!m_LoadedGradients.empty())
  {
    // the gradients of a loaded .ffp, until a gradient setting is changed
    m_Parameters.SetGradienDirections(m_LoadedGradients);
    m_Parameters.SetBvalue(m_LoadedBvalue);
  }
  else if (m_Controls->m_UseBvalsBvecsBox->isChecked())
  {
    double bval;
    m_Parameters.SetGradienDirections( mitk::gradients::ReadBvalsBvecs(m_Controls->m_LoadBvalsEdit->text().toStdString(), m_Controls->m_LoadBvecsEdit->text().toStdString(), bval) );
    m_Parameters.SetBvalue(bval);
  }
  else
  {
    m_Parameters.SetNumWeightedVolumes(m_Controls->m_NumGradientsBox->value());
    m_Parameters.SetBvalue(m_Controls->m_BvalueBox->value());
    m_Parameters.GenerateGradientHalfShell();
  }
}

void QmitkFiberfoxView::OnGradientSettingsChanged()
{
  m_LoadedGradients.clear();
  m_LoadedPrototypes.clear();   // their kernels belong to the file's shell
}

void QmitkFiberfoxView::UpdateParametersFromGui()
{
  m_Parameters.ClearSignalParameters();
  ReadOutputSettings();
  ReadImageGeometry();
  ReadArtifacts();
  ReadSequence();
  ReadNoise();
  ReadCompartments();
  SetResultProperties();
}

void QmitkFiberfoxView::ReadOutputSettings()
{
  m_Parameters.m_Misc.m_CheckAdvancedSignalOptionsBox = m_Controls->m_AdvancedOptionsBox_2->isChecked();
  m_Parameters.m_Misc.m_OutputAdditionalImages = m_Controls->m_VolumeFractionsBox->isChecked();

  // the log file is named by the result, like the image and parameters saved with it; the prefix
  // names the output of MitkFiberfox
  m_Parameters.m_Misc.m_OutputPrefix = "";
  std::string outputPath = m_Controls->m_SavePathEdit->text().toStdString();
  // "-" or an emptied box: no automatic save (an empty path plus "/" would be the root folder)
  if (!outputPath.empty() && outputPath.compare("-")!=0)
  {
    m_Parameters.m_Misc.m_OutputPath = outputPath;
    if (outputPath.back()!='/')
      m_Parameters.m_Misc.m_OutputPath += "/";
  }
  else {
    m_Parameters.m_Misc.m_OutputPath = "";
  }
}

void QmitkFiberfoxView::ReadImageGeometry()
{
  if (m_Controls->m_MaskComboBox->GetSelectedNode().IsNotNull())
  {
    mitk::Image::ConstPointer mitkMaskImage = mitk::GetMaskImage(m_Controls->m_MaskComboBox->GetSelectedNode()->GetData());
    mitk::CastToItkImage<ItkUcharImgType>(mitkMaskImage, m_Parameters.m_SignalGen.m_MaskImage);
    itk::ImageDuplicator<ItkUcharImgType>::Pointer duplicator = itk::ImageDuplicator<ItkUcharImgType>::New();
    duplicator->SetInputImage(m_Parameters.m_SignalGen.m_MaskImage);
    duplicator->Update();
    m_Parameters.m_SignalGen.m_MaskImage = duplicator->GetOutput();
  }

  if (m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull() && mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( m_Controls->m_TemplateComboBox->GetSelectedNode()))  // use parameters of selected DWI
  {
    mitk::Image::Pointer dwi = dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData());

    ItkDwiType::Pointer itkVectorImagePointer = ItkDwiType::New();
    mitk::CastToItkImage(dwi, itkVectorImagePointer);

    m_Parameters.m_SignalGen.m_ImageRegion = itkVectorImagePointer->GetLargestPossibleRegion();
    m_Parameters.m_SignalGen.m_ImageSpacing = itkVectorImagePointer->GetSpacing();
    m_Parameters.m_SignalGen.m_ImageOrigin = itkVectorImagePointer->GetOrigin();
    m_Parameters.m_SignalGen.m_ImageDirection = itkVectorImagePointer->GetDirection();
    // as FSL bvecs: the view applies the direction matrix before the simulation
    m_Parameters.SetGradientsFromImage(dwi);
  }
  else if (m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull())   // use geometry of selected image
  {
    mitk::Image::Pointer img = dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData());
    itk::Image< float, 3 >::Pointer itkImg = itk::Image< float, 3 >::New();
    CastToItkImage< itk::Image< float, 3 > >(img, itkImg);

    m_Parameters.m_SignalGen.m_ImageRegion = itkImg->GetLargestPossibleRegion();
    m_Parameters.m_SignalGen.m_ImageSpacing = itkImg->GetSpacing();
    m_Parameters.m_SignalGen.m_ImageOrigin = itkImg->GetOrigin();
    m_Parameters.m_SignalGen.m_ImageDirection = itkImg->GetDirection();

    SetGradientsFromGui();
  }
  else if (m_Parameters.m_SignalGen.m_MaskImage.IsNotNull())   // use geometry of mask image
  {
    ItkUcharImgType::Pointer itkImg = m_Parameters.m_SignalGen.m_MaskImage;
    m_Parameters.m_SignalGen.m_ImageRegion = itkImg->GetLargestPossibleRegion();
    m_Parameters.m_SignalGen.m_ImageSpacing = itkImg->GetSpacing();
    m_Parameters.m_SignalGen.m_ImageOrigin = itkImg->GetOrigin();
    m_Parameters.m_SignalGen.m_ImageDirection = itkImg->GetDirection();

    SetGradientsFromGui();
  }
  else    // use GUI parameters
  {
    m_Parameters.m_SignalGen.m_ImageRegion.SetSize(0, m_Controls->m_SizeX->value());
    m_Parameters.m_SignalGen.m_ImageRegion.SetSize(1, m_Controls->m_SizeY->value());
    m_Parameters.m_SignalGen.m_ImageRegion.SetSize(2, m_Controls->m_SizeZ->value());
    m_Parameters.m_SignalGen.m_ImageSpacing[0] = m_Controls->m_SpacingX->value();
    m_Parameters.m_SignalGen.m_ImageSpacing[1] = m_Controls->m_SpacingY->value();
    m_Parameters.m_SignalGen.m_ImageSpacing[2] = m_Controls->m_SpacingZ->value();
    m_Parameters.m_SignalGen.m_ImageOrigin[0] = m_Parameters.m_SignalGen.m_ImageSpacing[0]/2;
    m_Parameters.m_SignalGen.m_ImageOrigin[1] = m_Parameters.m_SignalGen.m_ImageSpacing[1]/2;
    m_Parameters.m_SignalGen.m_ImageOrigin[2] = m_Parameters.m_SignalGen.m_ImageSpacing[2]/2;
    m_Parameters.m_SignalGen.m_ImageDirection.SetIdentity();

    SetGradientsFromGui();
  }
}

void QmitkFiberfoxView::ReadArtifacts()
{
  // signal relaxation
  m_Parameters.m_SignalGen.m_DoSimulateRelaxation = false;
  if (m_Controls->m_RelaxationBox->isChecked())
  {
    m_Parameters.m_SignalGen.m_DoSimulateRelaxation = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_RELAX";
  }
  m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = m_Parameters.m_SignalGen.m_DoSimulateRelaxation;

  // N/2 ghosts
  m_Parameters.m_Misc.m_DoAddGhosts = m_Controls->m_AddGhosts->isChecked();
  m_Parameters.m_SignalGen.m_KspaceLineOffset = m_Controls->m_kOffsetBox->value();
  if (m_Controls->m_AddGhosts->isChecked())
  {
    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_GHOST";
  }

  // Aliasing
  m_Parameters.m_Misc.m_DoAddAliasing = m_Controls->m_AddAliasing->isChecked();
  m_Parameters.m_SignalGen.m_CroppingFactor = (100-m_Controls->m_WrapBox->value())/100;
  if (m_Controls->m_AddAliasing->isChecked())
  {
    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_ALIASING";
  }

  // Spikes
  m_Parameters.m_Misc.m_DoAddSpikes = m_Controls->m_AddSpikes->isChecked();
  m_Parameters.m_SignalGen.m_Spikes = m_Controls->m_SpikeNumBox->value();
  m_Parameters.m_SignalGen.m_SpikeAmplitude = m_Controls->m_SpikeScaleBox->value();
  if (m_Controls->m_AddSpikes->isChecked())
  {
    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_SPIKES";
  }

  // Drift
  m_Parameters.m_SignalGen.m_DoAddDrift = m_Controls->m_AddDrift->isChecked();
  m_Parameters.m_SignalGen.m_Drift = static_cast<float>(m_Controls->m_DriftFactor->value())/100;
  if (m_Controls->m_AddDrift->isChecked())
  {
    m_Parameters.m_Misc.m_ArtifactModelString += "_DRIFT";
  }

  // gibbs ringing
  m_Parameters.m_SignalGen.m_DoAddGibbsRinging = m_Controls->m_AddGibbsRinging->isChecked();
  m_Parameters.m_SignalGen.m_ZeroRinging = m_Controls->m_ZeroRinging->value();
  if (m_Controls->m_AddGibbsRinging->isChecked())
  {
    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_RINGING";
  }

  // add distortions
  m_Parameters.m_Misc.m_DoAddDistortions = m_Controls->m_AddDistortions->isChecked();
  if (m_Controls->m_AddDistortions->isChecked() && m_Controls->m_FrequencyMapBox->GetSelectedNode().IsNotNull())
  {
    mitk::DataNode::Pointer fMapNode = m_Controls->m_FrequencyMapBox->GetSelectedNode();
    mitk::Image* img = dynamic_cast<mitk::Image*>(fMapNode->GetData());
    ItkFloatImgType::Pointer itkImg = ItkFloatImgType::New();
    CastToItkImage< ItkFloatImgType >(img, itkImg);

    if (m_Controls->m_TemplateComboBox->GetSelectedNode().IsNull())   // use geometry of frequency map
    {
      m_Parameters.m_SignalGen.m_ImageRegion = itkImg->GetLargestPossibleRegion();
      m_Parameters.m_SignalGen.m_ImageSpacing = itkImg->GetSpacing();
      m_Parameters.m_SignalGen.m_ImageOrigin = itkImg->GetOrigin();
      m_Parameters.m_SignalGen.m_ImageDirection = itkImg->GetDirection();
    }

    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    itk::ImageDuplicator<ItkFloatImgType>::Pointer duplicator = itk::ImageDuplicator<ItkFloatImgType>::New();
    duplicator->SetInputImage(itkImg);
    duplicator->Update();
    m_Parameters.m_SignalGen.m_FrequencyMap = duplicator->GetOutput();
    m_Parameters.m_Misc.m_ArtifactModelString += "_DISTORTED";
  }

  m_Parameters.m_SignalGen.m_EddyStrength = m_Controls->m_EddyGradientStrength->value();
  m_Parameters.m_Misc.m_DoAddEddyCurrents = m_Controls->m_AddEddy->isChecked();
  if (m_Controls->m_AddEddy->isChecked())
  {
    m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_EDDY";
  }

  // phase of distortions and eddy currents, read also while greyed out, so a saved .ffp keeps the choice
  switch (m_Controls->m_KspacePhaseBox->currentIndex())
  {
  case 0:
    m_Parameters.m_SignalGen.m_KspaceTimeSegments = 0;   // exact sum
    break;
  case 1:
    m_Parameters.m_SignalGen.m_KspaceTimeSegments = mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS;
    break;
  default:
    m_Parameters.m_SignalGen.m_KspaceTimeSegments = m_Controls->m_KspaceTimeSegmentsBox->value();
  }

  // Motion
  m_Parameters.m_SignalGen.m_DoAddMotion = false;
  m_Parameters.m_SignalGen.m_DoRandomizeMotion = m_Controls->m_RandomMotion->isChecked();
  m_Parameters.m_SignalGen.m_Translation[0] = m_Controls->m_MaxTranslationBoxX->value();
  m_Parameters.m_SignalGen.m_Translation[1] = m_Controls->m_MaxTranslationBoxY->value();
  m_Parameters.m_SignalGen.m_Translation[2] = m_Controls->m_MaxTranslationBoxZ->value();
  m_Parameters.m_SignalGen.m_Rotation[0] = m_Controls->m_MaxRotationBoxX->value();
  m_Parameters.m_SignalGen.m_Rotation[1] = m_Controls->m_MaxRotationBoxY->value();
  m_Parameters.m_SignalGen.m_Rotation[2] = m_Controls->m_MaxRotationBoxZ->value();
  m_Parameters.m_SignalGen.m_MotionVolumes.clear();
  // trimmed: spaces alone count as blank, " random" as random
  m_Parameters.m_Misc.m_MotionVolumesBox = m_Controls->m_MotionVolumesBox->text().trimmed().toStdString();

  if ( m_Controls->m_AddMotion->isChecked())
  {
    m_Parameters.m_SignalGen.m_DoAddMotion = true;
    m_Parameters.m_Misc.m_ArtifactModelString += "_MOTION";

    auto randgen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    randgen->SetSeed();
    if (!m_Parameters.m_SignalGen.SetMotionVolumes(m_Parameters.m_Misc.m_MotionVolumesBox, randgen))
    {
      MITK_WARN << "Cannot read the motion volumes '" << m_Parameters.m_Misc.m_MotionVolumesBox << "'. No volume is moved.";
      QMessageBox::warning(nullptr, "Motion volumes", "Cannot read the motion volumes \""
                           + m_Controls->m_MotionVolumesBox->text()
                           + "\": give \"random\", the indices of the moved volumes, or the negative indices of the volumes that stay. No volume is moved.");
    }
  }
}

void QmitkFiberfoxView::ReadSequence()
{
  // other imaging parameters
  m_Parameters.m_SignalGen.m_AcquisitionType = (mitk::SignalGenerationParameters::AcquisitionType)m_Controls->m_AcquisitionTypeBox->currentIndex();
  m_Parameters.m_SignalGen.m_CoilSensitivityProfile = (mitk::SignalGenerationParameters::CoilSensitivityProfile)m_Controls->m_CoilSensBox->currentIndex();
  m_Parameters.m_SignalGen.m_NumberOfCoils = m_Controls->m_NumCoilsBox->value();
  m_Parameters.m_SignalGen.m_PartialFourier = m_Controls->m_PartialFourier->value();
  m_Parameters.m_SignalGen.m_ReversePhase = m_Controls->m_ReversePhaseBox->isChecked();
  m_Parameters.m_SignalGen.m_tLine = m_Controls->m_LineReadoutTimeBox->value();
  m_Parameters.m_SignalGen.m_tInhom = m_Controls->m_T2starBox->value();
  m_Parameters.m_SignalGen.m_EchoTrainLength = m_Controls->m_EtlBox->value();
  m_Parameters.m_SignalGen.m_tEcho = m_Controls->m_TEbox->value();
  m_Parameters.m_SignalGen.m_tRep = m_Controls->m_TRbox->value();
  m_Parameters.m_SignalGen.m_tInv = m_Controls->m_TIbox->value();
  m_Parameters.m_SignalGen.m_DoDisablePartialVolume = m_Controls->m_EnforcePureFiberVoxelsBox->isChecked();
  m_Parameters.m_SignalGen.m_AxonRadius = m_Controls->m_FiberRadius->value();
  m_Parameters.m_SignalGen.m_SignalScale = m_Controls->m_SignalScaleBox->value();

//  double voxelVolume = m_Parameters.m_SignalGen.m_ImageSpacing[0]
//      * m_Parameters.m_SignalGen.m_ImageSpacing[1]
//      * m_Parameters.m_SignalGen.m_ImageSpacing[2];
//  if ( m_Parameters.m_SignalGen.m_SignalScale*voxelVolume > itk::NumericTraits<short>::max()*0.75 )
//  {
//    m_Parameters.m_SignalGen.m_SignalScale = itk::NumericTraits<short>::max()*0.75/voxelVolume;
//    m_Controls->m_SignalScaleBox->setValue(m_Parameters.m_SignalGen.m_SignalScale);
//    QMessageBox::information( nullptr, "Warning",
//                              "Maximum signal exceeding data type limits. Automatically adjusted to "
//                              + QString::number(m_Parameters.m_SignalGen.m_SignalScale)
//                              + " to obtain a maximum  signal of 75% of the data type maximum."
//                                " Relaxation and other effects that affect the signal intensities are not accounted for.");
//  }
}

void QmitkFiberfoxView::ReadNoise()
{
  // Noise
  m_Parameters.m_Misc.m_DoAddNoise = m_Controls->m_AddNoise->isChecked();
  m_Parameters.m_SignalGen.m_NoiseVariance = m_Controls->m_NoiseLevel->value();
  if (m_Controls->m_AddNoise->isChecked())
  {
    switch (m_Controls->m_NoiseDistributionBox->currentIndex())
    {
    case 0:
    {
      if (m_Parameters.m_SignalGen.m_NoiseVariance>0)
      {
        m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
        m_Parameters.m_Misc.m_ArtifactModelString += "_COMPLEX-GAUSSIAN-";
      }
      break;
    }
    case 1:
    {
      if (m_Parameters.m_SignalGen.m_NoiseVariance>0)
      {
        m_Parameters.m_NoiseModel = std::make_shared< mitk::RicianNoiseModel<mitk::ScalarType> >();
        m_Parameters.m_Misc.m_ArtifactModelString += "_RICIAN-";
        m_Parameters.m_NoiseModel->SetNoiseVariance(m_Parameters.m_SignalGen.m_NoiseVariance);
      }
      break;
    }
    case 2:
    {
      if (m_Parameters.m_SignalGen.m_NoiseVariance>0)
      {
        // The magnitude noise of "Channels" coils, without simulating each coil. Its own box, not
        // "Number of Channels" (the coils simulated in k-space): .ffp files keep both counts.
        auto chi = std::make_shared< mitk::NoncentralChiNoiseModel<mitk::ScalarType> >();
        chi->SetNumberOfCoils(m_Controls->m_NoiseCoilsBox->value());
        m_Parameters.m_NoiseModel = chi;
        m_Parameters.m_Misc.m_ArtifactModelString += "_NONCENTRAL-CHI-";
        m_Parameters.m_NoiseModel->SetNoiseVariance(m_Parameters.m_SignalGen.m_NoiseVariance);
      }
      break;
    }
    default:
    {
      if (m_Parameters.m_SignalGen.m_NoiseVariance>0)
      {
        m_Parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
        m_Parameters.m_Misc.m_ArtifactModelString += "_COMPLEX-GAUSSIAN-";
      }
      break;
    }
    }

    if (m_Parameters.m_SignalGen.m_NoiseVariance>0)
    {
      m_Parameters.m_Misc.m_ArtifactModelString += QString::number(m_Parameters.m_SignalGen.m_NoiseVariance).toStdString();
    }
  }
}

std::vector<QmitkFiberfoxView::CompartmentControls> QmitkFiberfoxView::Compartments() const
{
  using M = CompartmentModel;
  return {
    {1, "Intra-axonal compartment", m_Controls->m_Compartment1Box, m_Controls->m_Comp1VolumeFraction, nullptr,
     {{M::Stick, m_Controls->m_StickWidget1}, {M::Zeppelin, m_Controls->m_ZeppelinWidget1},
      {M::Tensor, m_Controls->m_TensorWidget1}, {M::Prototype, m_Controls->m_PrototypeWidget1}}},
    {2, "Inter-axonal compartment", m_Controls->m_Compartment2Box, m_Controls->m_Comp2VolumeFraction, m_Controls->m_Comp2FractionFrame,
     {{M::None, nullptr}, {M::Stick, m_Controls->m_StickWidget2}, {M::Zeppelin, m_Controls->m_ZeppelinWidget2},
      {M::Tensor, m_Controls->m_TensorWidget2}}},
    {3, "Extra-axonal compartment 1", m_Controls->m_Compartment3Box, m_Controls->m_Comp3VolumeFraction, nullptr,
     {{M::Ball, m_Controls->m_BallWidget1}, {M::Astrosticks, m_Controls->m_AstrosticksWidget1},
      {M::Dot, m_Controls->m_DotWidget1}, {M::Prototype, m_Controls->m_PrototypeWidget3}}},
    {4, "Extra-axonal compartment 2", m_Controls->m_Compartment4Box, m_Controls->m_Comp4VolumeFraction, m_Controls->m_Comp4FractionFrame,
     {{M::None, nullptr}, {M::Ball, m_Controls->m_BallWidget2}, {M::Astrosticks, m_Controls->m_AstrosticksWidget2},
      {M::Dot, m_Controls->m_DotWidget2}, {M::Prototype, m_Controls->m_PrototypeWidget4}}}
  };
}

std::string QmitkFiberfoxView::ModelName(CompartmentModel type)
{
  switch (type)
  {
  case CompartmentModel::Stick:       return "Stick";
  case CompartmentModel::Zeppelin:    return "Zeppelin";
  case CompartmentModel::Tensor:      return "Tensor";
  case CompartmentModel::Ball:        return "Ball";
  case CompartmentModel::Astrosticks: return "Astrosticks";
  case CompartmentModel::Dot:         return "Dot";
  case CompartmentModel::Prototype:   return "Prototype";
  case CompartmentModel::None:        break;
  }
  return "";
}

QmitkFiberfoxView::CompartmentModel QmitkFiberfoxView::SelectedModel(const CompartmentControls& compartment)
{
  const int index = compartment.box->currentIndex();
  if (index<0 || index>=static_cast<int>(compartment.entries.size()))
    return CompartmentModel::None;
  return compartment.entries[index].first;
}

QWidget* QmitkFiberfoxView::ModelWidget(const CompartmentControls& compartment, CompartmentModel type)
{
  for (const auto& entry : compartment.entries)
    if (entry.first==type)
      return entry.second;
  return nullptr;
}

mitk::DiffusionSignalModel<>* QmitkFiberfoxView::NewCompartmentModel(const CompartmentControls& compartment)
{
  const CompartmentModel type = SelectedModel(compartment);
  QWidget* widget = ModelWidget(compartment, type);
  mitk::DiffusionSignalModel<>* model = nullptr;
  switch (type)
  {
  case CompartmentModel::None:
    return nullptr;
  case CompartmentModel::Prototype:
    return NewPrototypeModel(static_cast<QmitkPrototypeSignalParametersWidget*>(widget), compartment.id);
  case CompartmentModel::Stick:
  {
    auto parameters = static_cast<QmitkStickModelParametersWidget*>(widget);
    auto stick = new mitk::StickModel<>();
    stick->SetBvalue(m_Parameters.m_SignalGen.GetBvalue());
    stick->SetDiffusivity(parameters->GetD());
    stick->SetT2(parameters->GetT2());
    stick->SetT1(parameters->GetT1());
    model = stick;
    break;
  }
  case CompartmentModel::Zeppelin:
  {
    // a tensor with the second and third diffusivity equal
    auto parameters = static_cast<QmitkZeppelinModelParametersWidget*>(widget);
    auto tensor = new mitk::TensorModel<>();
    tensor->SetBvalue(m_Parameters.m_SignalGen.GetBvalue());
    tensor->SetDiffusivity1(parameters->GetD1());
    tensor->SetDiffusivity2(parameters->GetD2());
    tensor->SetDiffusivity3(parameters->GetD2());
    tensor->SetT2(parameters->GetT2());
    tensor->SetT1(parameters->GetT1());
    model = tensor;
    break;
  }
  case CompartmentModel::Tensor:
  {
    auto parameters = static_cast<QmitkTensorModelParametersWidget*>(widget);
    auto tensor = new mitk::TensorModel<>();
    tensor->SetBvalue(m_Parameters.m_SignalGen.GetBvalue());
    tensor->SetDiffusivity1(parameters->GetD1());
    tensor->SetDiffusivity2(parameters->GetD2());
    tensor->SetDiffusivity3(parameters->GetD3());
    tensor->SetT2(parameters->GetT2());
    tensor->SetT1(parameters->GetT1());
    model = tensor;
    break;
  }
  case CompartmentModel::Ball:
  {
    auto parameters = static_cast<QmitkBallModelParametersWidget*>(widget);
    auto ball = new mitk::BallModel<>();
    ball->SetBvalue(m_Parameters.m_SignalGen.GetBvalue());
    ball->SetDiffusivity(parameters->GetD());
    ball->SetT2(parameters->GetT2());
    ball->SetT1(parameters->GetT1());
    model = ball;
    break;
  }
  case CompartmentModel::Astrosticks:
  {
    auto parameters = static_cast<QmitkAstrosticksModelParametersWidget*>(widget);
    auto astrosticks = new mitk::AstroStickModel<>();
    astrosticks->SetBvalue(m_Parameters.m_SignalGen.GetBvalue());
    astrosticks->SetDiffusivity(parameters->GetD());
    astrosticks->SetT2(parameters->GetT2());
    astrosticks->SetT1(parameters->GetT1());
    astrosticks->SetRandomizeSticks(parameters->GetRandomizeSticks());
    model = astrosticks;
    break;
  }
  case CompartmentModel::Dot:
  {
    auto parameters = static_cast<QmitkDotModelParametersWidget*>(widget);
    auto dot = new mitk::DotModel<>();
    dot->SetT2(parameters->GetT2());
    dot->SetT1(parameters->GetT1());
    model = dot;
    break;
  }
  }
  model->SetGradientList(m_Parameters.m_SignalGen.GetGradientDirections());
  model->m_CompartmentId = compartment.id;
  return model;
}

void QmitkFiberfoxView::ShowCompartmentModel(const CompartmentControls& compartment, mitk::DiffusionSignalModel<>* model, bool zeppelin)
{
  CompartmentModel type = CompartmentModel::None;
  if (dynamic_cast<mitk::StickModel<>*>(model))
    type = CompartmentModel::Stick;
  else if (dynamic_cast<mitk::TensorModel<>*>(model))
    type = zeppelin && ModelWidget(compartment, CompartmentModel::Zeppelin)!=nullptr ? CompartmentModel::Zeppelin : CompartmentModel::Tensor;
  else if (dynamic_cast<mitk::BallModel<>*>(model))
    type = CompartmentModel::Ball;
  else if (dynamic_cast<mitk::AstroStickModel<>*>(model))
    type = CompartmentModel::Astrosticks;
  else if (dynamic_cast<mitk::DotModel<>*>(model))
    type = CompartmentModel::Dot;
  else if (dynamic_cast<mitk::RawShModel<>*>(model))
    type = CompartmentModel::Prototype;

  int index = -1;
  for (std::size_t i=0; i<compartment.entries.size(); ++i)
    if (compartment.entries[i].first==type)
      index = static_cast<int>(i);
  if (index<0 || type==CompartmentModel::None)
    return;   // the compartment does not offer this model

  QWidget* widget = compartment.entries[index].second;
  switch (type)
  {
  case CompartmentModel::Stick:
  {
    auto parameters = static_cast<QmitkStickModelParametersWidget*>(widget);
    auto stick = static_cast<mitk::StickModel<>*>(model);
    parameters->SetT2(stick->GetT2());
    parameters->SetT1(stick->GetT1());
    parameters->SetD(stick->GetDiffusivity());
    break;
  }
  case CompartmentModel::Zeppelin:
  {
    auto parameters = static_cast<QmitkZeppelinModelParametersWidget*>(widget);
    auto tensor = static_cast<mitk::TensorModel<>*>(model);
    parameters->SetT2(tensor->GetT2());
    parameters->SetT1(tensor->GetT1());
    parameters->SetD1(tensor->GetDiffusivity1());
    parameters->SetD2(tensor->GetDiffusivity2());
    break;
  }
  case CompartmentModel::Tensor:
  {
    auto parameters = static_cast<QmitkTensorModelParametersWidget*>(widget);
    auto tensor = static_cast<mitk::TensorModel<>*>(model);
    parameters->SetT2(tensor->GetT2());
    parameters->SetT1(tensor->GetT1());
    parameters->SetD1(tensor->GetDiffusivity1());
    parameters->SetD2(tensor->GetDiffusivity2());
    parameters->SetD3(tensor->GetDiffusivity3());
    break;
  }
  case CompartmentModel::Ball:
  {
    auto parameters = static_cast<QmitkBallModelParametersWidget*>(widget);
    auto ball = static_cast<mitk::BallModel<>*>(model);
    parameters->SetT2(ball->GetT2());
    parameters->SetT1(ball->GetT1());
    parameters->SetD(ball->GetDiffusivity());
    break;
  }
  case CompartmentModel::Astrosticks:
  {
    auto parameters = static_cast<QmitkAstrosticksModelParametersWidget*>(widget);
    auto astrosticks = static_cast<mitk::AstroStickModel<>*>(model);
    parameters->SetT2(astrosticks->GetT2());
    parameters->SetT1(astrosticks->GetT1());
    parameters->SetD(astrosticks->GetDiffusivity());
    parameters->SetRandomizeSticks(astrosticks->GetRandomizeSticks());
    break;
  }
  case CompartmentModel::Dot:
  {
    auto parameters = static_cast<QmitkDotModelParametersWidget*>(widget);
    parameters->SetT2(model->GetT2());
    parameters->SetT1(model->GetT1());
    break;
  }
  case CompartmentModel::Prototype:
  {
    auto parameters = static_cast<QmitkPrototypeSignalParametersWidget*>(widget);
    auto prototype = static_cast<mitk::RawShModel<>*>(model);
    parameters->SetNumberOfSamples(prototype->GetMaxNumKernels());
    parameters->SetMinFa(prototype->GetFaRange().first);
    parameters->SetMaxFa(prototype->GetFaRange().second);
    parameters->SetMinAdc(prototype->GetAdcRange().first);
    parameters->SetMaxAdc(prototype->GetAdcRange().second);
    break;
  }
  case CompartmentModel::None:
    break;
  }
  compartment.box->setCurrentIndex(index);
}

void QmitkFiberfoxView::ShowCompartmentWidget(int compartmentId, int index)
{
  const auto compartments = Compartments();
  const CompartmentControls& compartment = compartments.at(compartmentId-1);
  for (const auto& entry : compartment.entries)
    if (entry.second!=nullptr)
      entry.second->setVisible(false);

  const bool valid = index>=0 && index<static_cast<int>(compartment.entries.size());
  if (valid && compartment.entries[index].second!=nullptr)
    compartment.entries[index].second->setVisible(true);
  if (compartment.fractionFrame!=nullptr)
    compartment.fractionFrame->setVisible(valid && compartment.entries[index].first!=CompartmentModel::None);
}

void QmitkFiberfoxView::ReadCompartments()
{
  for (const auto& compartment : Compartments())
  {
    const CompartmentModel type = SelectedModel(compartment);
    mitk::DiffusionSignalModel<>* model = NewCompartmentModel(compartment);
    if (model==nullptr)
      continue;
    if (type==CompartmentModel::Prototype)
      m_Parameters.m_SignalGen.m_DoSimulateRelaxation = false;

    // compartments 1 and 2 hold the fibers, 3 and 4 the rest of the voxel
    if (compartment.id<=2)
      m_Parameters.m_FiberModelList.push_back(model);
    else
      m_Parameters.m_NonFiberModelList.push_back(model);
    m_Parameters.m_Misc.m_SignalModelString += ModelName(type);

    if (compartment.volumeFraction->GetSelectedNode().IsNotNull())
    {
      ItkDoubleImgType::Pointer volumeImage = ItkDoubleImgType::New();
      mitk::Image* img = dynamic_cast<mitk::Image*>(compartment.volumeFraction->GetSelectedNode()->GetData());
      CastToItkImage< ItkDoubleImgType >(img, volumeImage);
      model->SetVolumeFractionImage(volumeImage);
    }
  }
}

void QmitkFiberfoxView::SetResultProperties()
{
  auto& signalGen = m_Parameters.m_SignalGen;   // GetBvalue is not const
  const auto& misc = m_Parameters.m_Misc;
  mitk::DataNode* node = misc.m_ResultNode;

  node->SetProperty("Fiberfox.Relaxation", mitk::BoolProperty::New(signalGen.m_DoSimulateRelaxation));
  if (misc.m_DoAddGhosts)
    node->SetProperty("Fiberfox.Ghost", mitk::DoubleProperty::New(signalGen.m_KspaceLineOffset));
  if (misc.m_DoAddAliasing)
    node->SetProperty("Fiberfox.Aliasing", mitk::DoubleProperty::New(m_Controls->m_WrapBox->value()));  // in percent, as in the view
  if (misc.m_DoAddSpikes)
  {
    node->SetProperty("Fiberfox.Spikes.Number", mitk::IntProperty::New(signalGen.m_Spikes));
    node->SetProperty("Fiberfox.Spikes.Amplitude", mitk::DoubleProperty::New(signalGen.m_SpikeAmplitude));
  }
  if (signalGen.m_DoAddDrift)
    node->SetProperty("Fiberfox.Drift", mitk::FloatProperty::New(signalGen.m_Drift));
  if (signalGen.m_DoAddGibbsRinging)
    node->SetProperty("Fiberfox.Ringing", mitk::BoolProperty::New(true));
  if (misc.m_DoAddDistortions && signalGen.m_FrequencyMap.IsNotNull())
    node->SetProperty("Fiberfox.Distortions", mitk::BoolProperty::New(true));
  if (misc.m_DoAddEddyCurrents)
    node->SetProperty("Fiberfox.Eddy-strength", mitk::DoubleProperty::New(signalGen.m_EddyStrength));
  if ((misc.m_DoAddDistortions && signalGen.m_FrequencyMap.IsNotNull()) || misc.m_DoAddEddyCurrents)
  {
    const int segments = signalGen.m_KspaceTimeSegments;
    const std::string phase = segments==0 ? "Exact" : segments<0 ? "Time segments (automatic)" : std::to_string(segments)+" time segments";
    node->SetProperty("Fiberfox.Distortion-phase", mitk::StringProperty::New(phase));
  }
  if (signalGen.m_DoAddMotion)
  {
    node->SetProperty("Fiberfox.Motion.Random", mitk::BoolProperty::New(signalGen.m_DoRandomizeMotion));
    node->SetProperty("Fiberfox.Motion.Translation-x", mitk::DoubleProperty::New(signalGen.m_Translation[0]));
    node->SetProperty("Fiberfox.Motion.Translation-y", mitk::DoubleProperty::New(signalGen.m_Translation[1]));
    node->SetProperty("Fiberfox.Motion.Translation-z", mitk::DoubleProperty::New(signalGen.m_Translation[2]));
    node->SetProperty("Fiberfox.Motion.Rotation-x", mitk::DoubleProperty::New(signalGen.m_Rotation[0]));
    node->SetProperty("Fiberfox.Motion.Rotation-y", mitk::DoubleProperty::New(signalGen.m_Rotation[1]));
    node->SetProperty("Fiberfox.Motion.Rotation-z", mitk::DoubleProperty::New(signalGen.m_Rotation[2]));
  }
  if (misc.m_DoAddNoise && signalGen.m_NoiseVariance>0)
  {
    std::string distribution = "Complex Gaussian";   // in k-space, without a noise model
    if (auto chi = dynamic_cast<mitk::NoncentralChiNoiseModel<mitk::ScalarType>*>(m_Parameters.m_NoiseModel.get()))
    {
      distribution = "Noncentral chi";
      node->SetProperty("Fiberfox.Noise-Channels", mitk::IntProperty::New(static_cast<int>(chi->GetNumberOfCoils())));
    }
    else if (m_Parameters.m_NoiseModel)
      distribution = "Rician";
    node->SetProperty("Fiberfox.Noise-Distribution", mitk::StringProperty::New(distribution));
    node->SetProperty("Fiberfox.Noise-Variance", mitk::DoubleProperty::New(signalGen.m_NoiseVariance));
  }

  for (const auto& compartment : Compartments())
  {
    const CompartmentModel type = SelectedModel(compartment);
    mitk::DiffusionSignalModel<>* model = nullptr;
    for (auto list : {&m_Parameters.m_FiberModelList, &m_Parameters.m_NonFiberModelList})
      for (auto m : *list)
        if (static_cast<int>(m->m_CompartmentId)==compartment.id)
          model = m;
    if (model==nullptr)
      continue;

    const std::string prefix = "Fiberfox.Compartment" + std::to_string(compartment.id) + ".";
    node->SetProperty(prefix + "Description", mitk::StringProperty::New(compartment.description));
    node->SetProperty(prefix + "Model", mitk::StringProperty::New(ModelName(type)));
    switch (type)
    {
    case CompartmentModel::Stick:
      node->SetProperty(prefix + "D", mitk::DoubleProperty::New(static_cast<mitk::StickModel<>*>(model)->GetDiffusivity()));
      break;
    case CompartmentModel::Ball:
      node->SetProperty(prefix + "D", mitk::DoubleProperty::New(static_cast<mitk::BallModel<>*>(model)->GetDiffusivity()));
      break;
    case CompartmentModel::Astrosticks:
      node->SetProperty(prefix + "D", mitk::DoubleProperty::New(static_cast<mitk::AstroStickModel<>*>(model)->GetDiffusivity()));
      node->SetProperty(prefix + "RandomSticks", mitk::BoolProperty::New(static_cast<mitk::AstroStickModel<>*>(model)->GetRandomizeSticks()));
      break;
    case CompartmentModel::Zeppelin:
    case CompartmentModel::Tensor:
    {
      auto tensor = static_cast<mitk::TensorModel<>*>(model);
      node->SetProperty(prefix + "D1", mitk::DoubleProperty::New(tensor->GetDiffusivity1()));
      node->SetProperty(prefix + "D2", mitk::DoubleProperty::New(tensor->GetDiffusivity2()));
      if (type==CompartmentModel::Tensor)
        node->SetProperty(prefix + "D3", mitk::DoubleProperty::New(tensor->GetDiffusivity3()));
      break;
    }
    case CompartmentModel::Dot:
    case CompartmentModel::Prototype:
    case CompartmentModel::None:
      break;
    }
    if (type!=CompartmentModel::Prototype)
      node->SetProperty(prefix + "T2", mitk::DoubleProperty::New(model->GetT2()));
  }

  node->SetProperty("Fiberfox.SignalScale", mitk::DoubleProperty::New(signalGen.m_SignalScale));
  node->SetProperty("Fiberfox.FiberRadius", mitk::DoubleProperty::New(signalGen.m_AxonRadius));
  node->SetProperty("Fiberfox.Tinhom", mitk::DoubleProperty::New(signalGen.m_tInhom));
  node->SetProperty("Fiberfox.echoTrainLength", mitk::IntProperty::New(signalGen.m_EchoTrainLength));
  node->SetProperty("Fiberfox.Tline", mitk::DoubleProperty::New(signalGen.m_tLine));
  node->SetProperty("Fiberfox.TE", mitk::DoubleProperty::New(signalGen.m_tEcho));
  node->SetProperty("Fiberfox.b-value", mitk::DoubleProperty::New(signalGen.GetBvalue()));
  node->SetProperty("Fiberfox.NoPartialVolume", mitk::BoolProperty::New(signalGen.m_DoDisablePartialVolume));
  node->SetProperty("binary", mitk::BoolProperty::New(false));
}

std::vector< mitk::RawShModel<>* > QmitkFiberfoxView::GetPrototypeModelsWithoutKernels()
{
  std::vector< mitk::RawShModel<>* > models;
  for (auto list : {&m_Parameters.m_FiberModelList, &m_Parameters.m_NonFiberModelList})
    for (auto signalModel : *list)
    {
      auto model = dynamic_cast< mitk::RawShModel<>* >(signalModel);
      if (model!=nullptr && model->GetNumberOfKernels()<=0)
        models.push_back(model);
    }
  return models;
}

bool QmitkFiberfoxView::SamplePrototypeKernels(mitk::Image::Pointer diffImg, const std::vector< mitk::RawShModel<>* >& models)
{
  if (models.empty())
    return true;

  // tensor, Q-ball and ADC images of the template select the prototype voxels; computed once for all models
  ItkDwiType::Pointer itkVectorImagePointer = ItkDwiType::New();
  mitk::CastToItkImage(diffImg, itkVectorImagePointer);
  const double bValue = mitk::DiffusionPropertyHelper::GetReferenceBValue(diffImg);

  // ITK takes only exact zero vectors as b0 volumes and fits every other volume at the reference
  // b-value, so it gets a copy in which the b0 volumes of MITK's rule are zero vectors, as the
  // Q-ball and ADC filters below and RawShModel::SampleKernels take them.
  auto gradientsForItk = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
  for (const auto& g : *mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg))
    gradientsForItk->push_back(mitk::DiffusionPropertyHelper::IsB0(g, bValue) ? mitk::DiffusionPropertyHelper::GradientDirectionType(0.0) : g);

  typedef itk::DiffusionTensor3DReconstructionImageFilter< short, short, double > TensorReconstructionImageFilterType;
  TensorReconstructionImageFilterType::Pointer filter = TensorReconstructionImageFilterType::New();
  filter->SetBValue(bValue);
  filter->SetGradientImage(gradientsForItk, itkVectorImagePointer);
  filter->Update();

  typedef mitk::RawShModel<>::QballFilterType QballFilterType;
  QballFilterType::Pointer qballfilter = QballFilterType::New();
  qballfilter->SetBValue(bValue);
  qballfilter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg),
                                itkVectorImagePointer );
  qballfilter->SetLambda(0.006);
  qballfilter->SetNormalizationMethod(QballFilterType::QBAR_RAW_SIGNAL);
  qballfilter->Update();

  itk::AdcImageFilter< short, double >::Pointer adcFilter = itk::AdcImageFilter< short, double >::New();
  adcFilter->SetInput( itkVectorImagePointer );
  adcFilter->SetGradientDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(diffImg));
  adcFilter->SetB_value(bValue);
  adcFilter->Update();

  for (auto model : models)
  {
    if (!model->SampleKernels(diffImg, m_Parameters.m_SignalGen.m_MaskImage, filter->GetOutput(), qballfilter->GetCoefficientImage(), adcFilter->GetOutput()))
      return false;
  }
  return true;
}

void QmitkFiberfoxView::SaveParameters(QString filename)
{
  UpdateParametersFromGui();

  // Declining the sampling saves the prototype settings without kernels; a simulation samples them
  // from its template image then.
  auto prototypes = GetPrototypeModelsWithoutKernels();
  if ( !prototypes.empty()
       && QMessageBox::question(nullptr, "Prototype signal sampling",
                                "Do you want to sample prototype signals from the selected diffusion-weighted image and save them?",
                                QMessageBox::Yes, QMessageBox::No) == QMessageBox::Yes )
  {
    if ( m_Controls->m_TemplateComboBox->GetSelectedNode().IsNull()
         || !mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(
           dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData()) ) )
    {
      QMessageBox::information(nullptr, "Parameter file not saved", "No diffusion-weighted image selected to sample signal from.");
      return;
    }
    if (!SamplePrototypeKernels(dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData()), prototypes))
    {
      QMessageBox::information( nullptr, "Parameter file not saved", "No valid prototype signals could be sampled.");
      return;
    }
  }

  m_Parameters.SaveParameters(filename.toStdString());
  m_ParameterFile = filename;
}

void QmitkFiberfoxView::SaveParameters()
{
  QString filename = QFileDialog::getSaveFileName(
        0,
        tr("Save Parameters"),
        m_ParameterFile,
        tr("Fiberfox Parameters (*.ffp)") );
  if (filename.isEmpty())
    return;

  try
  {
    SaveParameters(filename);
  }
  catch (const itk::ExceptionObject& e)
  {
    QMessageBox::warning(nullptr, "Parameter file not saved", e.GetDescription());
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Parameter file not saved", e.what());
  }
}

void QmitkFiberfoxView::LoadParameters()
{
  QString filename = QFileDialog::getOpenFileName(0, tr("Load Parameters"), QString(itksys::SystemTools::GetFilenamePath(m_ParameterFile.toStdString()).c_str()), tr("Fiberfox Parameters (*.ffp)") );
  if(filename.isEmpty() || filename.isNull())
    return;

  try
  {
    m_Parameters.LoadParameters(filename.toStdString());
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Parameters not loaded", QString("Cannot read ") + filename + ":\n" + e.what());
    return;
  }
  m_ParameterFile = filename;

  if (m_Parameters.m_MissingTags.size()>0)
  {
    QString missing("Parameter file might be corrupted. The following parameters could not be read: ");
    missing += QString(m_Parameters.m_MissingTags.c_str());
    missing += "\nDefault values have been assigned to the missing parameters.";
    QMessageBox::information( nullptr, "Warning!", missing);
  }

  // image generation parameters
  m_Controls->m_SizeX->setValue(m_Parameters.m_SignalGen.m_ImageRegion.GetSize(0));
  m_Controls->m_SizeY->setValue(m_Parameters.m_SignalGen.m_ImageRegion.GetSize(1));
  m_Controls->m_SizeZ->setValue(m_Parameters.m_SignalGen.m_ImageRegion.GetSize(2));
  m_Controls->m_SpacingX->setValue(m_Parameters.m_SignalGen.m_ImageSpacing[0]);
  m_Controls->m_SpacingY->setValue(m_Parameters.m_SignalGen.m_ImageSpacing[1]);
  m_Controls->m_SpacingZ->setValue(m_Parameters.m_SignalGen.m_ImageSpacing[2]);
  m_Controls->m_NumGradientsBox->setValue(m_Parameters.m_SignalGen.GetNumWeightedVolumes());
  m_Controls->m_BvalueBox->setValue(m_Parameters.m_SignalGen.GetBvalue());
  m_Controls->m_SignalScaleBox->setValue(m_Parameters.m_SignalGen.m_SignalScale);
  m_Controls->m_TEbox->setValue(m_Parameters.m_SignalGen.m_tEcho);
  m_Controls->m_LineReadoutTimeBox->setValue(m_Parameters.m_SignalGen.m_tLine);
  m_Controls->m_T2starBox->setValue(m_Parameters.m_SignalGen.m_tInhom);
  m_Controls->m_EtlBox->setValue(m_Parameters.m_SignalGen.m_EchoTrainLength);
  m_Controls->m_FiberRadius->setValue(m_Parameters.m_SignalGen.m_AxonRadius);
  m_Controls->m_RelaxationBox->setChecked(m_Parameters.m_SignalGen.m_DoSimulateRelaxation);
  m_Controls->m_EnforcePureFiberVoxelsBox->setChecked(m_Parameters.m_SignalGen.m_DoDisablePartialVolume);
  m_Controls->m_ReversePhaseBox->setChecked(m_Parameters.m_SignalGen.m_ReversePhase);
  m_Controls->m_PartialFourier->setValue(m_Parameters.m_SignalGen.m_PartialFourier);
  m_Controls->m_TRbox->setValue(m_Parameters.m_SignalGen.m_tRep);
  m_Controls->m_TIbox->setValue(m_Parameters.m_SignalGen.m_tInv);
  m_Controls->m_NumCoilsBox->setValue(m_Parameters.m_SignalGen.m_NumberOfCoils);
  m_Controls->m_CoilSensBox->setCurrentIndex(m_Parameters.m_SignalGen.m_CoilSensitivityProfile);
  m_Controls->m_AcquisitionTypeBox->setCurrentIndex(m_Parameters.m_SignalGen.m_AcquisitionType);

  if (!m_Parameters.m_Misc.m_BvalsFile.empty())
  {
    m_Controls->m_UseBvalsBvecsBox->setChecked(true);
    m_Controls->m_LoadBvalsEdit->setText(QString(m_Parameters.m_Misc.m_BvalsFile.c_str()));
  }
  else
    m_Controls->m_LoadBvalsEdit->setText("-");

  if (!m_Parameters.m_Misc.m_BvecsFile.empty())
  {
    m_Controls->m_UseBvalsBvecsBox->setChecked(true);
    m_Controls->m_LoadBvecsEdit->setText(QString(m_Parameters.m_Misc.m_BvecsFile.c_str()));
  }
  else
    m_Controls->m_LoadBvecsEdit->setText("-");

  // The view keeps the gradients as loaded (set after the boxes, whose changes discard them): without
  // side files they come from the .ffp itself and the boxes show only their count and b-value; with
  // them, reading the side files again would give them rounded (FiberfoxParameters::LoadParameters),
  // and saving would change the file.
  if (m_Parameters.m_Misc.m_BvalsFile.empty() || m_Parameters.m_Misc.m_BvecsFile.empty())
    m_Controls->m_UseBvalsBvecsBox->setChecked(false);
  m_LoadedGradients = m_Parameters.m_SignalGen.GetGradientDirections();
  m_LoadedBvalue = m_Parameters.m_SignalGen.GetBvalue();

  // m_NoiseDistributionBox offers 0 "Complex Gaussian" (no noise model), 1 "Rician" and 2 "Noncentral chi"
  if (m_Parameters.m_NoiseModel!=nullptr)
  {
    m_Controls->m_AddNoise->setChecked(m_Parameters.m_Misc.m_DoAddNoise);
    m_Controls->m_NoiseDistributionBox->setCurrentIndex(1);
    if (auto chi = dynamic_cast<mitk::NoncentralChiNoiseModel<double>*>(m_Parameters.m_NoiseModel.get()))
    {
      m_Controls->m_NoiseDistributionBox->setCurrentIndex(2);
      m_Controls->m_NoiseCoilsBox->setValue(chi->GetNumberOfCoils());
    }
    m_Controls->m_NoiseLevel->setValue(m_Parameters.m_NoiseModel->GetNoiseVariance());
  }
  else
  {
    m_Controls->m_AddNoise->setChecked(m_Parameters.m_Misc.m_DoAddNoise);
    m_Controls->m_NoiseDistributionBox->setCurrentIndex(0);
    m_Controls->m_NoiseLevel->setValue(m_Parameters.m_SignalGen.m_NoiseVariance);
  }

  m_Controls->m_VolumeFractionsBox->setChecked(m_Parameters.m_Misc.m_OutputAdditionalImages);
  m_Controls->m_AdvancedOptionsBox_2->setChecked(m_Parameters.m_Misc.m_CheckAdvancedSignalOptionsBox);
  m_Controls->m_AddGhosts->setChecked(m_Parameters.m_Misc.m_DoAddGhosts);
  m_Controls->m_AddAliasing->setChecked(m_Parameters.m_Misc.m_DoAddAliasing);
  m_Controls->m_AddDistortions->setChecked(m_Parameters.m_Misc.m_DoAddDistortions);
  m_Controls->m_AddSpikes->setChecked(m_Parameters.m_Misc.m_DoAddSpikes);
  m_Controls->m_AddEddy->setChecked(m_Parameters.m_Misc.m_DoAddEddyCurrents);
  m_Controls->m_AddDrift->setChecked(m_Parameters.m_SignalGen.m_DoAddDrift);

  m_Controls->m_kOffsetBox->setValue(m_Parameters.m_SignalGen.m_KspaceLineOffset);
  m_Controls->m_WrapBox->setValue(100*(1-m_Parameters.m_SignalGen.m_CroppingFactor));
  m_Controls->m_DriftFactor->setValue(100*m_Parameters.m_SignalGen.m_Drift);
  m_Controls->m_SpikeNumBox->setValue(m_Parameters.m_SignalGen.m_Spikes);
  m_Controls->m_SpikeScaleBox->setValue(m_Parameters.m_SignalGen.m_SpikeAmplitude);
  m_Controls->m_EddyGradientStrength->setValue(m_Parameters.m_SignalGen.m_EddyStrength);
  ShowKspaceTimeSegments(m_Parameters.m_SignalGen.m_KspaceTimeSegments);
  m_Controls->m_AddGibbsRinging->setChecked(m_Parameters.m_SignalGen.m_DoAddGibbsRinging);
  m_Controls->m_ZeroRinging->setValue(m_Parameters.m_SignalGen.m_ZeroRinging);
  m_Controls->m_AddMotion->setChecked(m_Parameters.m_SignalGen.m_DoAddMotion);
  m_Controls->m_RandomMotion->setChecked(m_Parameters.m_SignalGen.m_DoRandomizeMotion);
  m_Controls->m_MotionVolumesBox->setText(QString(m_Parameters.m_Misc.m_MotionVolumesBox.c_str()));

  m_Controls->m_MaxTranslationBoxX->setValue(m_Parameters.m_SignalGen.m_Translation[0]);
  m_Controls->m_MaxTranslationBoxY->setValue(m_Parameters.m_SignalGen.m_Translation[1]);
  m_Controls->m_MaxTranslationBoxZ->setValue(m_Parameters.m_SignalGen.m_Translation[2]);
  m_Controls->m_MaxRotationBoxX->setValue(m_Parameters.m_SignalGen.m_Rotation[0]);
  m_Controls->m_MaxRotationBoxY->setValue(m_Parameters.m_SignalGen.m_Rotation[1]);
  m_Controls->m_MaxRotationBoxZ->setValue(m_Parameters.m_SignalGen.m_Rotation[2]);

  m_Controls->m_Compartment1Box->setCurrentIndex(0);
  m_Controls->m_Compartment2Box->setCurrentIndex(0);
  m_Controls->m_Compartment3Box->setCurrentIndex(0);
  m_Controls->m_Compartment4Box->setCurrentIndex(0);

  const auto compartments = Compartments();
  // A zeppelin is saved as a tensor; the file's signal model string, one name per model in list
  // order, tells them apart, so a saved view loads back as it was.
  const auto modelNames = ModelNames(m_Parameters.m_Misc.m_SignalModelString);
  const bool namesMatch = modelNames.size()==m_Parameters.m_FiberModelList.size()+m_Parameters.m_NonFiberModelList.size();
  for (unsigned int i=0; i<m_Parameters.m_FiberModelList.size()+m_Parameters.m_NonFiberModelList.size(); i++)
  {
    mitk::DiffusionSignalModel<mitk::ScalarType>* signalModel = nullptr;
    if (i<m_Parameters.m_FiberModelList.size())
      signalModel = m_Parameters.m_FiberModelList.at(i);
    else
      signalModel = m_Parameters.m_NonFiberModelList.at(i-m_Parameters.m_FiberModelList.size());

    mitk::DataNode::Pointer compVolNode;
    if ( signalModel->GetVolumeFractionImage().IsNotNull() )
    {
      compVolNode = mitk::DataNode::New();
      mitk::Image::Pointer image = mitk::Image::New();
      image->InitializeByItk(signalModel->GetVolumeFractionImage().GetPointer());
      image->SetVolume(signalModel->GetVolumeFractionImage()->GetBufferPointer());

      compVolNode->SetData( image );
      compVolNode->SetName("Compartment volume "+QString::number(signalModel->m_CompartmentId).toStdString());
      GetDataStorage()->Add(compVolNode);
    }

    if (signalModel->m_CompartmentId<1 || signalModel->m_CompartmentId>compartments.size())
      continue;
    const CompartmentControls& compartment = compartments[signalModel->m_CompartmentId-1];
    if (compVolNode.IsNotNull())
      compartment.volumeFraction->SetSelectedNode(compVolNode);
    ShowCompartmentModel(compartment, signalModel, namesMatch && modelNames[i]=="Zeppelin");
  }

  // Prototype kernels of the file are reused by NewPrototypeModel while the widgets keep the settings
  // just set (read back, i.e. as the spin boxes round them).
  m_LoadedPrototypes.clear();
  for (auto list : {&m_Parameters.m_FiberModelList, &m_Parameters.m_NonFiberModelList})
    for (auto signalModel : *list)
    {
      auto model = dynamic_cast< mitk::RawShModel<>* >(signalModel);
      QmitkPrototypeSignalParametersWidget* widget = model!=nullptr ? GetPrototypeWidget(model->m_CompartmentId) : nullptr;
      if (widget!=nullptr && model->GetNumberOfKernels()>0)
      {
        LoadedPrototype& loaded = m_LoadedPrototypes[model->m_CompartmentId];
        loaded.m_Model.reset(model->Clone());
        loaded.m_WidgetSettings = PrototypeSettings(widget);
        loaded.m_Shells = Shells(m_Parameters.m_SignalGen);
      }
    }

  if ( m_Parameters.m_SignalGen.m_MaskImage )
  {
    mitk::Image::Pointer image = mitk::Image::New();
    image->InitializeByItk(m_Parameters.m_SignalGen.m_MaskImage.GetPointer());
    image->SetVolume(m_Parameters.m_SignalGen.m_MaskImage->GetBufferPointer());

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData( image );
    node->SetName("Tissue mask");
    GetDataStorage()->Add(node);
    m_Controls->m_MaskComboBox->SetSelectedNode(node);
  }

  if ( m_Parameters.m_SignalGen.m_FrequencyMap )
  {
    mitk::Image::Pointer image = mitk::Image::New();
    image->InitializeByItk(m_Parameters.m_SignalGen.m_FrequencyMap.GetPointer());
    image->SetVolume(m_Parameters.m_SignalGen.m_FrequencyMap->GetBufferPointer());

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData( image );
    node->SetName("Frequency map");
    GetDataStorage()->Add(node);
    m_Controls->m_FrequencyMapBox->SetSelectedNode(node);
  }
}

void QmitkFiberfoxView::ShowAdvancedOptions(int state)
{
  if (state)
  {
    m_Controls->m_AdvancedSignalOptionsFrame->setVisible(true);
    m_Controls->m_AdvancedOptionsBox_2->setChecked(true);
  }
  else
  {
    m_Controls->m_AdvancedSignalOptionsFrame->setVisible(false);
    m_Controls->m_AdvancedOptionsBox_2->setChecked(false);
  }
}

void QmitkFiberfoxView::Comp1ModelFrameVisibility(int index)
{
  ShowCompartmentWidget(1, index);
}

void QmitkFiberfoxView::Comp2ModelFrameVisibility(int index)
{
  ShowCompartmentWidget(2, index);
}

void QmitkFiberfoxView::Comp3ModelFrameVisibility(int index)
{
  ShowCompartmentWidget(3, index);
}

void QmitkFiberfoxView::Comp4ModelFrameVisibility(int index)
{
  ShowCompartmentWidget(4, index);
}

void QmitkFiberfoxView::OnAddMotion(int value)
{
  if (value>0)
    m_Controls->m_MotionArtifactFrame->setVisible(true);
  else
    m_Controls->m_MotionArtifactFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddDrift(int value)
{
  if (value>0)
    m_Controls->m_DriftFrame->setVisible(true);
  else
    m_Controls->m_DriftFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddAliasing(int value)
{
  if (value>0)
    m_Controls->m_AliasingFrame->setVisible(true);
  else
    m_Controls->m_AliasingFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddSpikes(int value)
{
  if (value>0)
    m_Controls->m_SpikeFrame->setVisible(true);
  else
    m_Controls->m_SpikeFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddEddy(int value)
{
  if (value>0)
    m_Controls->m_EddyFrame->setVisible(true);
  else
    m_Controls->m_EddyFrame->setVisible(false);
  EnableKspacePhase();
}

void QmitkFiberfoxView::OnAddDistortions(int value)
{
  if (value>0)
    m_Controls->m_DistortionsFrame->setVisible(true);
  else
    m_Controls->m_DistortionsFrame->setVisible(false);
  EnableKspacePhase();
}

void QmitkFiberfoxView::EnableKspacePhase()
{
  m_Controls->m_KspacePhaseFrame->setEnabled(m_Controls->m_AddDistortions->isChecked() || m_Controls->m_AddEddy->isChecked());
}

void QmitkFiberfoxView::ShowKspaceTimeSegments(int segments)
{
  // any negative count means automatic, as in itk::KspaceImageFilter
  const int index = segments==0 ? 0 : (segments<0 ? 1 : 2);
  if (index==2)
    m_Controls->m_KspaceTimeSegmentsBox->setValue(segments);
  m_Controls->m_KspacePhaseBox->setCurrentIndex(index);
  // also without a change of the index, which emits no signal
  m_Controls->m_KspaceTimeSegmentsLabel->setVisible(index==2);
  m_Controls->m_KspaceTimeSegmentsBox->setVisible(index==2);
}

void QmitkFiberfoxView::OnAddGhosts(int value)
{
  if (value>0)
    m_Controls->m_GhostFrame->setVisible(true);
  else
    m_Controls->m_GhostFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddNoise(int value)
{
  if (value>0)
    m_Controls->m_NoiseFrame->setVisible(true);
  else
    m_Controls->m_NoiseFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddRinging(int value)
{
  if (value>0)
    m_Controls->m_ZeroRinging->setVisible(true);
  else
    m_Controls->m_ZeroRinging->setVisible(false);
}

void QmitkFiberfoxView::GenerateImage()
{
  // set before the thread starts and cleared in AfterThread: a second click before the queued
  // BeforeThread hides the button must not replace the running filter
  if (m_ThreadIsRunning)
    return;

  // bvals/bvecs files and other inputs are read here; a bad one stops with its message
  try
  {
    StartSimulation();
  }
  catch (const itk::ExceptionObject& e)
  {
    QMessageBox::warning(nullptr, "No image generated", e.GetDescription());
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "No image generated", e.what());
  }
}

void QmitkFiberfoxView::StartSimulation()
{
  if (m_Controls->m_FiberBundleComboBox->GetSelectedNode().IsNull() && !mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( m_Controls->m_TemplateComboBox->GetSelectedNode()))
  {
    mitk::Image::Pointer image = mitk::ImageGenerator::GenerateGradientImage<unsigned int>(
          m_Controls->m_SizeX->value(),
          m_Controls->m_SizeY->value(),
          m_Controls->m_SizeZ->value(),
          m_Controls->m_SpacingX->value(),
          m_Controls->m_SpacingY->value(),
          m_Controls->m_SpacingZ->value());
    mitk::Point3D origin;
    origin[0] = m_Controls->m_SpacingX->value()/2;
    origin[1] = m_Controls->m_SpacingY->value()/2;
    origin[2] = m_Controls->m_SpacingZ->value()/2;
    image->SetOrigin(origin);

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData( image );
    node->SetName("Dummy");
    GetDataStorage()->Add(node);
    m_SelectedImageNode = node;

    mitk::BaseData::Pointer basedata = node->GetData();
    if (basedata.IsNotNull())
    {
      mitk::RenderingManager::GetInstance()->InitializeViews( basedata->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true );
      mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    }
    UpdateGui();
    QMessageBox::information(nullptr, "Template image generated", "You have selected no fiber bundle or diffusion-weighted image, which can be used to simulate a new diffusion-weighted image. A template image with the specified geometry has been generated that can be used to draw artificial fibers (see view 'Fiber Generator').");
  }
  else if (m_Controls->m_FiberBundleComboBox->GetSelectedNode().IsNotNull())
    SimulateImageFromFibers(m_Controls->m_FiberBundleComboBox->GetSelectedNode());
  else if ( mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( m_Controls->m_TemplateComboBox->GetSelectedNode()) )
    SimulateForExistingDwi(m_Controls->m_TemplateComboBox->GetSelectedNode());
  else
    QMessageBox::information(nullptr, "No image generated", "You have selected no fiber bundle or diffusion-weighted image, which can be used to simulate a new diffusion-weighted image.");
}

void QmitkFiberfoxView::SetFocus()
{
}

void QmitkFiberfoxView::SimulateForExistingDwi(mitk::DataNode* imageNode)
{
  bool isDiffusionImage( mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( dynamic_cast<mitk::Image *>(imageNode->GetData())) );
  if ( !isDiffusionImage )
  {
    return;
  }

  // no TI check: a DWI input is simulated without relaxation, so TI is not used
  UpdateParametersFromGui();

  mitk::Image::Pointer diffImg = dynamic_cast<mitk::Image*>(imageNode->GetData());
  ItkDwiType::Pointer itkVectorImagePointer = ItkDwiType::New();
  mitk::CastToItkImage(diffImg, itkVectorImagePointer);

  m_TractsToDwiFilter = itk::TractsToDWIImageFilter< short >::New();
  m_Parameters.m_Misc.m_ParentNode = imageNode;
  m_Parameters.m_SignalGen.m_SignalScale = 1;

  m_Parameters.m_Misc.m_ResultNode->SetName("b"+QString::number(m_Parameters.m_SignalGen.GetBvalue()).toStdString()
                                            +"_"+m_Parameters.m_Misc.m_SignalModelString
                                            +m_Parameters.m_Misc.m_ArtifactModelString);

  m_SimulatedParameters = SimulatedParameters(m_Parameters);
  m_Parameters.ApplyDirectionMatrix();
  m_TractsToDwiFilter->SetParameters(m_Parameters);
  m_TractsToDwiFilter->SetInputImage(itkVectorImagePointer);
  if (!ConfirmMemory())
  {
    m_TractsToDwiFilter = nullptr;
    return;
  }
  m_ThreadIsRunning = true;
  m_Thread.start(QThread::LowestPriority);
}

void QmitkFiberfoxView::SimulateImageFromFibers(mitk::DataNode* fiberNode)
{
  mitk::FiberBundle::Pointer fiberBundle = dynamic_cast<mitk::FiberBundle*>(fiberNode->GetData());
  if (fiberBundle->GetNumFibers()<=0) { return; }

  UpdateParametersFromGui();
  // TI >= TR with relaxation stops here with its message (GenerateImage), before the simulation thread starts
  m_Parameters.m_SignalGen.CheckInversionTime();

  m_TractsToDwiFilter = itk::TractsToDWIImageFilter< short >::New();
  m_Parameters.m_Misc.m_ParentNode = fiberNode;

  m_Parameters.m_Misc.m_ResultNode->SetName("b"+QString::number(m_Parameters.m_SignalGen.GetBvalue()).toStdString()
                                            +"_"+m_Parameters.m_Misc.m_SignalModelString
                                            +m_Parameters.m_Misc.m_ArtifactModelString);

  if ( m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull()
       && mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( dynamic_cast<mitk::Image*>
                                                                   (m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData()) ) )
  {
    mitk::Image::Pointer diffImg = dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData());
    const bool ok = SamplePrototypeKernels(diffImg, GetPrototypeModelsWithoutKernels());

    if (!ok)
    {
      QMessageBox::information( nullptr, "Simulation cancelled", "No valid prototype signals could be sampled.");
      return;
    }
  }
  else if ( !GetPrototypeModelsWithoutKernels().empty() )
  {
    QMessageBox::information( nullptr, "Simulation cancelled",
                              "Prototype signal but no diffusion-weighted image selected to sample signal from.");
    return;
  }

  m_SimulatedParameters = SimulatedParameters(m_Parameters);
  m_Parameters.ApplyDirectionMatrix();
  m_TractsToDwiFilter->SetParameters(m_Parameters);
  m_TractsToDwiFilter->SetFiberBundle(fiberBundle);
  if (!ConfirmMemory())
  {
    m_TractsToDwiFilter = nullptr;
    return;
  }
  m_ThreadIsRunning = true;
  m_Thread.start(QThread::LowestPriority);
}

bool QmitkFiberfoxView::ConfirmMemory()
{
  // a simulation beyond the physical memory swaps until the whole computer stops responding
  const double gigabyte = 1024.0*1024.0*1024.0;
  const double estimate = m_TractsToDwiFilter->EstimateMemoryBytes()/gigabyte;
  const double available = mitk::MemoryUtilities::GetAvailableSizeOfPhysicalRam()/gigabyte;
  if (available <= 0 || estimate <= available)   // 0: unknown
    return true;

  QMessageBox box(QMessageBox::Warning, "Not enough memory",
                  QString("The simulation needs about %1 GB of memory, but only %2 GB are available. "
                          "It can make the computer stop responding.\n\n"
                          "A smaller image, fewer volumes, no additional output images or ringing by zeroing "
                          "k-space frequencies need less memory.").arg(estimate, 0, 'f', 1).arg(available, 0, 'f', 1),
                  QMessageBox::Yes | QMessageBox::Cancel);
  box.button(QMessageBox::Yes)->setText("Start anyway");
  box.setDefaultButton(QMessageBox::Cancel);
  return box.exec() == QMessageBox::Yes;
}

void QmitkFiberfoxView::SetBvalsEdit()
{
  // SELECT FOLDER DIALOG
  std::string filename;
  filename = QFileDialog::getOpenFileName(nullptr, "Select bvals file", QString(filename.c_str())).toStdString();

  if (filename.empty())
    m_Controls->m_LoadBvalsEdit->setText("-");
  else
    m_Controls->m_LoadBvalsEdit->setText(QString(filename.c_str()));
}

void QmitkFiberfoxView::SetBvecsEdit()
{
  // SELECT FOLDER DIALOG
  std::string filename;
  filename = QFileDialog::getOpenFileName(nullptr, "Select bvecs file", QString(filename.c_str())).toStdString();

  if (filename.empty())
    m_Controls->m_LoadBvecsEdit->setText("-");
  else
    m_Controls->m_LoadBvecsEdit->setText(QString(filename.c_str()));
}

void QmitkFiberfoxView::SetOutputPath()
{
  // SELECT FOLDER DIALOG
  std::string outputPath;
  outputPath = QFileDialog::getExistingDirectory(nullptr, "Save images to...", QString(outputPath.c_str())).toStdString();

  if (outputPath.empty())
    m_Controls->m_SavePathEdit->setText("-");
  else
  {
    outputPath += "/";
    m_Controls->m_SavePathEdit->setText(QString(outputPath.c_str()));
  }
}

void QmitkFiberfoxView::UpdateGui()
{
  OnTlineChanged();
  m_Controls->m_GeometryFrame->setEnabled(true);
  m_Controls->m_GeometryMessage->setVisible(false);
  m_Controls->m_DiffusionPropsMessage->setVisible(false);

  m_Controls->m_LoadGradientsFrame->setVisible(false);
  m_Controls->m_GenerateGradientsFrame->setVisible(false);
  if (m_Controls->m_UseBvalsBvecsBox->isChecked())
    m_Controls->m_LoadGradientsFrame->setVisible(true);
  else
    m_Controls->m_GenerateGradientsFrame->setVisible(true);

  // Signal generation gui
  if (m_Controls->m_MaskComboBox->GetSelectedNode().IsNotNull() || m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull())
  {
    m_Controls->m_GeometryMessage->setVisible(true);
    m_Controls->m_GeometryFrame->setEnabled(false);
  }

  if ( m_Controls->m_TemplateComboBox->GetSelectedNode().IsNotNull()
       && mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(
         dynamic_cast<mitk::Image*>(m_Controls->m_TemplateComboBox->GetSelectedNode()->GetData()) ) )
  {
    m_Controls->m_DiffusionPropsMessage->setVisible(true);
    m_Controls->m_GeometryMessage->setVisible(true);
    m_Controls->m_GeometryFrame->setEnabled(false);

    m_Controls->m_LoadGradientsFrame->setVisible(false);
    m_Controls->m_GenerateGradientsFrame->setVisible(false);
  }
}
