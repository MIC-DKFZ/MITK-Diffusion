/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkIVIMView.h"

// qt
#include "qmessagebox.h"
#include "qclipboard.h"

// qwt
#include <qwt_plot_renderer.h>
#include <qwt_scale_engine.h>

// mitk
#include "mitkImage.h"
#include "mitkImageCast.h"

#include "mitkLookupTable.h"
#include "mitkLookupTableProperty.h"
#include <mitkSliceNavigationController.h>
#include <QmitkRenderWindow.h>

// itk
#include "itkScalarImageToHistogramGenerator.h"
#include <itkVectorIndexSelectionCastImageFilter.h>

// itk/mitk
#include <mitkDiffusionMasks.h>
#include <mitkDiffusionSignalSampling.h>
#include <mitkImageStatisticsHolder.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkLevelWindowProperty.h>

const std::string QmitkIVIMView::VIEW_ID = "org.mitk.views.ivim";

QmitkIVIMView::QmitkIVIMView() = default;

QmitkIVIMView::~QmitkIVIMView() = default;

void QmitkIVIMView::CreateQtPartControl( QWidget *parent )
{
  // hold update untill all elements are set
  this->m_HoldUpdate = true;

  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = std::make_unique<Ui::QmitkIVIMViewControls>();
    m_Controls->setupUi( parent );

    connect( m_Controls->m_ButtonStart, &QAbstractButton::clicked, this, &QmitkIVIMView::FittIVIMStart );
    connect( m_Controls->m_ButtonAutoThres, &QAbstractButton::clicked, this, &QmitkIVIMView::AutoThreshold );

    connect( m_Controls->m_MethodCombo, &QComboBox::currentIndexChanged, this, &QmitkIVIMView::OnIvimFitChanged );

    connect( m_Controls->m_DStarSlider, &QAbstractSlider::valueChanged, this, &QmitkIVIMView::DStarSlider );
    connect( m_Controls->m_BThreshSlider, &QAbstractSlider::valueChanged, this, &QmitkIVIMView::BThreshSlider );
    connect( m_Controls->m_S0ThreshSlider, &QAbstractSlider::valueChanged, this, &QmitkIVIMView::S0ThreshSlider );

    connect( m_Controls->m_CheckDStar, &QAbstractButton::clicked, this, &QmitkIVIMView::OnSliceChanged );
    connect( m_Controls->m_CheckD, &QAbstractButton::clicked, this, &QmitkIVIMView::OnSliceChanged );
    connect( m_Controls->m_Checkf, &QAbstractButton::clicked, this, &QmitkIVIMView::OnSliceChanged );

    connect( m_Controls->m_CurveClipboard, &QAbstractButton::clicked, this, &QmitkIVIMView::ClipboardCurveButtonClicked );
    connect( m_Controls->m_ValuesClipboard, &QAbstractButton::clicked, this, &QmitkIVIMView::ClipboardStatisticsButtonClicked );
    connect( m_Controls->m_SavePlot, &QAbstractButton::clicked, this, &QmitkIVIMView::SavePlotButtonClicked );

    // all kurtosis parameters refit the preview
    connect( m_Controls->m_KurtosisRangeWidget, &ctkRangeWidget::rangeChanged, this, &QmitkIVIMView::OnSliceChanged );
    connect( m_Controls->m_OmitBZeroCB, &QAbstractButton::toggled, this, &QmitkIVIMView::OnSliceChanged );
    connect( m_Controls->m_KurtosisFitScale, &QComboBox::currentIndexChanged, this, &QmitkIVIMView::OnSliceChanged );
    connect( m_Controls->m_UseKurtosisBoundsCB, &QAbstractButton::clicked, this, &QmitkIVIMView::OnSliceChanged );

    m_Controls->m_DwiBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_DwiBox->SetPredicate( mitk::NodePredicateIsDWI::New() );
    connect( m_Controls->m_DwiBox, &QComboBox::currentIndexChanged, this, &QmitkIVIMView::UpdateGui );

    m_Controls->m_MaskBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_MaskBox->SetZeroEntryText("--");
    m_Controls->m_MaskBox->SetPredicate( mitk::GetBinaryMaskPredicate() );
    connect( m_Controls->m_MaskBox, &QComboBox::currentIndexChanged, this, &QmitkIVIMView::UpdateGui );

    connect( m_Controls->m_ModelTabSelectionWidget, &QTabWidget::currentChanged, this, &QmitkIVIMView::OnModelTabChanged );
  }

  m_Controls->m_DStarLabel->setText(QString::number(this->GetDStar()));
  m_Controls->m_BThreshLabel->setText(QString::number(this->GetBThreshold()));
  m_Controls->m_S0ThreshLabel->setText(QString::number(this->GetSignalThreshold()));

  m_Controls->m_Warning->setVisible(false);

  OnIvimFitChanged(m_Controls->m_MethodCombo->currentIndex());

  m_Controls->m_KurtosisRangeWidget->setSingleStep(0.1);
  m_Controls->m_KurtosisRangeWidget->setRange( 0.0, 10.0 );
  m_Controls->m_KurtosisRangeWidget->setMaximumValue( 5.0 );

  // release update block after the UI-elements were all set
  this->m_HoldUpdate = false;

  InitChart("S/S0");

  this->ActivateSliceListener();
}

void QmitkIVIMView::ActivateSliceListener()
{
  if (m_ListenerActive || !this->GetRenderWindowPart())
    return;

  m_SliceChangeListener.RenderWindowPartActivated(this->GetRenderWindowPart());
  connect(&m_SliceChangeListener, &QmitkSliceNavigationListener::SliceChanged, this, &QmitkIVIMView::OnSliceChanged);
  m_ListenerActive = true;
}

double QmitkIVIMView::GetDStar() const
{
  return m_Controls->m_DStarSlider->value() / 1000.0;
}

double QmitkIVIMView::GetBThreshold() const
{
  return m_Controls->m_BThreshSlider->value() * 5.0;
}

double QmitkIVIMView::GetSignalThreshold() const
{
  return m_Controls->m_S0ThreshSlider->value() * 0.5;
}

void QmitkIVIMView::OnModelTabChanged(int tab)
{
  if (tab==0)
    InitChart("S/S0");
  else if (tab==1)
    InitChart("S");

  UpdateGui();
}

void QmitkIVIMView::AddSecondFitPlot()
{
  m_SecondFitSignalCurve->setVisible(true);
  m_SecondFitSignalCurve->setItemAttribute(QwtPlotItem::Legend, true);
}

// QmitkPlotWidget cannot remove a single curve, so the second-fit curve is hidden instead.
void QmitkIVIMView::RemoveSecondFitPlot()
{
  m_SecondFitSignalCurve->setVisible(false);
  m_SecondFitSignalCurve->setItemAttribute(QwtPlotItem::Legend, false);
}

void QmitkIVIMView::InitChart(const char* yAxisTitle)
{
  QmitkPlotWidget* plot = m_Controls->m_ChartWidget;
  plot->Clear();

  // The canvas is white, so dark colors replace the white of the former dark chart theme.
  m_DPartCurveId = plot->InsertCurve("D-part of fitted model");
  plot->SetCurvePen(m_DPartCurveId, QPen(Qt::blue, 2, Qt::DashLine));

  m_FitCurveId = plot->InsertCurve("fitted model");
  plot->SetCurvePen(m_FitCurveId, QPen(Qt::red, 2));

  m_SignalCurveId = plot->InsertCurve("signal values");
  plot->SetCurvePen(m_SignalCurveId, QPen(Qt::NoPen));
  plot->SetCurveSymbol(m_SignalCurveId, new QwtSymbol(QwtSymbol::Diamond, QBrush(Qt::black), QPen(Qt::black), QSize(8, 8)));

  m_SecondFitSignalCurveId = plot->InsertCurve("signal values (used for second fit)");
  plot->SetCurvePen(m_SecondFitSignalCurveId, QPen(Qt::NoPen));
  plot->SetCurveSymbol(m_SecondFitSignalCurveId, new QwtSymbol(QwtSymbol::XCross, QBrush(), QPen(Qt::darkGray), QSize(8, 8)));
  m_SecondFitSignalCurve = static_cast<QwtPlotCurve*>(plot->GetPlot()->itemList(QwtPlotItem::Rtti_PlotCurve).back());
  RemoveSecondFitPlot();

  plot->SetAxisTitle(QwtPlot::xBottom, "b-value");
  plot->SetAxisTitle(QwtPlot::yLeft, yAxisTitle);
  plot->GetPlot()->setAxisScaleEngine(QwtPlot::yLeft, new QwtLogScaleEngine);
  plot->SetLegend(new QwtLegend, QwtPlot::BottomLegend);
  plot->Replot();
}

void QmitkIVIMView::SetFocus()
{
  m_Controls->m_ButtonAutoThres->setFocus();
}

void QmitkIVIMView::OnIvimFitChanged(int val)
{
  switch(val)
  {
  case 0:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(false);
    break;
  case 1:
    m_Controls->m_DstarFrame->setVisible(true);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(false);
    break;
  case 2:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(true);
    break;
  case 3:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(true);
    break;
  }

  OnSliceChanged();
}

void QmitkIVIMView::DStarSlider()
{
  m_Controls->m_DStarLabel->setText(QString::number(this->GetDStar()));
  OnSliceChanged();
}

void QmitkIVIMView::BThreshSlider()
{
  m_Controls->m_BThreshLabel->setText(QString::number(this->GetBThreshold()));
  OnSliceChanged();
}

void QmitkIVIMView::S0ThreshSlider()
{
  m_Controls->m_S0ThreshLabel->setText(QString::number(this->GetSignalThreshold()));
  OnSliceChanged();
}

void QmitkIVIMView::UpdateGui()
{
  m_Controls->m_FittedParamsLabel->setText("");
  if (m_Controls->m_DwiBox->GetSelectedNode().IsNotNull())
  {
    m_Controls->m_ChartWidget->setVisible(true);
    m_HoldUpdate = false;
  }
  else
  {
    m_Controls->m_ChartWidget->setVisible(false);
  }

  m_Controls->m_ButtonStart->setEnabled( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() );
  m_Controls->m_ButtonAutoThres->setEnabled( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() );

  m_Controls->m_ControlsFrame->setEnabled( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() );
  m_Controls->m_BottomControlsFrame->setEnabled( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() );

  OnSliceChanged();
}

void QmitkIVIMView::AutoThreshold()
{
  if (m_Controls->m_DwiBox->GetSelectedNode().IsNull())
  {
    // Nothing selected. Inform the user and return
    QMessageBox::information( nullptr, "IVIM", "Please load and select a diffusion image before starting image processing.");
    return;
  }

  mitk::Image* dimg = dynamic_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());

  if (!dimg)
  {
    // Nothing selected. Inform the user and return
    QMessageBox::information( nullptr, "IVIM", "No valid diffusion image was found.");
    return;
  }

  // find bzero index
  int index = -1;
  auto directions = mitk::DiffusionPropertyHelper::GetGradientContainer(dimg);
  for(DirContainerType::ConstIterator it = directions->Begin();
      it != directions->End(); ++it)
  {
    index++;
    GradientDirectionType g = it.Value();
    if(g[0] == 0 && g[1] == 0 && g[2] == 0 )
      break;
  }

  VecImgType::Pointer vecimg = VecImgType::New();
  mitk::CastToItkImage(dimg, vecimg);

  int vecLength = vecimg->GetVectorLength();
  index = index > vecLength-1 ? vecLength-1 : index;

  MITK_INFO << "Performing Histogram Analysis on Channel" << index;

  // A DWI cannot be cast to a scalar image, so the b=0 channel is extracted.
  using ImgType = itk::Image<short, 3>;
  using ChannelFilterType = itk::VectorIndexSelectionCastImageFilter<VecImgType, ImgType>;
  ChannelFilterType::Pointer channelFilter = ChannelFilterType::New();
  channelFilter->SetInput(vecimg);
  channelFilter->SetIndex(index);
  channelFilter->Update();
  ImgType::Pointer img = channelFilter->GetOutput();

  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImgType>;
  using HistogramType = HistogramGeneratorType::HistogramType;

  HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput( img );
  histogramGenerator->SetMarginalScale( 10 ); // Defines y-margin width of histogram
  histogramGenerator->SetNumberOfBins( 100 ); // CT range [-1024, +2048] --> bin size 4 values
  histogramGenerator->SetHistogramMin(  dimg->GetStatistics()->GetScalarValueMin() );
  histogramGenerator->SetHistogramMax(  dimg->GetStatistics()->GetScalarValueMax() * .5 );
  histogramGenerator->Compute();

  HistogramType::ConstIterator iter = histogramGenerator->GetOutput()->Begin();
  float maxFreq = 0;
  float maxValue = 0;
  while ( iter != histogramGenerator->GetOutput()->End() )
  {
    if(iter.GetFrequency() > maxFreq)
    {
      maxFreq = iter.GetFrequency();
      maxValue = iter.GetMeasurementVector()[0];
    }
    ++iter;
  }

  maxValue *= 2;

  // The slider position is twice the threshold; the slider grows when the threshold is beyond its
  // range, so that the slider always holds the threshold the fit uses.
  const int sliderPos = maxValue * 2;
  if (sliderPos > m_Controls->m_S0ThreshSlider->maximum())
    m_Controls->m_S0ThreshSlider->setMaximum(sliderPos);
  m_Controls->m_S0ThreshSlider->setValue(sliderPos);
}

void QmitkIVIMView::FittIVIMStart()
{
  auto dwiNode = m_Controls->m_DwiBox->GetSelectedNode();
  if (dwiNode.IsNull())
  {
    QMessageBox::information( nullptr, "IVIM", "No valid diffusion-weighted image selected.");
    return;
  }
  const auto* dwi = dynamic_cast<mitk::Image*>(dwiNode->GetData());

  try
  {
    mitk::Image::ConstPointer mask;
    if (m_Controls->m_MaskBox->GetSelectedNode().IsNotNull())
      mask = mitk::GetMaskImage(m_Controls->m_MaskBox->GetSelectedNode()->GetData());

    if( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
    {
      const auto maps = mitk::kurtosis::ComputeMaps(dwi, mask, this->GetKurtosisFitConfiguration(), m_Controls->m_SigmaSpinBox->value());
      const QString method = "_Method-" + m_Controls->m_KurtosisFitScale->currentText();
      if (m_Controls->m_CheckKurtD->isChecked())
        this->AddMapNode(maps.d, "Kurtosis_DMap" + method, std::nullopt, dwiNode);
      if (m_Controls->m_CheckKurtK->isChecked())
        this->AddMapNode(maps.k, "Kurtosis_KMap" + method, std::nullopt, dwiNode);
    }
    else
    {
      const auto maps = mitk::ivim::ComputeMaps(dwi, mask, this->GetIvimFitParameters());
      const QString method = QString::number(m_Controls->m_MethodCombo->currentIndex() + 1);
      if (m_Controls->m_CheckDStar->isChecked())
        this->AddMapNode(maps.dStar, "D* Fit-" + method, 0.2, dwiNode);
      if (m_Controls->m_CheckD->isChecked())
        this->AddMapNode(maps.d, "D Fit-" + method, 0.003, dwiNode);
      if (m_Controls->m_Checkf->isChecked())
        this->AddMapNode(maps.f, "f Fit-" + method, 1, dwiNode);
    }
  }
  catch (const itk::ExceptionObject& e)
  {
    QMessageBox::warning( nullptr, "IVIM", QString("Fit failed: ") + e.GetDescription() );
    return;
  }

  if (this->GetRenderWindowPart())
    this->GetRenderWindowPart()->RequestUpdate();
}

void QmitkIVIMView::AddMapNode(mitk::Image* map, const QString& name, std::optional<double> levelWindowMax, mitk::DataNode* parent)
{
  auto lut = mitk::LookupTable::New();
  lut->SetType(mitk::LookupTable::JET);

  auto node = mitk::DataNode::New();
  node->SetData(map);
  node->SetName(name.toStdString());
  node->SetProperty("LookupTable", mitk::LookupTableProperty::New(lut));
  if (levelWindowMax)
  {
    mitk::LevelWindow levelWindow;
    levelWindow.SetRangeMinMax(0, *levelWindowMax);
    node->SetProperty("levelwindow", mitk::LevelWindowProperty::New(levelWindow));
  }
  this->GetDataStorage()->Add(node, parent);
}

void QmitkIVIMView::OnSliceChanged()
{
  if(m_HoldUpdate || !m_Visible)
    return;

  m_Controls->m_Warning->setVisible(false);
  if(m_Controls->m_DwiBox->GetSelectedNode().IsNull())
    return;

  mitk::Image::Pointer diffusionImg = dynamic_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());
  mitk::Image::ConstPointer maskImg;
  if (m_Controls->m_MaskBox->GetSelectedNode().IsNotNull())
    maskImg = mitk::GetMaskImage(m_Controls->m_MaskBox->GetSelectedNode()->GetData());

  if (!this->GetRenderWindowPart())
    return;

  this->ActivateSliceListener();

  vnl_vector<double> signal;
  if(maskImg.IsNull())
  {
    const auto voxelSignal = mitk::GetDiffusionSignal(diffusionImg, this->GetRenderWindowPart()->GetSelectedPosition());
    if (!voxelSignal)
    {
      m_Controls->m_Warning->setText(QString("Crosshair position not inside of selected diffusion weighted image. Reinit needed!"));
      m_Controls->m_Warning->setVisible(true);
      return;
    }
    signal = *voxelSignal;
  }
  else
  {
    // The mean is fitted in double, not rounded to the pixel type of the DWI.
    const auto roi = mitk::AverageDiffusionSignalInMask(diffusionImg, maskImg);
    if (roi.numberOfVoxels == 0)
    {
      m_Controls->m_Warning->setText(QString("The mask does not overlap the diffusion-weighted image."));
      m_Controls->m_Warning->setVisible(true);
      return;
    }

    m_Controls->m_Warning->setText(QString("Averaging ")+QString::number(roi.numberOfVoxels)+QString(" voxels!"));
    m_Controls->m_Warning->setVisible(true);
    signal = roi.mean;
  }

  const auto dirs = mitk::DiffusionPropertyHelper::GetGradientContainer(diffusionImg);
  const auto bval = mitk::DiffusionPropertyHelper::GetReferenceBValue(diffusionImg);
  if( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
    FitKurtosisVoxel(signal, dirs, bval);
  else
    FitIvimVoxel(signal, dirs, bval);

  // do not update until selection changed, the values will remain the same as long as the mask is selected!
  if (maskImg.IsNotNull())
    m_HoldUpdate = true;
}

mitk::kurtosis::FitConfiguration QmitkIVIMView::GetKurtosisFitConfiguration() const
{
  mitk::kurtosis::FitConfiguration configuration;
  configuration.omitBZero = m_Controls->m_OmitBZeroCB->isChecked();
  configuration.scale = static_cast<mitk::kurtosis::FitScale>(m_Controls->m_KurtosisFitScale->currentIndex());
  if (m_Controls->m_UseKurtosisBoundsCB->isChecked())
  {
    configuration.useKBounds = true;
    configuration.kLowerBound = m_Controls->m_KurtosisRangeWidget->minimumValue();
    configuration.kUpperBound = m_Controls->m_KurtosisRangeWidget->maximumValue();
  }
  return configuration;
}

void QmitkIVIMView::FitKurtosisVoxel(const vnl_vector<double>& signal, DirContainerType::ConstPointer dirs, float bval)
{
  m_KurtosisMeasurements = signal;

  try
  {
    m_KurtosisBValues = mitk::kurtosis::GetBValues(dirs, bval);
    m_KurtosisFit = mitk::kurtosis::FitVoxel(m_KurtosisMeasurements, m_KurtosisBValues, this->GetKurtosisFitConfiguration());
  }
  catch (const itk::ExceptionObject& ex)
  {
    m_Controls->m_Warning->setText(QString("Kurtosis fit not possible: ") + ex.GetDescription());
    m_Controls->m_Warning->setVisible(true);
    return;
  }

  QString param_label_text("K=%2, D=%1");
  param_label_text = param_label_text.arg( m_KurtosisFit.k, 4);
  param_label_text = param_label_text.arg( m_KurtosisFit.d, 4);
  m_Controls->m_FittedParamsLabel->setText(param_label_text);

  const double maxb = m_KurtosisBValues.max_value();
  const double S0 = m_KurtosisFit.bZeroFitted ? m_KurtosisFit.fittedBZero : m_KurtosisFit.meanBZero;

  std::vector< std::pair<double, double> > d_line;
  d_line.emplace_back(0, S0);
  d_line.emplace_back(maxb, S0*exp(-maxb * m_KurtosisFit.d));
  m_Controls->m_ChartWidget->SetCurveData(m_DPartCurveId, d_line);

  const unsigned int num_samples = 50;
  std::vector< std::pair<double, double> > y;

  for( unsigned int i=0; i<=num_samples; ++i)
  {
    double b = (((1.0)*i)/(1.0*num_samples))*maxb;
    y.emplace_back(b, S0 * exp( -b * m_KurtosisFit.d + b*b * m_KurtosisFit.d * m_KurtosisFit.d * m_KurtosisFit.k / 6.0 ));
  }
  m_Controls->m_ChartWidget->SetCurveData(m_FitCurveId, y);

  std::vector< std::pair<double, double> > y_meas;
  for (unsigned int i=0; i<m_KurtosisMeasurements.size(); ++i)
  {
    if (!m_Controls->m_OmitBZeroCB->isChecked() || m_KurtosisBValues[i] > 0.01)
      y_meas.emplace_back(m_KurtosisBValues[i], m_KurtosisMeasurements[i]);
  }

  m_Controls->m_ChartWidget->SetCurveData(m_SignalCurveId, y_meas);
  m_Controls->m_ChartWidget->Replot();
}

mitk::ivim::FitParameters QmitkIVIMView::GetIvimFitParameters() const
{
  mitk::ivim::FitParameters parameters;
  parameters.method = static_cast<mitk::ivim::Method>(m_Controls->m_MethodCombo->currentIndex());
  parameters.signalThreshold = this->GetSignalThreshold();
  parameters.dStar = this->GetDStar();
  parameters.bThreshold = this->GetBThreshold();
  parameters.fitDStar = m_Controls->m_CheckDStar->isChecked();
  return parameters;
}

void QmitkIVIMView::FitIvimVoxel(const vnl_vector<double>& signal, DirContainerType::ConstPointer dirs, float bval)
{
  // The preview always shows D*, also when no D* map is generated.
  auto parameters = this->GetIvimFitParameters();
  parameters.fitDStar = true;

  try
  {
    m_IvimProtocol = mitk::ivim::MakeProtocol(dirs, bval);
    m_IvimFit = mitk::ivim::FitVoxel(signal, m_IvimProtocol, parameters);
  }
  catch (const itk::ExceptionObject& ex)
  {
    MITK_INFO << ex;
    m_Controls->m_Warning->setText(QString("IVIM fit not possible: ") + ex.GetDescription());
    m_Controls->m_Warning->setVisible(true);
    return;
  }

  QString param_label_text("f=%1, D=%2, D*=%3");
  param_label_text = param_label_text.arg(m_IvimFit.f, 4);
  param_label_text = param_label_text.arg(m_IvimFit.d, 4);
  param_label_text = param_label_text.arg(m_IvimFit.dStar, 4);
  m_Controls->m_FittedParamsLabel->setText(param_label_text);

  const double maxb = m_IvimProtocol.bValues.max_value();

  std::vector< std::pair<double, double> > d_line;
  d_line.emplace_back(0, 1 - m_IvimFit.unclampedF);
  d_line.emplace_back(maxb, d_line[0].second * exp(-maxb * m_IvimFit.d));
  m_Controls->m_ChartWidget->SetCurveData(m_DPartCurveId, d_line);

  std::vector< std::pair<double, double> > y;
  const int nsampling = 50;
  const double f = 1 - m_IvimFit.unclampedF;
  for (int i = 0; i <= nsampling; i++)
  {
    const double x = (((1.0) * i) / (1.0 * nsampling)) * maxb;
    y.emplace_back(x, f * exp(-x * m_IvimFit.d) + (1 - f) * exp(-x * (m_IvimFit.d + m_IvimFit.dStar)));
  }
  m_Controls->m_ChartWidget->SetCurveData(m_FitCurveId, y);

  std::vector< std::pair<double, double> > y_meas;
  for (unsigned int i = 0; i < m_IvimFit.firstFitMeasurements.size(); ++i)
    y_meas.emplace_back(m_IvimFit.firstFitBValues[i], m_IvimFit.firstFitMeasurements[i]);
  m_Controls->m_ChartWidget->SetCurveData(m_SignalCurveId, y_meas);

  if (m_IvimFit.secondFitBValues.size() > 0)
  {
    AddSecondFitPlot();
    std::vector< std::pair<double, double> > additonal_meas;
    for (unsigned int i = 0; i < m_IvimFit.secondFitBValues.size(); ++i)
      additonal_meas.emplace_back(m_IvimFit.secondFitBValues[i], m_IvimFit.secondFitMeasurements[i]);
    m_Controls->m_ChartWidget->SetCurveData(m_SecondFitSignalCurveId, additonal_meas);
  }
  else
  {
    RemoveSecondFitPlot();
  }
  m_Controls->m_ChartWidget->Replot();
}

void QmitkIVIMView::ClipboardCurveButtonClicked()
{
  // Kurtosis
  if ( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
  {
    std::stringstream ss;
    QString clipboard("Measurement Points\n");

    ss << m_KurtosisBValues << "\n" << m_KurtosisMeasurements << "\n\n";
    ss << "Fitted Values ( D  K  [b_0] ) \n" << m_KurtosisFit.d << " " << m_KurtosisFit.k;

    if( m_KurtosisFit.bZeroFitted )
      ss << " " << m_KurtosisFit.fittedBZero;

    ss << "\n\n";
    clipboard.append( QString( ss.str().c_str() ));

    ss.str( std::string() );
    ss.clear();

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
  else
  {
    std::stringstream ss;
    QString clipboard("Normalized Measurement Points\n");
    ss << m_IvimProtocol.bValues << "\n" << m_IvimFit.allMeasurements << "\n\n";
    ss << "Fitted Values ( f  D  D* ) \n" << m_IvimFit.f << " " << m_IvimFit.d << " " << m_IvimFit.dStar;


    ss << "\n\n";
    clipboard.append( QString( ss.str().c_str() ));

    ss.str( std::string() );
    ss.clear();

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
}

void QmitkIVIMView::SavePlotButtonClicked()
{
  // Opens a file dialog to choose the file name and format.
  QwtPlotRenderer renderer;
  renderer.exportTo(m_Controls->m_ChartWidget->GetPlot(), "ivim_plot.png");
}


void QmitkIVIMView::ClipboardStatisticsButtonClicked()
{
  // Kurtosis
  if ( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
  {

    QString clipboard( "D \t K \n" );
    clipboard = clipboard.append( "%L1 \t %L2" )
        .arg( m_KurtosisFit.d, 0, 'f', 10 )
        .arg( m_KurtosisFit.k, 0, 'f', 10 ) ;

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
  else
  {
    QString clipboard( "f \t D \t D* \n" );
    clipboard = clipboard.append( "%L1 \t %L2 \t %L3" )
        .arg( m_IvimFit.f, 0, 'f', 10 )
        .arg( m_IvimFit.d, 0, 'f', 10 )
        .arg( m_IvimFit.dStar, 0, 'f', 10 ) ;

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
}

void QmitkIVIMView::Activated()
{
}

void QmitkIVIMView::Deactivated()
{
}

void QmitkIVIMView::Visible()
{
  m_Visible = true;
  this->ActivateSliceListener();
}

void QmitkIVIMView::Hidden()
{
  m_Visible = false;
}
