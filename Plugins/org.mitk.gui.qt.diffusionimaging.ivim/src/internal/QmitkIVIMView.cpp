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
#include "QmitkIVIMView.h"

// qt
#include "qmessagebox.h"
#include "qclipboard.h"

// mitk
#include "mitkImage.h"
#include "mitkImageCast.h"

#include "mitkLookupTable.h"
#include "mitkLookupTableProperty.h"
#include <mitkSliceNavigationController.h>
#include <QmitkRenderWindow.h>

// itk
#include "itkScalarImageToHistogramGenerator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// itk/mitk
#include "itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.h"
#include "itkRegularizedIVIMReconstructionFilter.h"
#include "mitkImageCast.h"
#include <mitkImageStatisticsHolder.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkLevelWindowProperty.h>

const std::string QmitkIVIMView::VIEW_ID = "org.mitk.views.ivim";

QmitkIVIMView::QmitkIVIMView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_Active(false)
  , m_Visible(false)
  , m_HoldUpdate(false)
{
}

QmitkIVIMView::~QmitkIVIMView()
{
}

void QmitkIVIMView::CreateQtPartControl( QWidget *parent )
{
  // hold update untill all elements are set
  this->m_HoldUpdate = true;

  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkIVIMViewControls;
    m_Controls->setupUi( parent );

    connect( m_Controls->m_ButtonStart, SIGNAL(clicked()), this, SLOT(FittIVIMStart()) );
    connect( m_Controls->m_ButtonAutoThres, SIGNAL(clicked()), this, SLOT(AutoThreshold()) );

    connect( m_Controls->m_MethodCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(OnIvimFitChanged(int)) );

    connect( m_Controls->m_DStarSlider, SIGNAL(valueChanged(int)), this, SLOT(DStarSlider(int)) );
    connect( m_Controls->m_BThreshSlider, SIGNAL(valueChanged(int)), this, SLOT(BThreshSlider(int)) );
    connect( m_Controls->m_S0ThreshSlider, SIGNAL(valueChanged(int)), this, SLOT(S0ThreshSlider(int)) );
    connect( m_Controls->m_NumItSlider, SIGNAL(valueChanged(int)), this, SLOT(NumItsSlider(int)) );
    connect( m_Controls->m_LambdaSlider, SIGNAL(valueChanged(int)), this, SLOT(LambdaSlider(int)) );

    connect( m_Controls->m_CheckDStar, SIGNAL(clicked()), this, SLOT(Checkbox()) );
    connect( m_Controls->m_CheckD, SIGNAL(clicked()), this, SLOT(Checkbox()) );
    connect( m_Controls->m_Checkf, SIGNAL(clicked()), this, SLOT(Checkbox()) );

    connect( m_Controls->m_CurveClipboard, SIGNAL(clicked()), this, SLOT(ClipboardCurveButtonClicked()) );
    connect( m_Controls->m_ValuesClipboard, SIGNAL(clicked()), this, SLOT(ClipboardStatisticsButtonClicked()) );
    connect( m_Controls->m_SavePlot, SIGNAL(clicked()), this, SLOT(SavePlotButtonClicked()) );

    // connect all kurtosis actions to a recompute
    connect( m_Controls->m_KurtosisRangeWidget, SIGNAL( rangeChanged(double, double)), this, SLOT(OnKurtosisParamsChanged() ) );
    connect( m_Controls->m_OmitBZeroCB, SIGNAL( stateChanged(int) ), this, SLOT( OnKurtosisParamsChanged() ) );
    connect( m_Controls->m_KurtosisFitScale, SIGNAL( currentIndexChanged(int)), this, SLOT( OnKurtosisParamsChanged() ) );
    connect( m_Controls->m_UseKurtosisBoundsCB, SIGNAL(clicked() ), this, SLOT( OnKurtosisParamsChanged() ) );

    m_Controls->m_DwiBox->SetDataStorage(this->GetDataStorage());
    mitk::NodePredicateIsDWI::Pointer isDwi = mitk::NodePredicateIsDWI::New();
    m_Controls->m_DwiBox->SetPredicate( isDwi );
    connect( (QObject*)(m_Controls->m_DwiBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));

    m_Controls->m_MaskBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_MaskBox->SetZeroEntryText("--");
    mitk::TNodePredicateDataType<mitk::Image>::Pointer isImagePredicate = mitk::TNodePredicateDataType<mitk::Image>::New();
    mitk::NodePredicateProperty::Pointer isBinaryPredicate = mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
    mitk::NodePredicateDimension::Pointer is3D = mitk::NodePredicateDimension::New(3);
    m_Controls->m_MaskBox->SetPredicate( mitk::NodePredicateAnd::New(isBinaryPredicate, mitk::NodePredicateAnd::New(isImagePredicate, is3D)) );
    connect( (QObject*)(m_Controls->m_MaskBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));

    connect( (QObject*)(m_Controls->m_ModelTabSelectionWidget), SIGNAL(currentChanged(int)), this, SLOT(OnModelTabChanged(int)));
  }

  QString dstar = QString::number(m_Controls->m_DStarSlider->value()/1000.0);
  m_Controls->m_DStarLabel->setText(dstar);

  QString bthresh = QString::number(m_Controls->m_BThreshSlider->value()*5.0);
  m_Controls->m_BThreshLabel->setText(bthresh);

  QString s0thresh = QString::number(m_Controls->m_S0ThreshSlider->value()*0.5);
  m_Controls->m_S0ThreshLabel->setText(s0thresh);

  QString numits = QString::number(m_Controls->m_NumItSlider->value());
  m_Controls->m_NumItsLabel->setText(numits);

  QString lambda = QString::number(m_Controls->m_LambdaSlider->value()*.00001);
  m_Controls->m_LambdaLabel->setText(lambda);

  m_Controls->m_Warning->setVisible(false);

  OnIvimFitChanged(m_Controls->m_MethodCombo->currentIndex());

  m_Controls->m_KurtosisRangeWidget->setSingleStep(0.1);
  m_Controls->m_KurtosisRangeWidget->setRange( 0.0, 10.0 );
  m_Controls->m_KurtosisRangeWidget->setMaximumValue( 5.0 );

  // LogScale not working yet, have to fix that first
  // m_Controls->m_KurtosisFitScale->setEnabled(false);


  //m_Controls->m_MaximalBValueWidget->setVisible( false );

  // release update block after the UI-elements were all set
  this->m_HoldUpdate = false;

  QmitkIVIMView::InitChartIvim();

  m_ListenerActive = false;

  if (this->GetRenderWindowPart())
  {
    m_SliceChangeListener.RenderWindowPartActivated(this->GetRenderWindowPart());
    connect(&m_SliceChangeListener, SIGNAL(SliceChanged()), this, SLOT(OnSliceChanged()));
    m_ListenerActive = true;
  }
}

void QmitkIVIMView::OnModelTabChanged(int tab)
{
  if (tab==0)
    InitChartIvim();
  else if (tab==1)
    InitChartKurtosis();

  UpdateGui();
}

void QmitkIVIMView::AddSecondFitPlot()
{
  if (m_Controls->m_ChartWidget->GetDataElementByLabel("signal values (used for second fit)") == nullptr)
  {
    std::vector< std::pair<double, double> > init_data;
    m_Controls->m_ChartWidget->AddData2D(init_data, "signal values (used for second fit)", QmitkChartWidget::ChartType::scatter);
    m_Controls->m_ChartWidget->SetColor("signal values (used for second fit)", "white");
    m_Controls->m_ChartWidget->SetMarkerSymbol("signal values (used for second fit)", QmitkChartWidget::MarkerSymbol::x_thin);

    m_Controls->m_ChartWidget->Show(true);
    m_Controls->m_ChartWidget->SetShowDataPoints(false);
    m_Controls->m_ChartWidget->SetShowSubchart(false);
  }
}

void QmitkIVIMView::RemoveSecondFitPlot()
{
  if (m_Controls->m_ChartWidget->GetDataElementByLabel("signal values (used for second fit)") != nullptr)
  {
    m_Controls->m_ChartWidget->RemoveData("signal values (used for second fit)");

    m_Controls->m_ChartWidget->Show(true);
    m_Controls->m_ChartWidget->SetShowDataPoints(false);
    m_Controls->m_ChartWidget->SetShowSubchart(false);
  }
}

void QmitkIVIMView::InitChartIvim()
{
  m_Controls->m_ChartWidget->Clear();
  std::vector< std::pair<double, double> > init_data;
  m_Controls->m_ChartWidget->AddData2D(init_data, "D-part of fitted model", QmitkChartWidget::ChartType::line);
  m_Controls->m_ChartWidget->AddData2D(init_data, "fitted model", QmitkChartWidget::ChartType::line);
  m_Controls->m_ChartWidget->AddData2D(init_data, "signal values", QmitkChartWidget::ChartType::scatter);

  m_Controls->m_ChartWidget->SetColor("fitted model", "red");
  m_Controls->m_ChartWidget->SetColor("signal values", "white");
  m_Controls->m_ChartWidget->SetColor("D-part of fitted model", "cyan");
  m_Controls->m_ChartWidget->SetYAxisScale(QmitkChartWidget::AxisScale::log);

  m_Controls->m_ChartWidget->SetYAxisLabel("S/S0");
  m_Controls->m_ChartWidget->SetXAxisLabel("b-value");

  m_Controls->m_ChartWidget->SetLineStyle("fitted model", QmitkChartWidget::LineStyle::solid);
  m_Controls->m_ChartWidget->SetLineStyle("D-part of fitted model", QmitkChartWidget::LineStyle::dashed);

  m_Controls->m_ChartWidget->SetMarkerSymbol("signal values", QmitkChartWidget::MarkerSymbol::diamond);

  m_Controls->m_ChartWidget->Show(true);
  m_Controls->m_ChartWidget->SetShowDataPoints(false);
  m_Controls->m_ChartWidget->SetShowSubchart(false);

}

void QmitkIVIMView::InitChartKurtosis()
{
  m_Controls->m_ChartWidget->Clear();
  std::vector< std::pair<double, double> > init_data;
  m_Controls->m_ChartWidget->AddData2D(init_data, "D-part of fitted model", QmitkChartWidget::ChartType::line);
  m_Controls->m_ChartWidget->AddData2D(init_data, "fitted model", QmitkChartWidget::ChartType::line);
  m_Controls->m_ChartWidget->AddData2D(init_data, "signal values", QmitkChartWidget::ChartType::scatter);

  m_Controls->m_ChartWidget->SetColor("fitted model", "red");
  m_Controls->m_ChartWidget->SetColor("signal values", "white");
  m_Controls->m_ChartWidget->SetColor("D-part of fitted model", "cyan");
  m_Controls->m_ChartWidget->SetYAxisScale(QmitkChartWidget::AxisScale::log);

  m_Controls->m_ChartWidget->SetYAxisLabel("S");
  m_Controls->m_ChartWidget->SetXAxisLabel("b-value");

  m_Controls->m_ChartWidget->SetLineStyle("fitted model", QmitkChartWidget::LineStyle::solid);
  m_Controls->m_ChartWidget->SetLineStyle("D-part of fitted model", QmitkChartWidget::LineStyle::dashed);

  m_Controls->m_ChartWidget->SetMarkerSymbol("signal values", QmitkChartWidget::MarkerSymbol::diamond);

  m_Controls->m_ChartWidget->Show(true);
  m_Controls->m_ChartWidget->SetShowDataPoints(false);
  m_Controls->m_ChartWidget->SetShowSubchart(false);
}

void QmitkIVIMView::SetFocus()
{
  m_Controls->m_ButtonAutoThres->setFocus();
}

void QmitkIVIMView::Checkbox()
{
  OnSliceChanged();
}

void QmitkIVIMView::OnIvimFitChanged(int val)
{
  switch(val)
  {
  case 0:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(false);
    m_Controls->m_IterationsFrame->setVisible(false);
    m_Controls->m_LambdaFrame->setVisible(false);
    break;
  case 1:
    m_Controls->m_DstarFrame->setVisible(true);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(false);
    m_Controls->m_IterationsFrame->setVisible(false);
    m_Controls->m_LambdaFrame->setVisible(false);
    break;
  case 2:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(true);
    m_Controls->m_IterationsFrame->setVisible(false);
    m_Controls->m_LambdaFrame->setVisible(false);
    break;
  case 3:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(true);
    m_Controls->m_NeglBframe->setVisible(true);
    m_Controls->m_IterationsFrame->setVisible(false);
    m_Controls->m_LambdaFrame->setVisible(false);
    break;
  case 4:
    m_Controls->m_DstarFrame->setVisible(false);
    m_Controls->m_NeglSiFrame->setVisible(false);
    m_Controls->m_NeglBframe->setVisible(false);
    m_Controls->m_IterationsFrame->setVisible(false);
    m_Controls->m_LambdaFrame->setVisible(false);
    break;
  }

  OnSliceChanged();
}

void QmitkIVIMView::DStarSlider (int val)
{
  QString sval = QString::number(val/1000.0);
  m_Controls->m_DStarLabel->setText(sval);

  OnSliceChanged();
}

void QmitkIVIMView::BThreshSlider (int val)
{
  QString sval = QString::number(val*5.0);
  m_Controls->m_BThreshLabel->setText(sval);

  OnSliceChanged();
}

void QmitkIVIMView::S0ThreshSlider (int val)
{
  QString sval = QString::number(val*0.5);
  m_Controls->m_S0ThreshLabel->setText(sval);

  OnSliceChanged();
}

void QmitkIVIMView::NumItsSlider (int val)
{
  QString sval = QString::number(val);
  m_Controls->m_NumItsLabel->setText(sval);

  OnSliceChanged();
}

void QmitkIVIMView::LambdaSlider (int val)
{
  QString sval = QString::number(val*.00001);
  m_Controls->m_LambdaLabel->setText(sval);

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

void QmitkIVIMView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& )
{
//  UpdateGui();
}

void QmitkIVIMView::AutoThreshold()
{
  if (m_Controls->m_DwiBox->GetSelectedNode().IsNull())
  {
    // Nothing selected. Inform the user and return
    QMessageBox::information( nullptr, "Template", "Please load and select a diffusion image before starting image processing.");
    return;
  }

  mitk::Image* dimg = dynamic_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());

  if (!dimg)
  {
    // Nothing selected. Inform the user and return
    QMessageBox::information( nullptr, "Template", "No valid diffusion image was found.");
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

  typedef itk::Image<short,3> ImgType;
  ImgType::Pointer img = ImgType::New();
  mitk::CastToItkImage(dimg, img);

  itk::ImageRegionIterator<ImgType> itw (img, img->GetLargestPossibleRegion() );
  itw.GoToBegin();

  itk::ImageRegionConstIterator<VecImgType> itr (vecimg, vecimg->GetLargestPossibleRegion() );
  itr.GoToBegin();

  while(!itr.IsAtEnd())
  {
    itw.Set(itr.Get().GetElement(index));
    ++itr;
    ++itw;
  }

  typedef itk::Statistics::ScalarImageToHistogramGenerator< ImgType >
      HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;

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

  int sliderPos = maxValue * 2;
  m_Controls->m_S0ThreshSlider->setValue(sliderPos);
  S0ThreshSlider(sliderPos);
}

void QmitkIVIMView::FittIVIMStart()
{
  if (m_Controls->m_DwiBox->GetSelectedNode().IsNull())
  {
    QMessageBox::information( nullptr, "Template", "No valid diffusion-weighted image selected.");
    return;
  }
  mitk::Image* img = dynamic_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());


  VecImgType::Pointer vecimg = VecImgType::New();
  mitk::CastToItkImage(img, vecimg);

  OutImgType::IndexType dummy;

  if( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
  {
    // KURTOSIS
    KurtosisFilterType::Pointer filter = KurtosisFilterType::New();
    filter->SetInput(vecimg);
    filter->SetReferenceBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(img));
    filter->SetGradientDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(img));
    filter->SetSmoothingSigma( m_Controls->m_SigmaSpinBox->value() );

    if( m_Controls->m_UseKurtosisBoundsCB->isChecked() )
      filter->SetBoundariesForKurtosis( m_Controls->m_KurtosisRangeWidget->minimumValue(), m_Controls->m_KurtosisRangeWidget->maximumValue() );

    filter->SetFittingScale( static_cast<itk::FitScale>(m_Controls->m_KurtosisFitScale->currentIndex() ) );

    if( m_Controls->m_MaskBox->GetSelectedNode().IsNotNull() )
    {
      mitk::Image::Pointer maskImg = dynamic_cast<mitk::Image*>(m_Controls->m_MaskBox->GetSelectedNode()->GetData());
      typedef itk::Image<short, 3> MaskImgType;

      MaskImgType::Pointer maskItk;
      CastToItkImage( maskImg, maskItk );

      filter->SetImageMask( maskItk );
    }

    filter->Update();

    mitk::LookupTable::Pointer kurt_map_lut = mitk::LookupTable::New();
    kurt_map_lut->SetType( mitk::LookupTable::JET );
    mitk::LookupTableProperty::Pointer kurt_lut_prop =
        mitk::LookupTableProperty::New();
    kurt_lut_prop->SetLookupTable( kurt_map_lut );

    mitk::Image::Pointer dimage = mitk::Image::New();
    dimage->InitializeByItk( filter->GetOutput(0) );
    dimage->SetVolume( filter->GetOutput(0)->GetBufferPointer());

    mitk::Image::Pointer kimage = mitk::Image::New();
    kimage->InitializeByItk( filter->GetOutput(1) );
    kimage->SetVolume( filter->GetOutput(1)->GetBufferPointer());

    QString new_dname = "Kurtosis_DMap";
    new_dname.append("_Method-"+m_Controls->m_KurtosisFitScale->currentText());
    QString new_kname = "Kurtosis_KMap";
    new_kname.append("_Method-"+m_Controls->m_KurtosisFitScale->currentText());

    if( m_Controls->m_CheckKurtD->isChecked() )
    {
      mitk::DataNode::Pointer dnode = mitk::DataNode::New();
      dnode->SetData( dimage );
      dnode->SetName(new_dname.toLatin1());
      dnode->SetProperty("LookupTable", kurt_lut_prop );
      GetDataStorage()->Add(dnode, m_Controls->m_DwiBox->GetSelectedNode());
    }

    if( m_Controls->m_CheckKurtK->isChecked() )
    {
      mitk::DataNode::Pointer knode = mitk::DataNode::New();
      knode->SetData( kimage );
      knode->SetName(new_kname.toLatin1());
      knode->SetProperty("LookupTable", kurt_lut_prop );
      GetDataStorage()->Add(knode, m_Controls->m_DwiBox->GetSelectedNode());
    }

  }
  else
  {
    FittIVIM(vecimg, mitk::DiffusionPropertyHelper::GetGradientContainer(img), mitk::DiffusionPropertyHelper::GetReferenceBValue(img), true, dummy);

    OutputToDatastorage(m_Controls->m_DwiBox->GetSelectedNode());
  }
}

void QmitkIVIMView::OnKurtosisParamsChanged()
{
  OnSliceChanged();
}

void QmitkIVIMView::OnSliceChanged()
{
  if(m_HoldUpdate || !m_Visible)
    return;

  m_Controls->m_Warning->setVisible(false);
  if(m_Controls->m_DwiBox->GetSelectedNode().IsNull())
    return;

  mitk::Image::Pointer diffusionImg = dynamic_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());
  mitk::Image::Pointer maskImg = nullptr;
  if (m_Controls->m_MaskBox->GetSelectedNode().IsNotNull())
    maskImg = dynamic_cast<mitk::Image*>(m_Controls->m_MaskBox->GetSelectedNode()->GetData());

  if (!this->GetRenderWindowPart())
    return;

  if (!m_ListenerActive)
  {
    m_SliceChangeListener.RenderWindowPartActivated(this->GetRenderWindowPart());
    connect(&m_SliceChangeListener, SIGNAL(SliceChanged()), this, SLOT(OnSliceChanged()));
    m_ListenerActive = true;
  }

  VecImgType::Pointer vecimg = VecImgType::New();
  mitk::CastToItkImage(diffusionImg, vecimg);

  VecImgType::Pointer roiImage = VecImgType::New();

  if(maskImg.IsNull())
  {
    int roisize = 0;
    if(m_Controls->m_MethodCombo->currentIndex() == 4)
      roisize = 3;

    mitk::Point3D pos = this->GetRenderWindowPart()->GetSelectedPosition();

    VecImgType::IndexType crosspos;
    diffusionImg->GetGeometry()->WorldToIndex(pos, crosspos);
    if (!vecimg->GetLargestPossibleRegion().IsInside(crosspos))
    {
      m_Controls->m_Warning->setText(QString("Crosshair position not inside of selected diffusion weighted image. Reinit needed!"));
      m_Controls->m_Warning->setVisible(true);
      return;
    }
    else
      m_Controls->m_Warning->setVisible(false);

    VecImgType::IndexType index;
    index[0] = crosspos[0] - roisize; index[0] = index[0] < 0 ? 0 : index[0];
    index[1] = crosspos[1] - roisize; index[1] = index[1] < 0 ? 0 : index[1];
    index[2] = crosspos[2] - roisize; index[2] = index[2] < 0 ? 0 : index[2];

    VecImgType::SizeType size;
    size[0] = roisize*2+1;
    size[1] = roisize*2+1;
    size[2] = roisize*2+1;

    VecImgType::SizeType maxSize = vecimg->GetLargestPossibleRegion().GetSize();
    size[0] = index[0]+size[0] > maxSize[0] ? maxSize[0]-index[0] : size[0];
    size[1] = index[1]+size[1] > maxSize[1] ? maxSize[1]-index[1] : size[1];
    size[2] = index[2]+size[2] > maxSize[2] ? maxSize[2]-index[2] : size[2];

    VecImgType::RegionType region;
    region.SetSize( size );
    region.SetIndex( index );
    vecimg->SetRequestedRegion( region );

    VecImgType::IndexType  newstart;
    newstart.Fill(0);

    VecImgType::RegionType newregion;
    newregion.SetSize( size );
    newregion.SetIndex( newstart );

    roiImage->CopyInformation( vecimg );
    roiImage->SetRegions( newregion );
    roiImage->SetOrigin( pos );
    roiImage->Allocate();
//    roiImage->SetPixel(newstart, vecimg->GetPixel(index));

    typedef itk::ImageRegionIterator<VecImgType> VectorIteratorType;
    VectorIteratorType vecit(vecimg, vecimg->GetRequestedRegion() );
    vecit.GoToBegin();

    typedef itk::ImageRegionIterator<VecImgType> VectorIteratorType;
    VectorIteratorType vecit2(roiImage, roiImage->GetLargestPossibleRegion() );
    vecit2.GoToBegin();

    while( !vecit.IsAtEnd() )
    {
      vecit2.Set(vecit.Get());

      ++vecit;
      ++vecit2;
    }

    if( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
    {
      FitKurtosis(roiImage,
                  mitk::DiffusionPropertyHelper::GetGradientContainer(diffusionImg),
                  mitk::DiffusionPropertyHelper::GetReferenceBValue(diffusionImg),
                  newstart);
    }
    else
    {
      FittIVIM(roiImage,
               mitk::DiffusionPropertyHelper::GetGradientContainer(diffusionImg),
               mitk::DiffusionPropertyHelper::GetReferenceBValue(diffusionImg),
               false,
               crosspos);
    }
  }
  else
  {
    typedef itk::Image<float,3> MaskImgType;

    MaskImgType::Pointer maskItk;
    CastToItkImage( maskImg, maskItk );

    mitk::Point3D pos;
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;

    VecImgType::IndexType index;
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;

    VecImgType::SizeType size;
    size[0] = 1;
    size[1] = 1;
    size[2] = 1;

    VecImgType::RegionType region;
    region.SetSize( size );
    region.SetIndex( index );
    vecimg->SetRequestedRegion( region );

    // iterators over output and input
    itk::ImageRegionConstIteratorWithIndex<VecImgType>
        vecit(vecimg, vecimg->GetLargestPossibleRegion());

    itk::VariableLengthVector<double> avg(vecimg->GetVectorLength());
    avg.Fill(0);

    float numPixels = 0;
    while ( ! vecit.IsAtEnd() )
    {
      VecImgType::PointType point;
      vecimg->TransformIndexToPhysicalPoint(vecit.GetIndex(), point);

      MaskImgType::IndexType index;
      maskItk->TransformPhysicalPointToIndex(point, index);

      if(maskItk->GetPixel(index) != 0)
      {
        avg += vecit.Get();
        numPixels += 1.0;
      }

      // update iterators
      ++vecit;

    }

    avg /= numPixels;

    m_Controls->m_Warning->setText(QString("Averaging ")+QString::number((int)numPixels)+QString(" voxels!"));
    m_Controls->m_Warning->setVisible(true);

    roiImage->CopyInformation( vecimg );
    roiImage->SetRegions( region );
    roiImage->SetOrigin( pos );
    roiImage->Allocate();
    roiImage->SetPixel(index, avg);

    if( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
    {
      FitKurtosis(roiImage,
                  mitk::DiffusionPropertyHelper::GetGradientContainer(diffusionImg),
                  mitk::DiffusionPropertyHelper::GetReferenceBValue(diffusionImg),
                  index);
    }
    else
    {
      FittIVIM(roiImage,
               mitk::DiffusionPropertyHelper::GetGradientContainer(diffusionImg),
               mitk::DiffusionPropertyHelper::GetReferenceBValue(diffusionImg),
               false,
               index);
    }

    // do not update until selection changed, the values will remain the same as long as the mask is selected!
    m_HoldUpdate = true;
  }

  vecimg->SetRegions( vecimg->GetLargestPossibleRegion() );
}

bool QmitkIVIMView::FitKurtosis( itk::VectorImage<short, 3> *vecimg, DirContainerType::ConstPointer dirs, float bval, OutImgType::IndexType &crosspos )
{
  KurtosisFilterType::Pointer filter = KurtosisFilterType::New();

  itk::KurtosisFitConfiguration fit_config;
  fit_config.omit_bzero =  m_Controls->m_OmitBZeroCB->isChecked();
  if( m_Controls->m_UseKurtosisBoundsCB->isChecked() )
  {
    fit_config.use_K_limits = true;
    vnl_vector_fixed<double, 2> k_limits;
    k_limits[0] = m_Controls->m_KurtosisRangeWidget->minimumValue();
    k_limits[1] = m_Controls->m_KurtosisRangeWidget->maximumValue();

    fit_config.K_limits = k_limits;
  }
  fit_config.fit_scale = static_cast<itk::FitScale>(m_Controls->m_KurtosisFitScale->currentIndex() );

  m_KurtosisSnap = filter->GetSnapshot( vecimg->GetPixel( crosspos ), dirs, bval, fit_config );

  QString param_label_text("K=%2, D=%1");
  param_label_text = param_label_text.arg( m_KurtosisSnap.m_K, 4);
  param_label_text = param_label_text.arg( m_KurtosisSnap.m_D, 4);
  m_Controls->m_FittedParamsLabel->setText(param_label_text);

  const double maxb = m_KurtosisSnap.bvalues.max_value();
  double S0 = m_KurtosisSnap.m_Bzero;
  if (m_KurtosisSnap.m_fittedBZero)
    S0 = m_KurtosisSnap.m_BzeroFit;

  std::vector< std::pair<double, double> > d_line;
  d_line.emplace_back(0, S0);
  d_line.emplace_back(maxb, S0*exp(-maxb * m_KurtosisSnap.m_D));
  m_Controls->m_ChartWidget->UpdateData2D(d_line, "D-part of fitted model");

  const unsigned int num_samples = 50;
  std::vector< std::pair<double, double> > y;

  for( unsigned int i=0; i<=num_samples; ++i)
  {
    double b = (((1.0)*i)/(1.0*num_samples))*maxb;
    y.emplace_back(b, S0 * exp( -b * m_KurtosisSnap.m_D + b*b * m_KurtosisSnap.m_D * m_KurtosisSnap.m_D * m_KurtosisSnap.m_K / 6.0 ));
  }
  m_Controls->m_ChartWidget->UpdateData2D(y, "fitted model");

  std::vector< std::pair<double, double> > y_meas;
  for (unsigned int i=0; i<m_KurtosisSnap.measurements.size(); ++i)
  {
    if (!m_Controls->m_OmitBZeroCB->isChecked() || m_KurtosisSnap.bvalues[i] > 0.01)
      y_meas.emplace_back(m_KurtosisSnap.bvalues[i], m_KurtosisSnap.measurements[i]);
  }

  m_Controls->m_ChartWidget->UpdateData2D(y_meas, "signal values");

  return true;
}


bool QmitkIVIMView::FittIVIM(itk::VectorImage<short,3>* vecimg, DirContainerType::ConstPointer dirs, float bval, bool multivoxel, OutImgType::IndexType &crosspos)
{
  IVIMFilterType::Pointer filter = IVIMFilterType::New();
  filter->SetInput(vecimg);
  filter->SetGradientDirections(dirs);
  filter->SetBValue(bval);

  switch(m_Controls->m_MethodCombo->currentIndex())
  {

  case 0:
    filter->SetMethod(IVIMFilterType::IVIM_FIT_ALL);
    filter->SetS0Thres(m_Controls->m_S0ThreshLabel->text().toDouble());
    break;

  case 1:
    filter->SetMethod(IVIMFilterType::IVIM_DSTAR_FIX);
    filter->SetDStar(m_Controls->m_DStarLabel->text().toDouble());
    filter->SetS0Thres(m_Controls->m_S0ThreshLabel->text().toDouble());
    break;

  case 2:
    filter->SetMethod(IVIMFilterType::IVIM_D_THEN_DSTAR);
    filter->SetBThres(m_Controls->m_BThreshLabel->text().toDouble());
    filter->SetS0Thres(m_Controls->m_S0ThreshLabel->text().toDouble());
    filter->SetFitDStar(m_Controls->m_CheckDStar->isChecked());
    break;

  case 3:
    filter->SetMethod(IVIMFilterType::IVIM_LINEAR_D_THEN_F);
    filter->SetBThres(m_Controls->m_BThreshLabel->text().toDouble());
    filter->SetS0Thres(m_Controls->m_S0ThreshLabel->text().toDouble());
    filter->SetFitDStar(m_Controls->m_CheckDStar->isChecked());
    break;

  case 4:
    filter->SetMethod(IVIMFilterType::IVIM_REGULARIZED);
    filter->SetBThres(m_Controls->m_BThreshLabel->text().toDouble());
    filter->SetS0Thres(m_Controls->m_S0ThreshLabel->text().toDouble());
    filter->SetNumberIterations(m_Controls->m_NumItsLabel->text().toInt());
    filter->SetLambda(m_Controls->m_LambdaLabel->text().toDouble());
    filter->SetFitDStar(m_Controls->m_CheckDStar->isChecked());
    break;
  }

  if(!multivoxel)
  {
    filter->SetFitDStar(true);
  }

  filter->SetNumberOfThreads(1);
  filter->SetVerbose(false);
  filter->SetCrossPosition(crosspos);

  try{
    filter->Update();
    m_IvimSnap = filter->GetSnapshot();
    m_DStarMap = filter->GetOutput(2);
    m_DMap = filter->GetOutput(1);
    m_fMap = filter->GetOutput();

    QString param_label_text("f=%1, D=%2, D*=%3");
    param_label_text = param_label_text.arg(m_IvimSnap.currentF,4);
    param_label_text = param_label_text.arg(m_IvimSnap.currentD,4);
    param_label_text = param_label_text.arg(m_IvimSnap.currentDStar,4);
    m_Controls->m_FittedParamsLabel->setText(param_label_text);

    double maxb = m_IvimSnap.bvalues.max_value();

    std::vector< std::pair<double, double> > d_line;
    d_line.emplace_back(0, 1-m_IvimSnap.currentFunceiled);
    d_line.emplace_back(maxb, d_line[0].second*exp(-maxb * m_IvimSnap.currentD));
    m_Controls->m_ChartWidget->UpdateData2D(d_line, "D-part of fitted model");

    std::vector< std::pair<double, double> > y;
    int nsampling = 50;
    double f = 1-m_IvimSnap.currentFunceiled;
    for(int i=0; i<=nsampling; i++)
    {
      double x = (((1.0)*i)/(1.0*nsampling))*maxb;
      y.emplace_back(x, f*exp(- x * m_IvimSnap.currentD) + (1-f)*exp(- x * (m_IvimSnap.currentD+m_IvimSnap.currentDStar)));
    }
    m_Controls->m_ChartWidget->UpdateData2D(y, "fitted model");

    std::vector< std::pair<double, double> > y_meas;
    for (unsigned int i=0; i<m_IvimSnap.meas1.size(); ++i)
      y_meas.emplace_back(m_IvimSnap.bvals1[i], m_IvimSnap.meas1[i]);

    m_Controls->m_ChartWidget->UpdateData2D(y_meas, "signal values");

    if(m_IvimSnap.bvals2.size() > 0)
    {
      AddSecondFitPlot();
      std::vector< std::pair<double, double> > additonal_meas;
      for (int i=0; i<static_cast<int>(m_IvimSnap.bvals2.size()); ++i)
        additonal_meas.emplace_back(m_IvimSnap.bvals2[i], m_IvimSnap.meas2[i]);

      m_Controls->m_ChartWidget->UpdateData2D(additonal_meas, "signal values (used for second fit)");
    }
    else
    {
      RemoveSecondFitPlot();
    }
  }
  catch (itk::ExceptionObject &ex)
  {
    MITK_INFO << ex ;
    m_Controls->m_Warning->setText(QString("IVIM fit not possible: ")+ex.GetDescription());
    m_Controls->m_Warning->setVisible(true);
    return false;
  }
  return true;
}

void QmitkIVIMView::OutputToDatastorage(mitk::DataNode::Pointer node)
{
  mitk::LookupTable::Pointer lut = mitk::LookupTable::New();
  lut->SetType( mitk::LookupTable::JET );
  mitk::LookupTableProperty::Pointer lut_prop = mitk::LookupTableProperty::New();
  lut_prop->SetLookupTable( lut );

  if(m_Controls->m_CheckDStar->isChecked())
  {
    mitk::Image::Pointer dstarimage = mitk::Image::New();
    dstarimage->InitializeByItk(m_DStarMap.GetPointer());
    dstarimage->SetVolume(m_DStarMap->GetBufferPointer());
    QString newname2 = ""; newname2 = newname2.append("D* Fit-%1").arg(m_Controls->m_MethodCombo->currentIndex()+1);
    mitk::DataNode::Pointer node2=mitk::DataNode::New();
    node2->SetData( dstarimage );
    node2->SetName(newname2.toLatin1());
    node2->SetProperty("LookupTable", lut_prop );

    mitk::LevelWindowProperty::Pointer levWinProp = mitk::LevelWindowProperty::New();
    mitk::LevelWindow levelwindow;
    levelwindow.SetRangeMinMax(0, 0.2);
    levWinProp->SetLevelWindow(levelwindow);
    node2->SetProperty("levelwindow", levWinProp);

    GetDataStorage()->Add(node2, node);
  }

  if(m_Controls->m_CheckD->isChecked())
  {
    mitk::Image::Pointer dimage = mitk::Image::New();
    dimage->InitializeByItk(m_DMap.GetPointer());
    dimage->SetVolume(m_DMap->GetBufferPointer());
    QString newname1 = ""; newname1 = newname1.append("D Fit-%1").arg(m_Controls->m_MethodCombo->currentIndex()+1);
    mitk::DataNode::Pointer node1=mitk::DataNode::New();
    node1->SetData( dimage );
    node1->SetName(newname1.toLatin1());
    node1->SetProperty("LookupTable", lut_prop );

    mitk::LevelWindowProperty::Pointer levWinProp = mitk::LevelWindowProperty::New();
    mitk::LevelWindow levelwindow;
    levelwindow.SetRangeMinMax(0, 0.003);
    levWinProp->SetLevelWindow(levelwindow);
    node1->SetProperty("levelwindow", levWinProp);

    GetDataStorage()->Add(node1, node);
  }

  if(m_Controls->m_Checkf->isChecked())
  {
    mitk::Image::Pointer image = mitk::Image::New();
    image->InitializeByItk(m_fMap.GetPointer());
    image->SetVolume(m_fMap->GetBufferPointer());
    QString newname0 = ""; newname0 = newname0.append("f Fit-%1").arg(m_Controls->m_MethodCombo->currentIndex()+1);
    mitk::DataNode::Pointer node3=mitk::DataNode::New();
    node3->SetData( image );
    node3->SetName(newname0.toLatin1());
    node3->SetProperty("LookupTable", lut_prop );

    mitk::LevelWindowProperty::Pointer levWinProp = mitk::LevelWindowProperty::New();
    mitk::LevelWindow levelwindow;
    levelwindow.SetRangeMinMax(0, 1);
    levWinProp->SetLevelWindow(levelwindow);
    node3->SetProperty("levelwindow", levWinProp);

    GetDataStorage()->Add(node3, node);
  }

  this->GetRenderWindowPart()->RequestUpdate();
}

void QmitkIVIMView::ClipboardCurveButtonClicked()
{
  // Kurtosis
  if ( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
  {
    std::stringstream ss;
    QString clipboard("Measurement Points\n");

    ss << m_KurtosisSnap.bvalues << "\n" << m_KurtosisSnap.measurements << "\n\n";
    ss << "Fitted Values ( D  K  [b_0] ) \n" << m_KurtosisSnap.m_D << " " << m_KurtosisSnap.m_K;

    if( m_KurtosisSnap.m_fittedBZero )
      ss << " " << m_KurtosisSnap.m_BzeroFit;

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
    ss << m_IvimSnap.bvalues << "\n" << m_IvimSnap.allmeas << "\n\n";
    ss << "Fitted Values ( f  D  D* ) \n" << m_IvimSnap.currentF << " " << m_IvimSnap.currentD << " " << m_IvimSnap.currentDStar;


    ss << "\n\n";
    clipboard.append( QString( ss.str().c_str() ));

    ss.str( std::string() );
    ss.clear();

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
}

void QmitkIVIMView::SavePlotButtonClicked()
{
  m_Controls->m_ChartWidget->SavePlotAsImage();
}


void QmitkIVIMView::ClipboardStatisticsButtonClicked()
{
  // Kurtosis
  if ( m_Controls->m_ModelTabSelectionWidget->currentIndex() )
  {

    QString clipboard( "D \t K \n" );
    clipboard = clipboard.append( "%L1 \t %L2" )
        .arg( m_KurtosisSnap.m_D, 0, 'f', 10 )
        .arg( m_KurtosisSnap.m_K, 0, 'f', 10 ) ;

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
  else
  {
    QString clipboard( "f \t D \t D* \n" );
    clipboard = clipboard.append( "%L1 \t %L2 \t %L3" )
        .arg( m_IvimSnap.currentF, 0, 'f', 10 )
        .arg( m_IvimSnap.currentD, 0, 'f', 10 )
        .arg( m_IvimSnap.currentDStar, 0, 'f', 10 ) ;

    QApplication::clipboard()->setText(clipboard, QClipboard::Clipboard );
  }
}

void QmitkIVIMView::Activated()
{
  m_Active = true;
}

void QmitkIVIMView::Deactivated()
{
  m_Active = false;
}

void QmitkIVIMView::Visible()
{
  m_Visible = true;

  if (this->GetRenderWindowPart() && !m_ListenerActive)
  {
    m_SliceChangeListener.RenderWindowPartActivated(this->GetRenderWindowPart());
    connect(&m_SliceChangeListener, SIGNAL(SliceChanged()), this, SLOT(OnSliceChanged()));
    m_ListenerActive = true;
  }
}

void QmitkIVIMView::Hidden()
{
  m_Visible = false;
}
