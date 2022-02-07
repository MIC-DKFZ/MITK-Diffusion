/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

#include "QmitkTractometryView.h"

#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateAnd.h>
#include <mitkImageCast.h>
#include <mitkImageAccessByItk.h>
#include <mitkImage.h>
#include <mitkFiberBundle.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
#include <mitkPixelTypeMultiplex.h>
#include <mitkImagePixelReadAccessor.h>
#include <berryIPreferences.h>
#include <berryWorkbenchPlugin.h>
#include <berryQtPreferences.h>
#include <vtkLookupTable.h>
#include <QClipboard>
#include <mitkLexicalCast.h>
#include <QmitkChartWidget.h>
#include <mitkLookupTable.h>
#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricEuclideanStd.h>
#include <itkTractDensityImageFilter.h>
#include <itkTractParcellationFilter.h>
#include <mitkLookupTableProperty.h>
#include <mitkLevelWindowProperty.h>
#include <itkMaskedStatisticsImageFilter.h>
#include <mitkDiffusionFunctionCollection.h>
#include <mitkTractometry.h>


const std::string QmitkTractometryView::VIEW_ID = "org.mitk.views.tractometry";
using namespace mitk;

QmitkTractometryView::QmitkTractometryView()
  : QmitkAbstractView()
  , m_Controls( nullptr )
  , m_Visible(false)
{

}

// Destructor
QmitkTractometryView::~QmitkTractometryView()
{

}

void QmitkTractometryView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkTractometryViewControls;
    m_Controls->setupUi( parent );

    connect( m_Controls->m_MethodBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()) );
    connect( m_Controls->m_StartButton, SIGNAL(clicked()), this, SLOT(StartTractometry()) );

    mitk::TNodePredicateDataType<mitk::Image>::Pointer imageP = mitk::TNodePredicateDataType<mitk::Image>::New();
    mitk::NodePredicateDimension::Pointer dimP = mitk::NodePredicateDimension::New(3);

    m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateAnd::New(imageP, dimP));

    m_Controls->m_ChartWidget->SetXAxisLabel("Tract position");
    m_Controls->m_ChartWidget->SetYAxisLabel("Image Value");

    m_Controls->m_ChartWidget->SetTheme(QmitkChartWidget::ColorTheme::darkstyle);
  }
}

void QmitkTractometryView::SetFocus()
{
}

void QmitkTractometryView::UpdateGui()
{
  berry::IWorkbenchPart::Pointer nullPart;
  OnSelectionChanged(nullPart, QList<mitk::DataNode::Pointer>(m_CurrentSelection));
}

void QmitkTractometryView::StaticResamplingTractometry(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector<std::vector<double> > &data, std::string& clipboard_string)
{
  itk::Image<float, 3>::Pointer itkImage = itk::Image<float, 3>::New();
  CastToItkImage(image, itkImage);

  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  unsigned int num_points = m_NumSamplingPoints;
  mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();

  vnl_matrix<float> output = mitk::Tractometry::StaticResamplingTractometry(itkImage, working_fib, num_points, m_ReferenceFib);

  std::vector< double > std_values1;
  std::vector< double > std_values2;
  std::vector< double > mean_values;

  for (unsigned int i=0; i<output.rows(); ++i)
  {
    auto row = output.get_row(i);
    float mean = row.mean();
    double stdev = 0;

    for (unsigned int j=0; j<row.size(); ++j)
    {
      double diff = mean - row.get(j);
      diff *= diff;
      stdev += diff;
    }
    stdev /= row.size();
    stdev = std::sqrt(stdev);

    clipboard_string += boost::lexical_cast<std::string>(mean);
    clipboard_string += " ";
    clipboard_string += boost::lexical_cast<std::string>(stdev);
    clipboard_string += "\n";

    mean_values.push_back(mean);
    std_values1.push_back(mean + stdev);
    std_values2.push_back(mean - stdev);
  }
  clipboard_string += "\n";

  data.push_back(mean_values);
  data.push_back(std_values1);
  data.push_back(std_values2);

  if (m_Controls->m_ShowBinned->isChecked())
  {
    mitk::DataNode::Pointer new_node = mitk::DataNode::New();
    new_node->SetData(working_fib);
    new_node->SetName("binned_static");
    new_node->SetVisibility(true);
    node->SetVisibility(false);
    GetDataStorage()->Add(new_node, node);
  }
}

void QmitkTractometryView::NearestCentroidPointTractometry(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector< std::vector< double > >& data, std::string& clipboard_string)
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  unsigned int num_points = m_NumSamplingPoints;

  mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();
  working_fib->ResampleSpline(1.0);

  itk::Image<float, 3>::Pointer itkImage = itk::Image<float, 3>::New();
  CastToItkImage(image, itkImage);

  auto output = mitk::Tractometry::NearestCentroidPointTractometry(itkImage, working_fib, num_points, m_Controls->m_MaxCentroids->value(), m_Controls->m_ClusterSize->value(), m_ReferenceFib);

  std::vector< double > std_values1;
  std::vector< double > std_values2;
  std::vector< double > mean_values;

  for (auto row : output)
  {
    float mean = row.mean();
    double stdev = 0;

    for (unsigned int j=0; j<row.size(); ++j)
    {
      double diff = mean - row.get(j);
      diff *= diff;
      stdev += diff;
    }
    stdev /= row.size();
    stdev = std::sqrt(stdev);

    clipboard_string += boost::lexical_cast<std::string>(mean);
    clipboard_string += " ";
    clipboard_string += boost::lexical_cast<std::string>(stdev);
    clipboard_string += "\n";

    mean_values.push_back(mean);
    std_values1.push_back(mean + stdev);
    std_values2.push_back(mean - stdev);
  }
  clipboard_string += "\n";

  data.push_back(mean_values);
  data.push_back(std_values1);
  data.push_back(std_values2);

  if (m_Controls->m_ShowBinned->isChecked())
  {
    mitk::DataNode::Pointer new_node = mitk::DataNode::New();
    new_node->SetData(working_fib);
    new_node->SetName("binned_centroid");
    new_node->SetVisibility(true);
    node->SetVisibility(false);
    GetDataStorage()->Add(new_node, node);
  }
}

void QmitkTractometryView::Activated()
{

}

void QmitkTractometryView::Deactivated()
{

}

void QmitkTractometryView::Visible()
{
  m_Visible = true;
  QList<mitk::DataNode::Pointer> selection = GetDataManagerSelection();
  berry::IWorkbenchPart::Pointer nullPart;
  OnSelectionChanged(nullPart, selection);
}

void QmitkTractometryView::Hidden()
{
  m_Visible = false;
}

std::string QmitkTractometryView::RGBToHexString(double *rgb)
{
  std::ostringstream os;
  for (int i = 0; i < 3; ++i)
    {
    os << std::setw(2) << std::setfill('0') << std::hex
       << static_cast<int>(rgb[i] * 255);
    }
  return os.str();
}

void QmitkTractometryView::AlongTractRadiomicsPreprocessing(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector< std::vector< double > >& data, std::string& clipboard_string)
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  // calculate mask
  typedef itk::Image<unsigned char, 3> ParcellationImageType;
  ParcellationImageType::Pointer itkImage = ParcellationImageType::New();
  CastToItkImage(image, itkImage);

  itk::TractParcellationFilter< >::Pointer parcellator = itk::TractParcellationFilter< >::New();
  parcellator->SetInputImage(itkImage);
  parcellator->SetNumParcels(m_NumSamplingPoints);
  parcellator->SetInputTract(fib);
  parcellator->SetNumCentroids(m_Controls->m_MaxCentroids->value());
  parcellator->SetStartClusterSize(m_Controls->m_ClusterSize->value());
  parcellator->Update();
  ParcellationImageType::Pointer out_image = parcellator->GetOutput(0);
  ParcellationImageType::Pointer out_image_pp = parcellator->GetOutput(1);
  auto binary_masks = parcellator->GetBinarySplit(out_image_pp);

  mitk::Image::Pointer seg_img = mitk::Image::New();
  seg_img->InitializeByItk(out_image.GetPointer());
  seg_img->SetVolume(out_image->GetBufferPointer());

  mitk::Image::Pointer seg_img_pp = mitk::Image::New();
  seg_img_pp->InitializeByItk(out_image_pp.GetPointer());
  seg_img_pp->SetVolume(out_image_pp->GetBufferPointer());

  std::vector< double > std_values1;
  std::vector< double > std_values2;
  std::vector< double > mean_values;
  for (auto mask : binary_masks)
  {
    itk::Image<float, 3>::Pointer data_image = itk::Image<float, 3>::New();
    CastToItkImage(image, data_image);
    itk::MaskedStatisticsImageFilter<itk::Image<float, 3>>::Pointer statisticsImageFilter = itk::MaskedStatisticsImageFilter<itk::Image<float, 3>>::New();
    statisticsImageFilter->SetInput(data_image);
    statisticsImageFilter->SetMask(mask);
    statisticsImageFilter->Update();
    double mean = statisticsImageFilter->GetMean();
    double stdev = std::sqrt(statisticsImageFilter->GetVariance());

    std_values1.push_back(mean + stdev);
    std_values2.push_back(mean - stdev);
    mean_values.push_back(mean);

    clipboard_string += boost::lexical_cast<std::string>(mean);
    clipboard_string += " ";
    clipboard_string += boost::lexical_cast<std::string>(stdev);
    clipboard_string += "\n";
  }
  clipboard_string += "\n";

  data.push_back(mean_values);
  data.push_back(std_values1);
  data.push_back(std_values2);

  mitk::LookupTable::Pointer lut = mitk::LookupTable::New();
  lut->SetType( mitk::LookupTable::MULTILABEL );
  mitk::LookupTableProperty::Pointer lut_prop = mitk::LookupTableProperty::New();
  lut_prop->SetLookupTable( lut );

  mitk::LevelWindow lw;
  lw.SetRangeMinMax(0, parcellator->GetNumParcels());

//  mitk::DataNode::Pointer new_node = mitk::DataNode::New();
//  new_node->SetData(seg_img);
//  new_node->SetName("tract parcellation");
//  new_node->SetVisibility(true);
//  new_node->SetProperty("LookupTable", lut_prop );
//  new_node->SetProperty( "levelwindow", mitk::LevelWindowProperty::New( lw ) );
//  node->SetVisibility(false);
//  GetDataStorage()->Add(new_node, node);

  mitk::DataNode::Pointer new_node2 = mitk::DataNode::New();
  new_node2->SetData(seg_img_pp);
  new_node2->SetName("tract parcellation pp");
  new_node2->SetVisibility(false);
  new_node2->SetProperty("LookupTable", lut_prop );
  new_node2->SetProperty( "levelwindow", mitk::LevelWindowProperty::New( lw ) );
  GetDataStorage()->Add(new_node2, node);

  mitk::DataNode::Pointer new_node3 = mitk::DataNode::New();
  auto working_fib = fib->GetDeepCopy();
  working_fib->ColorFibersByScalarMap(seg_img, false, false, mitk::LookupTable::LookupTableType::MULTILABEL, 0.9);
  new_node3->SetData(working_fib);
  new_node3->SetName("centroids");
  GetDataStorage()->Add(new_node3, node);
}

void QmitkTractometryView::StartTractometry()
{
  m_ReferenceFib = dynamic_cast<mitk::FiberBundle*>(m_CurrentSelection.at(0)->GetData())->GetDeepCopy();

  mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(m_Controls->m_ImageBox->GetSelectedNode()->GetData());

  MITK_INFO << "Resanmpling reference fibers";
  if (m_Controls->m_SamplingPointsBox->value()<3)
  {
    typedef itk::Image<unsigned char, 3> ParcellationImageType;
    ParcellationImageType::Pointer itkImage = ParcellationImageType::New();
    CastToItkImage(image, itkImage);

    m_NumSamplingPoints = mitk::Tractometry::EstimateNumSamplingPoints(itkImage, m_ReferenceFib, 3);
  }
  else
    m_NumSamplingPoints = m_Controls->m_SamplingPointsBox->value();
  m_ReferenceFib->ResampleToNumPoints(m_NumSamplingPoints);

  double color[3] = {0,0,0};
  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  this->m_Controls->m_ChartWidget->Clear();
  std::string clipboardString = "";
  int c = 1;
  for (auto node : m_CurrentSelection)
  {
    clipboardString += node->GetName() + "\n";
    clipboardString += "mean stdev\n";

    std::vector< std::vector< double > > data;
    switch (m_Controls->m_MethodBox->currentIndex())
    {
    case 0:
    {
      StaticResamplingTractometry( image, node, data, clipboardString );
      break;
    }
    case 1:
    {
      NearestCentroidPointTractometry( image, node, data, clipboardString );
      break;
    }
    case 2:
    {
      AlongTractRadiomicsPreprocessing(image, node, data, clipboardString);
      break;
    }
    default:
    {
      StaticResamplingTractometry( image, node, data, clipboardString );
    }
    }

    m_Controls->m_ChartWidget->AddData1D(data.at(0), node->GetName() + " Mean", QmitkChartWidget::ChartType::line);
    m_Controls->m_ChartWidget->SetLineStyle(node->GetName() + " Mean", QmitkChartWidget::LineStyle::solid);
    if (m_Controls->m_StDevBox->isChecked())
    {
      m_Controls->m_ChartWidget->AddData1D(data.at(1), node->GetName() + " +STDEV", QmitkChartWidget::ChartType::line);
      m_Controls->m_ChartWidget->AddData1D(data.at(2), node->GetName() + " -STDEV", QmitkChartWidget::ChartType::line);
      m_Controls->m_ChartWidget->SetLineStyle(node->GetName() + " +STDEV", QmitkChartWidget::LineStyle::dashed);
      m_Controls->m_ChartWidget->SetLineStyle(node->GetName() + " -STDEV", QmitkChartWidget::LineStyle::dashed);
    }

    lookupTable->GetTableValue(c, color);
    this->m_Controls->m_ChartWidget->SetColor(node->GetName() + " Mean", RGBToHexString(color));

    if (m_Controls->m_StDevBox->isChecked())
    {
      color[0] *= 0.5;
      color[1] *= 0.5;
      color[2] *= 0.5;
      this->m_Controls->m_ChartWidget->SetColor(node->GetName() + " +STDEV", RGBToHexString(color));
      this->m_Controls->m_ChartWidget->SetColor(node->GetName() + " -STDEV", RGBToHexString(color));
    }

    this->m_Controls->m_ChartWidget->Show(true);
    this->m_Controls->m_ChartWidget->SetShowDataPoints(false);
    ++c;
  }

  QApplication::clipboard()->setText(clipboardString.c_str(), QClipboard::Clipboard);
}

void QmitkTractometryView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& nodes)
{
  m_Controls->m_StartButton->setEnabled(false);

  if (!m_Visible)
    return;

  if (m_Controls->m_MethodBox->currentIndex()==0)
    m_Controls->m_ClusterFrame->setVisible(false);
  else
    m_Controls->m_ClusterFrame->setVisible(true);

  m_CurrentSelection.clear();
  if(m_Controls->m_ImageBox->GetSelectedNode().IsNull())
    return;

  for (auto node: nodes)
    if ( dynamic_cast<mitk::FiberBundle*>(node->GetData()) )
      m_CurrentSelection.push_back(node);
  if (!m_CurrentSelection.empty())
    m_Controls->m_StartButton->setEnabled(true);
}
