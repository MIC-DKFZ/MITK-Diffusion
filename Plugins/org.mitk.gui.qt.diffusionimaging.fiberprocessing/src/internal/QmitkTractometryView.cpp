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
#include <itkImageRegionIterator.h>


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

bool QmitkTractometryView::Flip(vtkSmartPointer< vtkPolyData > polydata1, int i, vtkSmartPointer< vtkPolyData > ref_poly)
{
  double d_direct = 0;
  double d_flipped = 0;

  vtkCell* cell1 = polydata1->GetCell(0);
  if (ref_poly!=nullptr)
    cell1 = ref_poly->GetCell(0);
  auto numPoints1 = cell1->GetNumberOfPoints();
  vtkPoints* points1 = cell1->GetPoints();

  std::vector<itk::Point<double, 3>> ref_points;
  for (int j=0; j<numPoints1; ++j)
  {
    double* p1 = points1->GetPoint(j);
    itk::Point<double, 3> itk_p;
    itk_p[0] = p1[0];
    itk_p[1] = p1[1];
    itk_p[2] = p1[2];
    ref_points.push_back(itk_p);
  }

  vtkCell* cell2 = polydata1->GetCell(i);
  vtkPoints* points2 = cell2->GetPoints();

  for (int j=0; j<numPoints1; ++j)
  {
    auto p1 = ref_points.at(j);

    double* p2 = points2->GetPoint(j);
    d_direct = (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]);

    double* p3 = points2->GetPoint(numPoints1-j-1);
    d_flipped = (p1[0]-p3[0])*(p1[0]-p3[0]) + (p1[1]-p3[1])*(p1[1]-p3[1]) + (p1[2]-p3[2])*(p1[2]-p3[2]);
  }

  if (d_direct>d_flipped)
    return true;
  return false;
}

template <typename TPixel>
void QmitkTractometryView::StaticResamplingTractometry(const mitk::PixelType, mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector<std::vector<double> > &data, std::string& clipboard_string)
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  unsigned int num_points = m_Controls->m_SamplingPointsBox->value();
  mitk::ImagePixelReadAccessor<TPixel,3> readimage(image, image->GetVolumeData(0));
  mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();
  working_fib->ResampleToNumPoints(num_points);
  vtkSmartPointer< vtkPolyData > polydata = working_fib->GetFiberPolyData();

  double rgb[3] = {0,0,0};

  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  std::vector<std::vector<double> > all_values;
  std::vector< double > mean_values;
  for (unsigned int i=0; i<num_points; ++i)
    mean_values.push_back(0);

  double min = 100000.0;
  double max = 0;
  double mean = 0;
  for (unsigned int i=0; i<working_fib->GetNumFibers(); ++i)
  {
    vtkCell* cell = polydata->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    std::vector< double > fib_vals;

    bool flip = false;
    if (i>0)
      flip = Flip(polydata, i);
    else if (m_ReferencePolyData!=nullptr)
      flip = Flip(polydata, 0, m_ReferencePolyData);

    for (int j=0; j<numPoints; j++)
    {
      lookupTable->GetTableValue(j, rgb);

      double* p;
      if (flip)
      {
        auto p_idx = numPoints - j - 1;
        p = points->GetPoint(p_idx);

        working_fib->ColorSinglePoint(i, p_idx, rgb);
      }
      else
      {
        p = points->GetPoint(j);

        working_fib->ColorSinglePoint(i, j, rgb);
      }

      Point3D px;
      px[0] = p[0];
      px[1] = p[1];
      px[2] = p[2];
      double pixelValue = static_cast<double>(readimage.GetPixelByWorldCoordinates(px));
      fib_vals.push_back(pixelValue);
      mean += pixelValue;
      if (pixelValue<min)
        min = pixelValue;
      else if (pixelValue>max)
        max = pixelValue;

      mean_values.at(j) += pixelValue;

    }

    all_values.push_back(fib_vals);
  }

  if (m_ReferencePolyData==nullptr)
    m_ReferencePolyData = polydata;

  std::vector< double > std_values1;
  std::vector< double > std_values2;
  for (unsigned int i=0; i<num_points; ++i)
  {
    mean_values.at(i) /= working_fib->GetNumFibers();
    double stdev = 0;

    for (unsigned int j=0; j<all_values.size(); ++j)
    {
      double diff = mean_values.at(i) - all_values.at(j).at(i);
      diff *= diff;
      stdev += diff;
    }
    stdev /= all_values.size();
    stdev = std::sqrt(stdev);
    std_values1.push_back(mean_values.at(i) + stdev);
    std_values2.push_back(mean_values.at(i) - stdev);

    clipboard_string += boost::lexical_cast<std::string>(mean_values.at(i));
    clipboard_string += " ";
    clipboard_string += boost::lexical_cast<std::string>(stdev);
    clipboard_string += "\n";
  }
  clipboard_string += "\n";

  data.push_back(mean_values);
  data.push_back(std_values1);
  data.push_back(std_values2);

  MITK_INFO << "Min: " << min;
  MITK_INFO << "Max: " << max;
  MITK_INFO << "Mean: " << mean/working_fib->GetNumberOfPoints();

  if (m_Controls->m_ShowBinned->isChecked())
  {
    mitk::DataNode::Pointer new_node = mitk::DataNode::New();
    auto children = GetDataStorage()->GetDerivations(node);
    for (unsigned int i=0; i<children->size(); ++i)
    {
      if (children->at(i)->GetName() == "binned_static")
      {
        new_node = children->at(i);
        new_node->SetData(working_fib);
        new_node->SetVisibility(true);
        node->SetVisibility(false);
        mitk::RenderingManager::GetInstance()->RequestUpdateAll();
        return;
      }
    }

    new_node->SetData(working_fib);
    new_node->SetName("binned_static");
    new_node->SetVisibility(true);
    node->SetVisibility(false);
    GetDataStorage()->Add(new_node, node);
  }
}

template <typename TPixel>
void QmitkTractometryView::NearestCentroidPointTractometry(const mitk::PixelType, mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector< std::vector< double > >& data, std::string& clipboard_string)
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  unsigned int num_points = m_Controls->m_SamplingPointsBox->value();
  mitk::ImagePixelReadAccessor<TPixel,3> readimage(image, image->GetVolumeData(0));
  mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();
  working_fib->ResampleSpline(1.0);
  vtkSmartPointer< vtkPolyData > working_polydata = working_fib->GetFiberPolyData();

  // clustering
  std::vector< mitk::ClusteringMetric* > metrics;
  metrics.push_back({new mitk::ClusteringMetricEuclideanStd()});

  mitk::FiberBundle::Pointer fib_static_resampled = fib->GetDeepCopy();
  fib_static_resampled->ResampleToNumPoints(num_points);
  vtkSmartPointer< vtkPolyData > polydata_static_resampled = fib_static_resampled->GetFiberPolyData();

  std::vector<mitk::FiberBundle::Pointer> centroids;
  std::shared_ptr< mitk::TractClusteringFilter > clusterer = std::make_shared<mitk::TractClusteringFilter>();
  int c=0;
  while (c<30 && (centroids.empty() || centroids.size()>static_cast<unsigned long>(m_Controls->m_MaxCentroids->value())))
  {
    float cluster_size = m_Controls->m_ClusterSize->value() + m_Controls->m_ClusterSize->value()*c*0.2;
    float max_d = 0;
    int i=1;
    std::vector< float > distances;
    while (max_d < working_fib->GetGeometry()->GetDiagonalLength()/2)
    {
      distances.push_back(cluster_size*i);
      max_d = cluster_size*i;
      ++i;
    }

    clusterer->SetDistances(distances);
    clusterer->SetTractogram(fib_static_resampled);
    clusterer->SetMetrics(metrics);
    clusterer->SetMergeDuplicateThreshold(cluster_size);
    clusterer->SetDoResampling(false);
    clusterer->SetNumPoints(num_points);
  //  clusterer->SetMaxClusters(m_Controls->m_MaxCentroids->value());
    clusterer->SetMinClusterSize(1);
    clusterer->Update();
    centroids = clusterer->GetOutCentroids();
    ++c;
  }

  double rgb[3] = {0,0,0};
  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  std::vector<std::vector<double> > all_values;
  std::vector< double > mean_values;
  std::vector< unsigned int > value_count;
  for (unsigned int i=0; i<num_points; ++i)
  {
    mean_values.push_back(0);
    value_count.push_back(0);
  }

  double min = 100000.0;
  double max = 0;
  double mean = 0;
  for (unsigned int i=0; i<working_fib->GetNumFibers(); ++i)
  {
    vtkCell* cell = working_polydata->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    std::vector< double > fib_vals;
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);

      int min_bin = 0;
      float d=999999;
      for (auto centroid : centroids)
      {
        auto centroid_polydata = centroid->GetFiberPolyData();

        vtkCell* centroid_cell = centroid_polydata->GetCell(0);
        auto centroid_numPoints = centroid_cell->GetNumberOfPoints();
        vtkPoints* centroid_points = centroid_cell->GetPoints();

        bool centroid_flip = Flip(centroid_polydata, 0, centroids.at(0)->GetFiberPolyData());

        for (int bin=0; bin<centroid_numPoints; ++bin)
        {
          double* centroid_p;
          centroid_p = centroid_points->GetPoint(bin);
          float temp_d = std::sqrt((p[0]-centroid_p[0])*(p[0]-centroid_p[0]) + (p[1]-centroid_p[1])*(p[1]-centroid_p[1]) + (p[2]-centroid_p[2])*(p[2]-centroid_p[2]));
          if (temp_d<d)
          {
            d = temp_d;
            if (centroid_flip)
              min_bin = centroid_numPoints-bin-1;
            else
              min_bin = bin;
          }
        }

      }

      lookupTable->GetTableValue(min_bin, rgb);
      working_fib->ColorSinglePoint(i, j, rgb);

      Point3D px;
      px[0] = p[0];
      px[1] = p[1];
      px[2] = p[2];
      double pixelValue = static_cast<double>(readimage.GetPixelByWorldCoordinates(px));
      fib_vals.push_back(pixelValue);
      mean += pixelValue;
      if (pixelValue<min)
        min = pixelValue;
      else if (pixelValue>max)
        max = pixelValue;

      mean_values.at(min_bin) += pixelValue;
      value_count.at(min_bin) += 1;
    }

    all_values.push_back(fib_vals);
  }

  if (m_ReferencePolyData==nullptr)
    m_ReferencePolyData = working_polydata;

  std::vector< double > std_values1;
  std::vector< double > std_values2;
  for (unsigned int i=0; i<num_points; ++i)
  {
    mean_values.at(i) /= value_count.at(i);
    double stdev = 0;

    for (unsigned int j=0; j<all_values.size(); ++j)
    {
      double diff = mean_values.at(i) - all_values.at(j).at(i);
      diff *= diff;
      stdev += diff;
    }
    stdev /= all_values.size();
    stdev = std::sqrt(stdev);
    std_values1.push_back(mean_values.at(i) + stdev);
    std_values2.push_back(mean_values.at(i) - stdev);

    clipboard_string += boost::lexical_cast<std::string>(mean_values.at(i));
    clipboard_string += " ";
    clipboard_string += boost::lexical_cast<std::string>(stdev);
    clipboard_string += "\n";
  }
  clipboard_string += "\n";

  data.push_back(mean_values);
  data.push_back(std_values1);
  data.push_back(std_values2);

  MITK_INFO << "Min: " << min;
  MITK_INFO << "Max: " << max;
  MITK_INFO << "Mean: " << mean/working_fib->GetNumberOfPoints();

  if (m_Controls->m_ShowBinned->isChecked())
  {
    mitk::DataNode::Pointer new_node = mitk::DataNode::New();
//    mitk::DataNode::Pointer new_node2;
    auto children = GetDataStorage()->GetDerivations(node);
    for (unsigned int i=0; i<children->size(); ++i)
    {
      if (children->at(i)->GetName() == "binned_centroid")
      {
        new_node = children->at(i);
        new_node->SetData(working_fib);
        new_node->SetVisibility(true);
        node->SetVisibility(false);
        mitk::RenderingManager::GetInstance()->RequestUpdateAll();
        return;
      }
    }
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

void QmitkTractometryView::AlongTractRadiomicsPreprocessing(mitk::Image::Pointer image, mitk::DataNode::Pointer node)
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

  // calculate mask
  typedef unsigned int OutPixType;
  typedef itk::Image<OutPixType, 3> OutImageType;
  itk::TractDensityImageFilter< OutImageType >::Pointer generator = itk::TractDensityImageFilter< OutImageType >::New();
  generator->SetFiberBundle(fib);
  generator->SetMode(TDI_MODE::BINARY);
  OutImageType::Pointer itkImage = OutImageType::New();
  CastToItkImage(image, itkImage);
  generator->SetInputImage(itkImage);
  generator->SetUseImageGeometry(true);
  generator->Update();
  OutImageType::Pointer count_map = generator->GetOutput();

  unsigned int num_points = m_Controls->m_SamplingPointsBox->value();
  mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();
  working_fib->ResampleToNumPoints(num_points);
  vtkSmartPointer< vtkPolyData > polydata = working_fib->GetFiberPolyData();


  itk::ImageRegionIterator< OutImageType > it(count_map, count_map->GetLargestPossibleRegion());
  while( !it.IsAtEnd() )
  {
    if (it.Get()>0)
    {
      std::vector<unsigned int> seg_vote; seg_vote.resize(num_points, 0);
      typename OutImageType::PointType image_point;
      count_map->TransformIndexToPhysicalPoint(it.GetIndex(), image_point);

      for (unsigned int i=0; i<working_fib->GetNumFibers(); ++i)
      {
        vtkCell* cell = polydata->GetCell(i);
        auto numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        bool flip = false;
        if (i>0)
          flip = Flip(polydata, i);
        else if (m_ReferencePolyData!=nullptr)
          flip = Flip(polydata, 0, m_ReferencePolyData);

        float local_d = 99999999;
        int local_closest_seg = -1;

        for (int j=0; j<numPoints; j++)
        {
          double* p;
          int segment_id = -1;
          if (flip)
          {
            segment_id = numPoints - j - 1;
            p = points->GetPoint(segment_id);
          }
          else
          {
            p = points->GetPoint(j);
            segment_id = j;
          }

          float d = std::fabs( (p[0]-image_point[0]) ) + std::fabs( (p[1]-image_point[1]) ) + std::fabs( (p[2]-image_point[2]) );
          if (d<local_d)
          {
            local_d = d;
            local_closest_seg = segment_id;
          }
        }

        seg_vote[local_closest_seg] += 1;
      }

      unsigned char final_seg_id = -1;
      unsigned int max_count = 0;
      for (unsigned int i=0; i<seg_vote.size(); ++i)
      {
        if (seg_vote.at(i)>=max_count)
        {
          final_seg_id = i;
          max_count = seg_vote.at(i);
        }
      }

      it.Set(final_seg_id + 1);
    }
    ++it;
  }


  mitk::Image::Pointer seg_img = mitk::Image::New();
  seg_img->InitializeByItk(count_map.GetPointer());
  seg_img->SetVolume(count_map->GetBufferPointer());

  mitk::DataNode::Pointer new_node = mitk::DataNode::New();
  new_node->SetData(seg_img);
  new_node->SetName("segment image");
  new_node->SetVisibility(true);
  node->SetVisibility(false);
  GetDataStorage()->Add(new_node, node);
}

void QmitkTractometryView::StartTractometry()
{
  m_ReferencePolyData = nullptr;

  double color[3] = {0,0,0};
  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(m_Controls->m_ImageBox->GetSelectedNode()->GetData());

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
    case 1:
    {
      mitkPixelTypeMultiplex4( StaticResamplingTractometry, image->GetPixelType(), image, node, data, clipboardString );
      break;
    }
    case 2:
    {
      AlongTractRadiomicsPreprocessing(image, node);
      return;
    }
    default:
    {
      mitkPixelTypeMultiplex4( NearestCentroidPointTractometry, image->GetPixelType(), image, node, data, clipboardString );
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
