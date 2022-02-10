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
#include "QmitkFiberQuantificationView.h"

// Qt
#include <QMessageBox>

// MITK
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateAnd.h>
#include <mitkImageCast.h>
#include <mitkPeakImage.h>
#include <mitkLabelSetImage.h>
#include <mitkDICOMSegmentationConstants.h>
#include <mitkDICOMSegmentationPropertyHelper.cpp>
#include <mitkDICOMQIPropertyHelper.h>

// ITK
#include <itkTractDensityImageFilter.h>
#include <itkTractsToRgbaImageFilter.h>
#include <itkTractsToVectorImageFilter.h>
#include <itkTractsToFiberEndingsImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkBinaryThinningImageFilter.h>

#include <mitkLexicalCast.h>

const std::string QmitkFiberQuantificationView::VIEW_ID = "org.mitk.views.fiberquantification";
using namespace mitk;

QmitkFiberQuantificationView::QmitkFiberQuantificationView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_UpsamplingFactor(5)
  , m_Visible(false)
{

}

// Destructor
QmitkFiberQuantificationView::~QmitkFiberQuantificationView()
{

}

void QmitkFiberQuantificationView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkFiberQuantificationViewControls;
    m_Controls->setupUi( parent );

    connect( m_Controls->m_ProcessFiberBundleButton, SIGNAL(clicked()), this, SLOT(ProcessSelectedBundles()) );
    connect( m_Controls->m_ExtractFiberPeaks, SIGNAL(clicked()), this, SLOT(CalculateFiberDirections()) );

    m_Controls->m_TractBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isFib = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    m_Controls->m_TractBox->SetPredicate( isFib );

    m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_ImageBox->SetZeroEntryText("--");
    mitk::TNodePredicateDataType<mitk::Image>::Pointer isImagePredicate = mitk::TNodePredicateDataType<mitk::Image>::New();
    mitk::NodePredicateDimension::Pointer is3D = mitk::NodePredicateDimension::New(3);
    m_Controls->m_ImageBox->SetPredicate( mitk::NodePredicateAnd::New(isImagePredicate, is3D) );

    connect( (QObject*)(m_Controls->m_TractBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
    connect( (QObject*)(m_Controls->m_ImageBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
  }
}

void QmitkFiberQuantificationView::Activated()
{

}

void QmitkFiberQuantificationView::Deactivated()
{

}

void QmitkFiberQuantificationView::Visible()
{
  m_Visible = true;
}

void QmitkFiberQuantificationView::Hidden()
{
  m_Visible = false;
}

void QmitkFiberQuantificationView::SetFocus()
{
  m_Controls->m_ProcessFiberBundleButton->setFocus();
}

void QmitkFiberQuantificationView::CalculateFiberDirections()
{
  typedef itk::Image<unsigned char, 3>                                            ItkUcharImgType;

  // load fiber bundle
  mitk::FiberBundle::Pointer inputTractogram = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.back()->GetData());

  itk::TractsToVectorImageFilter<float>::Pointer fOdfFilter = itk::TractsToVectorImageFilter<float>::New();
  if (m_SelectedImage.IsNotNull())
  {
    ItkUcharImgType::Pointer itkMaskImage = ItkUcharImgType::New();
    mitk::CastToItkImage<ItkUcharImgType>(m_SelectedImage, itkMaskImage);
    fOdfFilter->SetMaskImage(itkMaskImage);
  }

  // extract directions from fiber bundle
  fOdfFilter->SetFiberBundle(inputTractogram);
  fOdfFilter->SetAngularThreshold(cos(m_Controls->m_AngularThreshold->value()*itk::Math::pi/180));
  switch (m_Controls->m_FiberDirNormBox->currentIndex())
  {
  case 0:
    fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::GLOBAL_MAX);
    break;
  case 1:
    fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::SINGLE_VEC_NORM);
    break;
  case 2:
    fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::MAX_VEC_NORM);
    break;
  }
  fOdfFilter->SetOnlyUseMaskGeometry(true);
  fOdfFilter->SetSizeThreshold(m_Controls->m_PeakThreshold->value());
  fOdfFilter->SetMaxNumDirections(m_Controls->m_MaxNumDirections->value());
  fOdfFilter->Update();

  QString name = m_SelectedFB.back()->GetName().c_str();

  if (m_Controls->m_NumDirectionsBox->isChecked())
  {
    mitk::Image::Pointer mitkImage = mitk::Image::New();
    mitkImage->InitializeByItk( fOdfFilter->GetNumDirectionsImage().GetPointer() );
    mitkImage->SetVolume( fOdfFilter->GetNumDirectionsImage()->GetBufferPointer() );

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(mitkImage);
    node->SetName((name+"_NUM_DIRECTIONS").toStdString().c_str());
    GetDataStorage()->Add(node, m_SelectedFB.back());
  }

  Image::Pointer mitkImage = dynamic_cast<Image*>(PeakImage::New().GetPointer());
  mitk::CastToMitkImage(fOdfFilter->GetDirectionImage(), mitkImage);
  mitkImage->SetVolume(fOdfFilter->GetDirectionImage()->GetBufferPointer());

  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(mitkImage);
  node->SetName( (name+"_DIRECTIONS").toStdString().c_str());
  GetDataStorage()->Add(node, m_SelectedFB.back());
}

void QmitkFiberQuantificationView::UpdateGui()
{
  m_SelectedFB.clear();
  if (m_Controls->m_TractBox->GetSelectedNode().IsNotNull())
    m_SelectedFB.push_back(m_Controls->m_TractBox->GetSelectedNode());

  m_SelectedImage = nullptr;
  if (m_Controls->m_ImageBox->GetSelectedNode().IsNotNull())
    m_SelectedImage = dynamic_cast<mitk::Image*>(m_Controls->m_ImageBox->GetSelectedNode()->GetData());

  m_Controls->m_ProcessFiberBundleButton->setEnabled(!m_SelectedFB.empty());
  m_Controls->m_ExtractFiberPeaks->setEnabled(!m_SelectedFB.empty());
}

void QmitkFiberQuantificationView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& )
{
  UpdateGui();
}

void QmitkFiberQuantificationView::ProcessSelectedBundles()
{
  if ( m_SelectedFB.empty() ){
    QMessageBox::information( nullptr, "Warning", "No fibe bundle selected!");
    MITK_WARN("QmitkFiberQuantificationView") << "no fibe bundle selected";
    return;
  }

  int generationMethod = m_Controls->m_GenerationBox->currentIndex();

  for( unsigned int i=0; i<m_SelectedFB.size(); i++ )
  {
    mitk::DataNode::Pointer node = m_SelectedFB[i];
    if (node.IsNotNull() && dynamic_cast<mitk::FiberBundle*>(node->GetData()))
    {
      mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
      QString name(node->GetName().c_str());
      DataNode::Pointer newNode = nullptr;
      switch(generationMethod){
      case 0:
        newNode = GenerateTractDensityImage(fib, TDI_MODE::DENSITY, true, node->GetName());
        name += "_TDI";
        break;
      case 1:
        newNode = GenerateTractDensityImage(fib, TDI_MODE::DENSITY, false, node->GetName());
        name += "_TDI";
        break;
      case 2:
        newNode = GenerateTractDensityImage(fib, TDI_MODE::BINARY, false, node->GetName());
        name += "_envelope";
        break;
      case 3:
        newNode = GenerateColorHeatmap(fib);
        break;
      case 4:
        newNode = GenerateFiberEndingsImage(fib);
        name += "_fiber_endings";
        break;
      case 5:
        newNode = GenerateFiberEndingsPointSet(fib);
        name += "_fiber_endings";
        break;
      case 6:
        newNode = GenerateDistanceMap(fib);
        name += "_distance_map";
        break;
      case 7:
        newNode = GenerateBinarySkeleton(fib);
        name += "_skeleton";
        break;
      case 8:
        newNode = GenerateTractDensityImage(fib, TDI_MODE::VISITATION_COUNT, true, node->GetName());
        name += "_visitations";
        break;
      }
      if (newNode.IsNotNull())
      {
        newNode->SetName(name.toStdString());
        GetDataStorage()->Add(newNode);
      }
    }
  }
}

// generate pointset displaying the fiber endings
mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateFiberEndingsPointSet(mitk::FiberBundle::Pointer fib)
{
  mitk::PointSet::Pointer pointSet = mitk::PointSet::New();
  vtkSmartPointer<vtkPolyData> fiberPolyData = fib->GetFiberPolyData();

  int count = 0;
  int numFibers = fib->GetNumFibers();
  for( int i=0; i<numFibers; i++ )
  {
    vtkCell* cell = fiberPolyData->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (numPoints>0)
    {
      double* point = points->GetPoint(0);
      itk::Point<double,3> itkPoint = mitk::imv::GetItkPoint(point);
      pointSet->InsertPoint(count, itkPoint);
      count++;
    }
    if (numPoints>2)
    {
      double* point = points->GetPoint(numPoints-1);
      itk::Point<double,3> itkPoint = mitk::imv::GetItkPoint(point);
      pointSet->InsertPoint(count, itkPoint);
      count++;
    }
  }

  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData( pointSet );
  return node;
}

// generate image displaying the fiber endings
mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateFiberEndingsImage(mitk::FiberBundle::Pointer fib)
{
  typedef unsigned int OutPixType;
  typedef itk::Image<OutPixType,3> OutImageType;

  typedef itk::TractsToFiberEndingsImageFilter< OutImageType > ImageGeneratorType;
  ImageGeneratorType::Pointer generator = ImageGeneratorType::New();
  generator->SetFiberBundle(fib);
  generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
  if (m_SelectedImage.IsNotNull())
  {
    OutImageType::Pointer itkImage = OutImageType::New();
    CastToItkImage(m_SelectedImage, itkImage);
    generator->SetInputImage(itkImage);
    generator->SetUseImageGeometry(true);
  }
  generator->Update();

  // get output image
  OutImageType::Pointer outImg = generator->GetOutput();
  mitk::Image::Pointer img = mitk::Image::New();
  img->InitializeByItk(outImg.GetPointer());
  img->SetVolume(outImg->GetBufferPointer());

  // init data node
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(img);
  return node;
}

// generate rgba heatmap from fiber bundle
mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateColorHeatmap(mitk::FiberBundle::Pointer fib)
{
  typedef itk::RGBAPixel<unsigned char> OutPixType;
  typedef itk::Image<OutPixType, 3> OutImageType;
  typedef itk::TractsToRgbaImageFilter< OutImageType > ImageGeneratorType;
  ImageGeneratorType::Pointer generator = ImageGeneratorType::New();
  generator->SetFiberBundle(fib);
  generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
  if (m_SelectedImage.IsNotNull())
  {
    itk::Image<unsigned char, 3>::Pointer itkImage = itk::Image<unsigned char, 3>::New();
    CastToItkImage(m_SelectedImage, itkImage);
    generator->SetInputImage(itkImage);
    generator->SetUseImageGeometry(true);
  }
  generator->Update();

  // get output image
  typedef itk::Image<OutPixType,3> OutType;
  OutType::Pointer outImg = generator->GetOutput();
  mitk::Image::Pointer img = mitk::Image::New();
  img->InitializeByItk(outImg.GetPointer());
  img->SetVolume(outImg->GetBufferPointer());

  // init data node
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(img);
  return node;
}

mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateBinarySkeleton(mitk::FiberBundle::Pointer fib)
{
  typedef itk::Image<unsigned char, 3> UcharImageType;

  itk::TractDensityImageFilter< UcharImageType >::Pointer envelope_generator = itk::TractDensityImageFilter< UcharImageType >::New();
  envelope_generator->SetFiberBundle(fib);
  envelope_generator->SetMode(TDI_MODE::BINARY);
  envelope_generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
  if (m_SelectedImage.IsNotNull())
  {
    UcharImageType::Pointer itkImage = UcharImageType::New();
    CastToItkImage(m_SelectedImage, itkImage);
    envelope_generator->SetInputImage(itkImage);
    envelope_generator->SetUseImageGeometry(true);

  }
  envelope_generator->Update();

  itk::BinaryThinningImageFilter<UcharImageType, UcharImageType>::Pointer map_generator = itk::BinaryThinningImageFilter<UcharImageType, UcharImageType>::New();
  map_generator->SetInput(envelope_generator->GetOutput());
  map_generator->Update();

  UcharImageType::Pointer outImg = map_generator->GetOutput();
  mitk::Image::Pointer img = mitk::Image::New();
  img->InitializeByItk(outImg.GetPointer());
  img->SetVolume(outImg->GetBufferPointer());

  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(img);
  return node;
}

mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateDistanceMap(mitk::FiberBundle::Pointer fib)
{
  typedef itk::Image<unsigned char, 3> UcharImageType;
  typedef itk::Image<float, 3> FloatImageType;

  itk::TractDensityImageFilter< UcharImageType >::Pointer envelope_generator = itk::TractDensityImageFilter< UcharImageType >::New();
  envelope_generator->SetFiberBundle(fib);
  envelope_generator->SetMode(TDI_MODE::BINARY);
  envelope_generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
  if (m_SelectedImage.IsNotNull())
  {
    UcharImageType::Pointer itkImage = UcharImageType::New();
    CastToItkImage(m_SelectedImage, itkImage);
    envelope_generator->SetInputImage(itkImage);
    envelope_generator->SetUseImageGeometry(true);

  }
  envelope_generator->Update();

  itk::SignedMaurerDistanceMapImageFilter<UcharImageType, FloatImageType>::Pointer map_generator = itk::SignedMaurerDistanceMapImageFilter<UcharImageType, FloatImageType>::New();
  map_generator->SetInput(envelope_generator->GetOutput());
  map_generator->SetUseImageSpacing(true);
  map_generator->SetSquaredDistance(false);
  map_generator->SetInsideIsPositive(true);
  map_generator->Update();

  FloatImageType::Pointer outImg = map_generator->GetOutput();
  mitk::Image::Pointer img = mitk::Image::New();
  img->InitializeByItk(outImg.GetPointer());
  img->SetVolume(outImg->GetBufferPointer());

  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(img);
  return node;
}

// generate tract density image from fiber bundle
mitk::DataNode::Pointer QmitkFiberQuantificationView::GenerateTractDensityImage(mitk::FiberBundle::Pointer fib, TDI_MODE mode, bool absolute, std::string name)
{
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  if (mode==TDI_MODE::BINARY)
  {
    typedef unsigned char OutPixType;
    typedef itk::Image<OutPixType, 3> OutImageType;

    itk::TractDensityImageFilter< OutImageType >::Pointer generator = itk::TractDensityImageFilter< OutImageType >::New();
    generator->SetFiberBundle(fib);
    generator->SetMode(mode);
    generator->SetOutputAbsoluteValues(absolute);
    generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
    if (m_SelectedImage.IsNotNull())
    {
      OutImageType::Pointer itkImage = OutImageType::New();
      CastToItkImage(m_SelectedImage, itkImage);
      generator->SetInputImage(itkImage);
      generator->SetUseImageGeometry(true);

    }
    generator->Update();

    // get output image
    typedef itk::Image<OutPixType,3> OutType;
    OutType::Pointer outImg = generator->GetOutput();
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg.GetPointer());
    img->SetVolume(outImg->GetBufferPointer());


    if (m_SelectedImage.IsNotNull())
    {
      mitk::LabelSetImage::Pointer multilabelImage = mitk::LabelSetImage::New();
      multilabelImage->InitializeByLabeledImage(img);
      multilabelImage->GetActiveLabelSet()->SetActiveLabel(1);
      mitk::Label::Pointer label = multilabelImage->GetActiveLabel();
      label->SetName("Tractogram");

      // Add Segmented Property Category Code Sequence tags (0062, 0003): Sequence defining the general category of this
      // segment.
      // (0008,0100) Code Value
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_CATEGORY_CODE_VALUE_PATH()).c_str(),
            TemporoSpatialStringProperty::New("T-D000A"));

      // (0008,0102) Coding Scheme Designator
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_CATEGORY_CODE_SCHEME_PATH()).c_str(),
            TemporoSpatialStringProperty::New("SRT"));

      // (0008,0104) Code Meaning
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_CATEGORY_CODE_MEANING_PATH()).c_str(),
            TemporoSpatialStringProperty::New("Anatomical Structure"));
      //------------------------------------------------------------
      // Add Segmented Property Type Code Sequence (0062, 000F): Sequence defining the specific property type of this
      // segment.
      // (0008,0100) Code Value
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_TYPE_CODE_VALUE_PATH()).c_str(),
            TemporoSpatialStringProperty::New("DUMMY"));

      // (0008,0102) Coding Scheme Designator
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_TYPE_CODE_SCHEME_PATH()).c_str(),
            TemporoSpatialStringProperty::New("SRT"));

      // (0008,0104) Code Meaning
      label->SetProperty(
            DICOMTagPathToPropertyName(DICOMSegmentationConstants::SEGMENT_TYPE_CODE_MEANING_PATH()).c_str(),
            TemporoSpatialStringProperty::New(name));

      //Error: is undeclared//      mitk::DICOMQIPropertyHandler::DeriveDICOMSourceProperties(m_SelectedImage, multilabelImage);

      // init data node
      node->SetData(multilabelImage);
    }
    else
    {
      // init data node
      node->SetData(img);
    }

  }
  else
  {
    typedef float OutPixType;
    typedef itk::Image<OutPixType, 3> OutImageType;

    itk::TractDensityImageFilter< OutImageType >::Pointer generator = itk::TractDensityImageFilter< OutImageType >::New();
    generator->SetFiberBundle(fib);
    generator->SetMode(mode);
    generator->SetOutputAbsoluteValues(absolute);
    generator->SetUpsamplingFactor(m_Controls->m_UpsamplingSpinBox->value());
    if (m_SelectedImage.IsNotNull())
    {
      OutImageType::Pointer itkImage = OutImageType::New();
      CastToItkImage(m_SelectedImage, itkImage);
      generator->SetInputImage(itkImage);
      generator->SetUseImageGeometry(true);

    }
    //generator->SetDoFiberResampling(false);
    generator->Update();

    // get output image
    typedef itk::Image<OutPixType,3> OutType;
    OutType::Pointer outImg = generator->GetOutput();
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg.GetPointer());
    img->SetVolume(outImg->GetBufferPointer());

    // init data node
    node->SetData(img);
  }
  return node;
}
