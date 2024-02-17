/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkTensorReconstructionView.h"
#include "mitkDiffusionImagingConfigure.h"

// qt includes
#include <QMessageBox>
#include <QImage>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QGraphicsLinearLayout>


// itk includes
#include "itkTimeProbe.h"
//#include "itkTensor.h"

// mitk includes
#include "mitkProgressBar.h"
#include "mitkStatusBar.h"

#include "mitkNodePredicateDataType.h"
#include "QmitkDataStorageComboBox.h"

#include "mitkTeemDiffusionTensor3DReconstructionImageFilter.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkTensorImageToDiffusionImageFilter.h"
#include "itkPointShell.h"
#include "itkVector.h"
#include "itkB0ImageExtractionImageFilter.h"
#include "itkTensorReconstructionWithEigenvalueCorrectionFilter.h"

#include "mitkImageCast.h"
#include "mitkImageAccessByItk.h"
#include <itkBinaryThresholdImageFilter.h>
#include <mitkImageVtkMapper2D.h>

#include "mitkProperties.h"
#include "mitkDataNodeObject.h"
#include "mitkOdfNormalizationMethodProperty.h"
#include "mitkOdfScaleByProperty.h"
#include "mitkLookupTableProperty.h"
#include "mitkLookupTable.h"
#include "mitkImageStatisticsHolder.h"
#include <mitkITKImageImport.h>

#include <itkTensorImageToOdfImageFilter.h>
#include <berryIWorkbenchWindow.h>
#include <berryISelectionService.h>
#include <mitkImageVtkMapper2D.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkNodePredicateDataType.h>

const std::string QmitkTensorReconstructionView::VIEW_ID = "org.mitk.views.tensorreconstruction";

typedef mitk::TensorImage::ScalarPixelType          TTensorPixelType;
typedef mitk::TensorImage::PixelType                TensorPixelType;
typedef mitk::TensorImage::ItkTensorImageType       TensorImageType;

using namespace berry;

QmitkTensorReconstructionView::QmitkTensorReconstructionView()
  : QmitkAbstractView(),
    m_Controls(nullptr)
{

}

QmitkTensorReconstructionView::~QmitkTensorReconstructionView()
{

}

void QmitkTensorReconstructionView::CreateQtPartControl(QWidget *parent)
{
  if (!m_Controls)
  {
    // create GUI widgets
    m_Controls = new Ui::QmitkTensorReconstructionViewControls;
    m_Controls->setupUi(parent);
    this->CreateConnections();

    Advanced1CheckboxClicked();
  }
}

void QmitkTensorReconstructionView::SetFocus()
{
  m_Controls->m_Advanced1->setFocus();
}

void QmitkTensorReconstructionView::CreateConnections()
{
  if ( m_Controls )
  {
    connect( (QObject*)(m_Controls->m_StartReconstruction), SIGNAL(clicked()), this, SLOT(Reconstruct()) );
    connect( (QObject*)(m_Controls->m_Advanced1), SIGNAL(clicked()), this, SLOT(Advanced1CheckboxClicked()) );
    connect( (QObject*)(m_Controls->m_TensorsToDWIButton), SIGNAL(clicked()), this, SLOT(TensorsToDWI()) );
    connect( (QObject*)(m_Controls->m_TensorsToOdfButton), SIGNAL(clicked()), this, SLOT(TensorsToOdf()) );
    connect( (QObject*)(m_Controls->m_TensorReconstructionThreshold), SIGNAL(valueChanged(int)), this, SLOT(PreviewThreshold(int)) );

    m_Controls->m_DwiBox->SetDataStorage(this->GetDataStorage());
    mitk::NodePredicateIsDWI::Pointer isDwi = mitk::NodePredicateIsDWI::New();
    m_Controls->m_DwiBox->SetPredicate( isDwi );
    connect( (QObject*)(m_Controls->m_DwiBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));

    m_Controls->m_OdfBox->SetDataStorage(this->GetDataStorage());
    mitk::TNodePredicateDataType<mitk::TensorImage>::Pointer isDti = mitk::TNodePredicateDataType<mitk::TensorImage>::New();
    m_Controls->m_OdfBox->SetPredicate( isDti );
    connect( (QObject*)(m_Controls->m_OdfBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));

    m_Controls->m_DtiBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_DtiBox->SetPredicate( isDti );
    connect( (QObject*)(m_Controls->m_DtiBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
  }
}

void QmitkTensorReconstructionView::Advanced1CheckboxClicked()
{
  bool check = m_Controls->
               m_Advanced1->isChecked();

  m_Controls->frame->setVisible(check);
}

void QmitkTensorReconstructionView::Activated()
{
}

void QmitkTensorReconstructionView::Deactivated()
{

  // Get all current nodes

  mitk::DataStorage::SetOfObjects::ConstPointer objects =  this->GetDataStorage()->GetAll();
  mitk::DataStorage::SetOfObjects::const_iterator itemiter( objects->begin() );
  mitk::DataStorage::SetOfObjects::const_iterator itemiterend( objects->end() );
  while ( itemiter != itemiterend ) // for all items
  {
    mitk::DataNode::Pointer node = *itemiter;
    if (node.IsNull())
      continue;

    // only look at interesting types
    if( mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dynamic_cast<mitk::Image*>(node->GetData())))
    {
      if (this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", *itemiter))
      {
        node = this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", *itemiter);
        this->GetDataStorage()->Remove(node);
      }
    }
    itemiter++;
  }


  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkTensorReconstructionView::Visible()
{
}

void QmitkTensorReconstructionView::Hidden()
{
}

void QmitkTensorReconstructionView::Reconstruct()
{
  int method = m_Controls->m_ReconctructionMethodBox->currentIndex();

  switch (method)
  {
  case 0:
    ItkTensorReconstruction();
    break;
  case 1:
    TensorReconstructionWithCorr();
    break;
  default:
    ItkTensorReconstruction();
  }
}

void QmitkTensorReconstructionView::TensorReconstructionWithCorr()
{
  try
  {
    if ( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() ) // for all items
    {
      mitk::Image* vols = static_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());
      std::string nodename = m_Controls->m_DwiBox->GetSelectedNode()->GetName();

      typedef itk::TensorReconstructionWithEigenvalueCorrectionFilter< DiffusionPixelType, TTensorPixelType > ReconstructionFilter;

      float b0Threshold = m_Controls->m_TensorReconstructionThreshold->value();
      mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer gradientContainerCopy = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
      for(auto it = mitk::DiffusionPropertyHelper::GetGradientContainer(vols)->Begin(); it != mitk::DiffusionPropertyHelper::GetGradientContainer(vols)->End(); it++)
        gradientContainerCopy->push_back(it.Value());

      ITKDiffusionImageType::Pointer itkVectorImagePointer = ITKDiffusionImageType::New();
      mitk::CastToItkImage(vols, itkVectorImagePointer);

      ReconstructionFilter::Pointer reconFilter = ReconstructionFilter::New();
      reconFilter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(vols));
      reconFilter->SetGradientImage( mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer(gradientContainerCopy), itkVectorImagePointer);
      reconFilter->SetB0Threshold(b0Threshold);
      reconFilter->Update();

      typedef mitk::TensorImage::ItkTensorImageType TensorImageType;
      TensorImageType::Pointer outputTensorImg = reconFilter->GetOutput();

      typedef itk::ImageRegionIterator<TensorImageType> TensorImageIteratorType;
      TensorImageIteratorType tensorIt(outputTensorImg, outputTensorImg->GetRequestedRegion());
      tensorIt.GoToBegin();

      int negatives = 0;
      while(!tensorIt.IsAtEnd())
      {
        typedef mitk::TensorImage::PixelType TensorType;
        TensorType tensor = tensorIt.Get();

        TensorType::EigenValuesArrayType ev;
        tensor.ComputeEigenValues(ev);

        for(unsigned int i=0; i<ev.Size(); i++)
        {
          if(ev[i] < 0.0)
          {
            tensor.Fill(0.0);
            tensorIt.Set(tensor);
            negatives++;
            break;
          }
        }
        ++tensorIt;
      }
      MITK_INFO << negatives << " tensors with negative eigenvalues" << std::endl;

      mitk::TensorImage::Pointer image = mitk::TensorImage::New();
      image->InitializeByItk( outputTensorImg.GetPointer() );
      image->SetVolume( outputTensorImg->GetBufferPointer() );
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(vols, image);
      mitk::DataNode::Pointer node=mitk::DataNode::New();
      node->SetData( image );
      node->SetName(nodename+"_EigenvalueCorrected_DT");
      GetDataStorage()->Add(node, m_Controls->m_DwiBox->GetSelectedNode());
    }

    this->GetRenderWindowPart()->RequestUpdate();
  }
  catch (itk::ExceptionObject &ex)
  {
    MITK_INFO << ex ;
    QMessageBox::information(0, "Reconstruction not possible:", ex.GetDescription());
  }
}

void QmitkTensorReconstructionView::ItkTensorReconstruction()
{
  try
  {
    if ( m_Controls->m_DwiBox->GetSelectedNode().IsNotNull() ) // for all items
    {
      mitk::Image* vols = static_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());

      std::string nodename = m_Controls->m_DwiBox->GetSelectedNode()->GetName();


      typedef itk::DiffusionTensor3DReconstructionImageFilter<DiffusionPixelType, DiffusionPixelType, TTensorPixelType > TensorReconstructionImageFilterType;
      TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = TensorReconstructionImageFilterType::New();

      GradientDirectionContainerType::Pointer gradientContainerCopy = GradientDirectionContainerType::New();
      for(GradientDirectionContainerType::ConstIterator it = mitk::DiffusionPropertyHelper::GetGradientContainer(vols)->Begin();
          it != mitk::DiffusionPropertyHelper::GetGradientContainer(vols)->End(); it++)
      {
        gradientContainerCopy->push_back(it.Value());
      }

      ITKDiffusionImageType::Pointer itkVectorImagePointer = ITKDiffusionImageType::New();
      mitk::CastToItkImage(vols, itkVectorImagePointer);

      tensorReconstructionFilter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(vols));
      tensorReconstructionFilter->SetGradientImage( gradientContainerCopy, itkVectorImagePointer );
      tensorReconstructionFilter->SetThreshold( m_Controls->m_TensorReconstructionThreshold->value() );
      tensorReconstructionFilter->Update();

      // TENSORS TO DATATREE
      mitk::TensorImage::Pointer image = mitk::TensorImage::New();

      typedef mitk::TensorImage::ItkTensorImageType TensorImageType;
      TensorImageType::Pointer tensorImage;
      tensorImage = tensorReconstructionFilter->GetOutput();

      // Check the tensor for negative eigenvalues
      if(m_Controls->m_CheckNegativeEigenvalues->isChecked())
      {
        typedef itk::ImageRegionIterator<TensorImageType> TensorImageIteratorType;
        TensorImageIteratorType tensorIt(tensorImage, tensorImage->GetRequestedRegion());
        tensorIt.GoToBegin();

        while(!tensorIt.IsAtEnd())
        {

          typedef mitk::TensorImage::PixelType TensorType;
          //typedef itk::Tensor<TTensorPixelType, 3> TensorType2;

          TensorType tensor = tensorIt.Get();

          TensorType::EigenValuesArrayType ev;
          tensor.ComputeEigenValues(ev);
          for(unsigned int i=0; i<ev.Size(); i++)
          {
            if(ev[i] < 0.0)
            {
              tensor.Fill(0.0);
              tensorIt.Set(tensor);
              break;
            }
          }
          ++tensorIt;
        }
      }

      tensorImage->SetDirection( itkVectorImagePointer->GetDirection() );
      image->InitializeByItk( tensorImage.GetPointer() );
      image->SetVolume( tensorReconstructionFilter->GetOutput()->GetBufferPointer() );
      mitk::DataNode::Pointer node=mitk::DataNode::New();
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(vols, image);
      node->SetData( image );
      node->SetName(nodename+"_LinearLeastSquares_DT");
      GetDataStorage()->Add(node, m_Controls->m_DwiBox->GetSelectedNode());
    }

    this->GetRenderWindowPart()->RequestUpdate();
  }
  catch (itk::ExceptionObject &ex)
  {
    MITK_INFO << ex ;
    QMessageBox::information(0, "Reconstruction not possible:", ex.GetDescription());
    return;
  }
}

void QmitkTensorReconstructionView::TensorsToDWI()
{
  DoTensorsToDWI();
}

void QmitkTensorReconstructionView::TensorsToOdf()
{
  if (m_Controls->m_OdfBox->GetSelectedNode().IsNotNull())
  {
    mitk::DataNode::Pointer tensorImageNode = m_Controls->m_OdfBox->GetSelectedNode();

    typedef mitk::TensorImage::ScalarPixelType    TTensorPixelType;
    typedef mitk::TensorImage::ItkTensorImageType TensorImageType;

    TensorImageType::Pointer itkvol = TensorImageType::New();
    mitk::CastToItkImage(dynamic_cast<mitk::TensorImage*>(tensorImageNode->GetData()), itkvol);

    typedef itk::TensorImageToOdfImageFilter< TTensorPixelType, TTensorPixelType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( itkvol );
    filter->Update();

    typedef itk::Vector<TTensorPixelType,ODF_SAMPLING_SIZE>  OutputPixelType;
    typedef itk::Image<OutputPixelType,3>                OutputImageType;

    mitk::OdfImage::Pointer image = mitk::OdfImage::New();
    OutputImageType::Pointer outimg = filter->GetOutput();
    image->InitializeByItk( outimg.GetPointer() );
    image->SetVolume( outimg->GetBufferPointer() );
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    mitk::DiffusionPropertyHelper::CopyDICOMProperties(tensorImageNode->GetData(), image);
    node->SetData( image );
    node->SetName(tensorImageNode->GetName()+"_Odf");
    GetDataStorage()->Add(node, tensorImageNode);
  }
}

void QmitkTensorReconstructionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& )
{
  UpdateGui();
}

void QmitkTensorReconstructionView::UpdateGui()
{
  m_Controls->m_StartReconstruction->setEnabled(m_Controls->m_DwiBox->GetSelectedNode().IsNotNull());

  m_Controls->m_TensorsToDWIButton->setEnabled(m_Controls->m_DtiBox->GetSelectedNode().IsNotNull());
  m_Controls->m_TensorsToOdfButton->setEnabled(m_Controls->m_OdfBox->GetSelectedNode().IsNotNull());
}

template<int ndirs>
itk::VectorContainer<unsigned int, vnl_vector_fixed<double, 3> >::Pointer
QmitkTensorReconstructionView::MakeGradientList()
{
  itk::VectorContainer<unsigned int, vnl_vector_fixed<double,3> >::Pointer retval =
      itk::VectorContainer<unsigned int, vnl_vector_fixed<double,3> >::New();
  vnl_matrix_fixed<double, 3, ndirs>* U =
      itk::PointShell<ndirs, vnl_matrix_fixed<double, 3, ndirs> >::DistributePointShell();

  for(int i=0; i<ndirs;i++)
  {
    vnl_vector_fixed<double,3> v;
    v[0] = U->get(0,i); v[1] = U->get(1,i); v[2] = U->get(2,i);
    retval->push_back(v);
  }
  // Add 0 vector for B0
  vnl_vector_fixed<double,3> v;
  v.fill(0.0);
  retval->push_back(v);

  return retval;
}

void QmitkTensorReconstructionView::DoTensorsToDWI()
{
  try
  {
    if (m_Controls->m_DtiBox->GetSelectedNode().IsNotNull())
    {
      std::string nodename = m_Controls->m_DtiBox->GetSelectedNode()->GetName();

      mitk::TensorImage* vol = static_cast<mitk::TensorImage*>(m_Controls->m_DtiBox->GetSelectedNode()->GetData());

      typedef mitk::TensorImage::ScalarPixelType      TTensorPixelType;
      typedef mitk::TensorImage::ItkTensorImageType   TensorImageType;


      TensorImageType::Pointer itkvol = TensorImageType::New();
      mitk::CastToItkImage(vol, itkvol);

      typedef itk::TensorImageToDiffusionImageFilter<
          TTensorPixelType, DiffusionPixelType > FilterType;

      FilterType::GradientListPointerType gradientList;

      switch(m_Controls->m_TensorsToDWINumDirsSelect->currentIndex())
      {
      case 0:
        gradientList = MakeGradientList<12>();
        break;
      case 1:
        gradientList = MakeGradientList<42>();
        break;
      case 2:
        gradientList = MakeGradientList<92>();
        break;
      case 3:
        gradientList = MakeGradientList<162>();
        break;
      case 4:
        gradientList = MakeGradientList<252>();
        break;
      case 5:
        gradientList = MakeGradientList<362>();
        break;
      case 6:
        gradientList = MakeGradientList<492>();
        break;
      case 7:
        gradientList = MakeGradientList<642>();
        break;
      case 8:
        gradientList = MakeGradientList<812>();
        break;
      case 9:
        gradientList = MakeGradientList<1002>();
        break;
      default:
        gradientList = MakeGradientList<92>();

      }
      double bVal = m_Controls->m_TensorsToDWIBValueEdit->text().toDouble();

      FilterType::Pointer filter = FilterType::New();
      filter->SetInput( itkvol );
      filter->SetBValue(bVal);
      filter->SetGradientList(gradientList);
      filter->Update();

      mitk::Image::Pointer image = mitk::GrabItkImageMemory( filter->GetOutput() );

      mitk::DiffusionPropertyHelper::SetGradientContainer(image, gradientList);
      mitk::DiffusionPropertyHelper::SetReferenceBValue(image, bVal);
      mitk::DiffusionPropertyHelper::InitializeImage( image );

      mitk::DataNode::Pointer node=mitk::DataNode::New();
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(vol, image);
      node->SetData( image );
      mitk::ImageVtkMapper2D::SetDefaultProperties(node);
      node->SetName(nodename+"_DWI");
      GetDataStorage()->Add(node, m_Controls->m_DtiBox->GetSelectedNode());
    }
    this->GetRenderWindowPart()->RequestUpdate();
  }
  catch (itk::ExceptionObject &ex)
  {
    MITK_INFO << ex ;
    QMessageBox::information(0, "DWI estimation failed:", ex.GetDescription());
    return ;
  }
}


void QmitkTensorReconstructionView::PreviewThreshold(int threshold)
{
  if (m_Controls->m_DwiBox->GetSelectedNode().IsNotNull())
  {
    mitk::Image* vols = static_cast<mitk::Image*>(m_Controls->m_DwiBox->GetSelectedNode()->GetData());

    ITKDiffusionImageType::Pointer itkVectorImagePointer = ITKDiffusionImageType::New();
    mitk::CastToItkImage(vols, itkVectorImagePointer);

    // Extract b0 image
    typedef itk::B0ImageExtractionImageFilter<short, short> FilterType;
    FilterType::Pointer filterB0 = FilterType::New();
    filterB0->SetInput(itkVectorImagePointer);
    filterB0->SetDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(vols));
    filterB0->Update();

    mitk::Image::Pointer mitkImage = mitk::Image::New();

    typedef itk::Image<short, 3> ImageType;
    typedef itk::Image<short, 3> SegmentationType;
    typedef itk::BinaryThresholdImageFilter<ImageType, SegmentationType> ThresholdFilterType;
    // apply threshold
    ThresholdFilterType::Pointer filterThreshold = ThresholdFilterType::New();
    filterThreshold->SetInput(filterB0->GetOutput());
    filterThreshold->SetLowerThreshold(threshold);
    filterThreshold->SetInsideValue(0);
    filterThreshold->SetOutsideValue(1); // mark cut off values red
    filterThreshold->Update();

    mitkImage->InitializeByItk( filterThreshold->GetOutput() );
    mitkImage->SetVolume( filterThreshold->GetOutput()->GetBufferPointer() );
    mitk::DataNode::Pointer node;
    if (this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", m_Controls->m_DwiBox->GetSelectedNode()))
    {
      node = this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", m_Controls->m_DwiBox->GetSelectedNode());
    }
    else
    {
      // create a new node, to show thresholded values
      node = mitk::DataNode::New();
      GetDataStorage()->Add( node, m_Controls->m_DwiBox->GetSelectedNode() );
      node->SetProperty( "name", mitk::StringProperty::New("ThresholdOverlay"));
      node->SetBoolProperty("helper object", true);
    }
    node->SetData( mitkImage );
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}
