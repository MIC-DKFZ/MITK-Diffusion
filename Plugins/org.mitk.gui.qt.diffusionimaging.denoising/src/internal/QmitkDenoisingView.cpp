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
#include "QmitkDenoisingView.h"
#include <mitkImageCast.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>

#include <mitkDiffusionPropertyHelper.h>

#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <itkMaskImageFilter.h>
#include <itkSpatialNeighborSubsampler.h>
#include <QMessageBox>
#include <mitkNodePredicateIsDWI.h>
#include <mitkDiffusionImageDenoising.h>

const std::string QmitkDenoisingView::VIEW_ID = "org.mitk.views.denoisingview";

QmitkDenoisingView::QmitkDenoisingView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_ImageNode(nullptr)
{
}

QmitkDenoisingView::~QmitkDenoisingView()
{

}

void QmitkDenoisingView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkDenoisingViewControls;
    m_Controls->setupUi( parent );

    connect( (QObject*)(m_Controls->m_ApplyButton), SIGNAL(clicked()), this, SLOT(StartDenoising()));
    connect( (QObject*)(m_Controls->m_SelectFilterComboBox), SIGNAL(activated(int)), this, SLOT(UpdateGui()));
    connect( (QObject*)(m_Controls->m_InputImageBox), SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));


    m_Controls->m_InputImageBox->SetDataStorage(this->GetDataStorage());
    mitk::NodePredicateIsDWI::Pointer isImagePredicate = mitk::NodePredicateIsDWI::New();
    m_Controls->m_InputImageBox->SetPredicate( isImagePredicate );

    UpdateGui();
  }
}

void QmitkDenoisingView::SetFocus()
{
  m_Controls->m_SelectFilterComboBox->setFocus();
}

void QmitkDenoisingView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& )
{

}

void QmitkDenoisingView::StartDenoising()
{
  m_ImageNode = m_Controls->m_InputImageBox->GetSelectedNode();
  if (m_ImageNode.IsNull())
    return;
  mitk::Image::Pointer input_image = dynamic_cast<mitk::Image*>(m_ImageNode->GetData());

  if (!mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(input_image))
  {
    QMessageBox::warning(nullptr, "Image not processed", "Input is not a diffusion-weighted image!");
    return;
  }

  mitk::DiffusionImageDenoising::Parameters parameters;
  parameters.tvLambda = m_Controls->m_TvLambdaBox->value();
  parameters.tvIterations = m_Controls->m_TvIterationsBox->value();
  parameters.gaussVariance = m_Controls->m_GaussVarianceBox->value();
  parameters.nlmSearchRadius = m_Controls->m_NlmSearchRadiusBox->value();
  parameters.nlmNumberOfPatches = m_Controls->m_NlmNumPatchesBox->value();
  parameters.nlmPatchRadius = m_Controls->m_NlmPatchSizeBox->value();
  parameters.nlmIterations = m_Controls->m_NlmIterationsBox->value();

  QString name;
  switch (m_Controls->m_SelectFilterComboBox->currentIndex())
  {
  case 0:
    parameters.method = mitk::DiffusionImageDenoising::Method::TotalVariation;
    name = "TotalVariation";
    break;
  case 1:
    parameters.method = mitk::DiffusionImageDenoising::Method::Gaussian;
    name = "Gauss";
    break;
  case 2:
    parameters.method = mitk::DiffusionImageDenoising::Method::NonLocalMeans;
    name = "NLM";
    break;
  default:
    return;
  }

  mitk::Image::Pointer denoised_image;
  try
  {
    denoised_image = mitk::DiffusionImageDenoising::Denoise(input_image, parameters);
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Image not processed", QString("Denoising failed: ") + e.what());
    return;
  }

  mitk::DataNode::Pointer denoised_image_node = mitk::DataNode::New();
  denoised_image_node->SetName( name.toStdString() );
  denoised_image_node->SetData( denoised_image );
  GetDataStorage()->Add(denoised_image_node, m_ImageNode);
}

void QmitkDenoisingView::UpdateGui()
{
  m_Controls->m_ApplyButton->setEnabled(false);
  m_Controls->m_TvFrame->setVisible(false);
  m_Controls->m_NlmFrame->setVisible(false);
  m_Controls->m_GaussFrame->setVisible(false);

  switch (m_Controls->m_SelectFilterComboBox->currentIndex())
  {
  case 0:
  {
    m_Controls->m_TvFrame->setVisible(true);
    break;
  }
  case 1:
  {
    m_Controls->m_GaussFrame->setVisible(true);
    break;
  }
  case 2:
  {
    m_Controls->m_NlmFrame->setVisible(true);
    break;
  }
  }

  if (m_Controls->m_InputImageBox->GetSelectedNode().IsNotNull())
    m_Controls->m_ApplyButton->setEnabled(true);
}

