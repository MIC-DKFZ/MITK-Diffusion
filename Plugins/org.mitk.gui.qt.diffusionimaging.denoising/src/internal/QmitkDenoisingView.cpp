/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkDenoisingView.h"

#include <mitkDiffusionImageDenoising.h>
#include <mitkNodePredicateIsDWI.h>

#include <QMessageBox>
#include <QtConcurrentRun>

const std::string QmitkDenoisingView::VIEW_ID = "org.mitk.views.denoisingview";

QmitkDenoisingView::QmitkDenoisingView()
  : QmitkAbstractView()
{
}

QmitkDenoisingView::~QmitkDenoisingView()
{
  // The worker must not outlive the view it reports to.
  m_Watcher.waitForFinished();
}

void QmitkDenoisingView::CreateQtPartControl( QWidget *parent )
{
  if (m_Controls)
    return;

  m_Controls = std::make_unique<Ui::QmitkDenoisingViewControls>();
  m_Controls->setupUi(parent);

  connect(m_Controls->m_ApplyButton, &QAbstractButton::clicked, this, &QmitkDenoisingView::StartDenoising);
  connect(m_Controls->m_SelectFilterComboBox, &QComboBox::activated, this, &QmitkDenoisingView::UpdateGui);
  connect(m_Controls->m_InputImageBox, &QComboBox::currentIndexChanged, this, &QmitkDenoisingView::UpdateGui);
  connect(&m_Watcher, &QFutureWatcher<Result>::finished, this, &QmitkDenoisingView::OnDenoisingFinished);

  m_Controls->m_InputImageBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_InputImageBox->SetPredicate(mitk::NodePredicateIsDWI::New());

  this->UpdateGui();
}

void QmitkDenoisingView::SetFocus()
{
  m_Controls->m_SelectFilterComboBox->setFocus();
}

void QmitkDenoisingView::StartDenoising()
{
  if (m_Watcher.isRunning())
    return;

  auto node = m_Controls->m_InputImageBox->GetSelectedNode();
  if (node.IsNull())
    return;
  mitk::Image::ConstPointer inputImage = dynamic_cast<mitk::Image*>(node->GetData());

  mitk::DiffusionImageDenoising::Parameters parameters;
  parameters.tvLambda = m_Controls->m_TvLambdaBox->value();
  parameters.tvIterations = m_Controls->m_TvIterationsBox->value();
  parameters.gaussVariance = m_Controls->m_GaussVarianceBox->value();
  parameters.nlmSearchRadius = m_Controls->m_NlmSearchRadiusBox->value();
  parameters.nlmNumberOfPatches = m_Controls->m_NlmNumPatchesBox->value();
  parameters.nlmPatchRadius = m_Controls->m_NlmPatchSizeBox->value();
  parameters.nlmIterations = m_Controls->m_NlmIterationsBox->value();

  QString method;
  switch (m_Controls->m_SelectFilterComboBox->currentIndex())
  {
  case 0:
    parameters.method = mitk::DiffusionImageDenoising::Method::TotalVariation;
    method = "TotalVariation";
    break;
  case 1:
    parameters.method = mitk::DiffusionImageDenoising::Method::Gaussian;
    method = "Gauss";
    break;
  case 2:
    parameters.method = mitk::DiffusionImageDenoising::Method::NonLocalMeans;
    method = "NLM";
    break;
  default:
    return;
  }

  m_ImageNode = node;
  m_ResultName = QString::fromStdString(node->GetName()) + "_" + method;
  m_Controls->m_ApplyButton->setEnabled(false);
  m_Controls->m_ApplyButton->setText("Denoising ...");

  // The worker only reads the input; exceptions are passed back as text.
  m_Watcher.setFuture(QtConcurrent::run([inputImage, parameters]()
  {
    Result result;
    try
    {
      result.image = mitk::DiffusionImageDenoising::Denoise(inputImage, parameters);
    }
    catch (const std::exception& e)
    {
      result.error = e.what();
    }
    return result;
  }));
}

void QmitkDenoisingView::OnDenoisingFinished()
{
  // Taken, not copied: the watcher's future would otherwise keep the whole image alive until the
  // next run, even after its node is removed.
  const Result result = m_Watcher.future().takeResult();
  m_Controls->m_ApplyButton->setText("Start");
  this->UpdateGui();

  if (result.image.IsNull())
  {
    QMessageBox::warning(nullptr, "Image not processed", QString("Denoising failed: ") + QString::fromStdString(result.error));
    return;
  }

  auto node = mitk::DataNode::New();
  node->SetName(m_ResultName.toStdString());
  node->SetData(result.image);
  // The input may have been removed while the denoising ran.
  auto dataStorage = this->GetDataStorage();
  dataStorage->Add(node, dataStorage->Exists(m_ImageNode) ? m_ImageNode.GetPointer() : nullptr);
  m_ImageNode = nullptr;
}

void QmitkDenoisingView::UpdateGui()
{
  const int method = m_Controls->m_SelectFilterComboBox->currentIndex();
  m_Controls->m_TvFrame->setVisible(method == 0);
  m_Controls->m_GaussFrame->setVisible(method == 1);
  m_Controls->m_NlmFrame->setVisible(method == 2);

  m_Controls->m_ApplyButton->setEnabled(!m_Watcher.isRunning() && m_Controls->m_InputImageBox->GetSelectedNode().IsNotNull());
}
