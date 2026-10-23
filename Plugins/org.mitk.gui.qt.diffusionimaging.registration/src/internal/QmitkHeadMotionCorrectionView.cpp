/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkHeadMotionCorrectionView.h"

#include <mitkDWIHeadMotionCorrectionFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkNodePredicateIsDWI.h>

#include <QMessageBox>
#include <QtConcurrentRun>

const std::string QmitkHeadMotionCorrectionView::VIEW_ID = "org.mitk.views.headmotioncorrectionview";

QmitkHeadMotionCorrectionView::QmitkHeadMotionCorrectionView()
  : QmitkAbstractView()
{
}

QmitkHeadMotionCorrectionView::~QmitkHeadMotionCorrectionView()
{
  // The worker must not outlive the view it reports to.
  m_Watcher.waitForFinished();
}

void QmitkHeadMotionCorrectionView::CreateQtPartControl( QWidget *parent )
{
  if (m_Controls)
    return;

  m_Controls = std::make_unique<Ui::QmitkHeadMotionCorrectionViewControls>();
  m_Controls->setupUi(parent);
  connect(m_Controls->m_ImageBox, &QComboBox::currentIndexChanged, this, &QmitkHeadMotionCorrectionView::UpdateGUI);
  connect(m_Controls->m_RegistrationStartButton, &QAbstractButton::clicked, this, &QmitkHeadMotionCorrectionView::StartCorrection);
  connect(&m_Watcher, &QFutureWatcher<Result>::finished, this, &QmitkHeadMotionCorrectionView::OnCorrectionFinished);

  m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateIsDWI::New());
  this->UpdateGUI();
}

void QmitkHeadMotionCorrectionView::UpdateGUI()
{
  const bool running = m_Watcher.isRunning();
  const bool selected = m_Controls->m_ImageBox->GetSelectedNode().IsNotNull();
  m_Controls->m_RegistrationStartButton->setEnabled(!running && selected);
  if (running)
    m_Controls->m_RegistrationStartButton->setToolTip("A correction is running.");
  else
    m_Controls->m_RegistrationStartButton->setToolTip(selected ? "Start Registration" : "No Diffusion image selected.");
}

void QmitkHeadMotionCorrectionView::SetFocus()
{
  this->UpdateGUI();
  m_Controls->m_RegistrationStartButton->setFocus();
}

void QmitkHeadMotionCorrectionView::StartCorrection()
{
  if (m_Watcher.isRunning())
    return;

  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  if (node.IsNull())
    return;
  mitk::Image::ConstPointer inImage = dynamic_cast<mitk::Image*>(node->GetData());
  if (!mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(inImage))
  {
    QMessageBox::warning(nullptr, "Image not processed", "Input is not a diffusion-weighted image!");
    return;
  }

  m_ImageNode = node;
  m_Controls->m_RegistrationStartButton->setText("Correction in progress ...");

  // The worker only reads the input; exceptions are passed back as text.
  m_Watcher.setFuture(QtConcurrent::run([inImage]()
  {
    Result result;
    try
    {
      auto registerer = mitk::DWIHeadMotionCorrectionFilter::New();
      registerer->SetInput(inImage);
      registerer->Update();
      result.image = registerer->GetCorrectedImage();
      result.image->DisconnectPipeline();
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(inImage, result.image);
    }
    catch (const std::exception& e)
    {
      result.error = e.what();
    }
    return result;
  }));
  // Only now isRunning() is true, so the button is disabled.
  this->UpdateGUI();
}

void QmitkHeadMotionCorrectionView::OnCorrectionFinished()
{
  // Taken, not copied: the watcher's future would otherwise keep the whole image alive until the
  // next run, even after its node is removed.
  const Result result = m_Watcher.future().takeResult();
  m_Controls->m_RegistrationStartButton->setText("Start Head Motion Correction");
  this->UpdateGUI();

  if (result.image.IsNull())
  {
    QMessageBox::warning(nullptr, "Image not processed", QString("Head-motion correction failed: ") + QString::fromStdString(result.error));
    return;
  }

  auto correctedNode = mitk::DataNode::New();
  correctedNode->SetData(result.image);
  correctedNode->SetName(m_ImageNode->GetName() + "_Corrected");
  // The input may have been removed while the correction ran.
  auto dataStorage = this->GetDataStorage();
  dataStorage->Add(correctedNode, dataStorage->Exists(m_ImageNode) ? m_ImageNode.GetPointer() : nullptr);
  m_ImageNode = nullptr;
}
