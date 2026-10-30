/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkSimpleRegistrationView.h"

#include <QmitkDiffusionViewControlsGuard.h>

#include <mitkNodePredicateDataType.h>
#include <mitkDiffusionPropertyHelper.h>
// MatchPoint's ITKInitializedMultiResImageRegistrationAlgorithm::doSetProperty() assigns a bool
// that unwrapMetaProperty() may leave unset; GCC's optimizer reports that in Release builds
// despite -isystem. MatchPoint is a MITK external, so the warning is silenced for its headers.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <mitkMultiModalRigidDefaultRegistrationAlgorithm.h>
#include <mitkMultiModalAffineDefaultRegistrationAlgorithm.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#include <mitkFiberBundle.h>
#include <mitkRegistrationHelper.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionImageMapping.h>

#include <QMessageBox>
#include <QThreadPool>

const std::string QmitkSimpleRegistrationView::VIEW_ID = "org.mitk.views.simpleregistrationview";

QmitkSimpleRegistrationView::QmitkSimpleRegistrationView()
  : QmitkAbstractView()
{
}

QmitkSimpleRegistrationView::~QmitkSimpleRegistrationView() = default;

void QmitkSimpleRegistrationView::StartRegistration()
{
  // The result slot uses the moving image and type of the running job.
  if (m_RegistrationRunning)
    return;

  // Of a DWI, the volume chosen in the channel box is registered. The inputs are prepared before
  // the job is allocated, so a failure neither leaks the job nor leaves this slot with an exception.
  auto movingNode = m_Controls->m_MovingImageBox->GetSelectedNode();
  auto fixedNode = m_Controls->m_FixedImageBox->GetSelectedNode();
  if (movingNode.IsNull() || fixedNode.IsNull())
    return;
  mitk::Image::Pointer movingData, targetData;
  mitk::NodeUIDType movingUID, targetUID;
  try
  {
    mitk::Image::Pointer movingImage = dynamic_cast<mitk::Image*>(movingNode->GetData());
    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(movingImage))
      movingData = mitk::DiffusionImageEditing::ExtractVolume(movingImage, m_Controls->m_MovingChannelBox->value());
    else
      movingData = movingImage;

    mitk::Image::Pointer fixedImage = dynamic_cast<mitk::Image*>(fixedNode->GetData());
    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(fixedImage))
      targetData = mitk::DiffusionImageEditing::ExtractVolume(fixedImage, m_Controls->m_FixedChannelBox->value());
    else
      targetData = fixedImage;

    targetUID = mitk::EnsureUID(fixedNode->GetData());
    movingUID = mitk::EnsureUID(movingNode->GetData());
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Registration failed", QString("The images could not be prepared: ") + e.what());
    return;
  }

  QmitkRegistrationJob* pJob;

  if (m_Controls->m_RegBox->currentIndex()==0)
  {
    mitk::MultiModalRigidDefaultRegistrationAlgorithm< ItkFloatImageType >::Pointer algo = mitk::MultiModalRigidDefaultRegistrationAlgorithm< ItkFloatImageType >::New();
    pJob = new QmitkRegistrationJob(algo);
    m_RegistrationType = 0;
  }
  else
  {
    mitk::MultiModalAffineDefaultRegistrationAlgorithm< ItkFloatImageType >::Pointer algo = mitk::MultiModalAffineDefaultRegistrationAlgorithm< ItkFloatImageType >::New();
    pJob = new QmitkRegistrationJob(algo);
    m_RegistrationType = 1;
  }

  pJob->setAutoDelete(true);

  m_MovingImageNode = movingNode;
  pJob->m_spMovingData = movingData;
  pJob->m_spTargetData = targetData;
  pJob->m_TargetDataUID = targetUID;
  pJob->m_MovingDataUID = movingUID;

  connect(pJob, &QmitkRegistrationJob::RegResultIsAvailable, this, &QmitkSimpleRegistrationView::OnRegResultIsAvailable, Qt::BlockingQueuedConnection);
  connect(pJob, &QmitkRegistrationJob::Error, this, &QmitkSimpleRegistrationView::OnRegistrationError);
  connect(pJob, &QmitkRegistrationJob::Finished, this, &QmitkSimpleRegistrationView::OnRegistrationFinished);

  m_RegistrationRunning = true;
  m_Controls->m_RegistrationStartButton->setEnabled(false);
  m_Controls->m_RegistrationStartButton->setText("Registration in progress ...");

  QThreadPool* threadPool = QThreadPool::globalInstance();
  threadPool->start(pJob);
}

void QmitkSimpleRegistrationView::OnRegistrationError(QString message)
{
  QMessageBox::warning(nullptr, "Registration failed", message);
}

void QmitkSimpleRegistrationView::OnRegistrationFinished()
{
  m_RegistrationRunning = false;
  m_MovingImageNode = nullptr;
  m_Controls->m_RegistrationStartButton->setText("Start Registration");
  UpdateStartButton();
}

void QmitkSimpleRegistrationView::OnRegResultIsAvailable(mitk::MAPRegistrationWrapper::Pointer spResultRegistration, const QmitkRegistrationJob* job)
{
  mitk::Image::Pointer movingImage = dynamic_cast<mitk::Image*>(m_MovingImageNode->GetData());
  mitk::Image::Pointer image;

  try
  {
    // A rigid result can be applied to the geometry alone, without resampling.
    if (m_RegistrationType==0 && !m_Controls->m_ResampleBox->isChecked())
      image = mitk::DiffusionImageMapping::RefineGeometry(movingImage, spResultRegistration);
    else
      image = mitk::DiffusionImageMapping::MapImage(movingImage, spResultRegistration, job->m_spTargetData->GetGeometry(), mitk::ImageMappingInterpolator::BSpline_3);
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Registration failed", QString("The registration could not be applied: ") + e.what());
    return;
  }
  mitk::DiffusionPropertyHelper::CopyDICOMProperties(movingImage, image);

  mitk::DataNode::Pointer resultNode = mitk::DataNode::New();
  resultNode->SetData(image);

  m_MovingImageNode->SetVisibility(false);
  QString name = m_MovingImageNode->GetName().c_str();
  if (m_RegistrationType==0)
    resultNode->SetName((name+"_registered (rigid)").toStdString().c_str());
  else
    resultNode->SetName((name+"_registered (affine)").toStdString().c_str());
  GetDataStorage()->Add(resultNode);

  mitk::RenderingManager::GetInstance()->InitializeViews( resultNode->GetData()->GetTimeGeometry(),
                                                          mitk::RenderingManager::REQUEST_UPDATE_ALL,
                                                          true);

  if (m_Controls->m_RegOutputBox->isChecked())
  {
    mitk::DataNode::Pointer registration_node = mitk::DataNode::New();
    registration_node->SetData(spResultRegistration);
    if (m_RegistrationType==0)
      registration_node->SetName("Registration Object (rigid)");
    else
      registration_node->SetName("Registration Object (affine)");
    GetDataStorage()->Add(registration_node, resultNode);
  }

  // Also without an open render window part, where GetRenderWindowPart() is nullptr.
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();

  TractoChanged();
}

void QmitkSimpleRegistrationView::CreateQtPartControl( QWidget *parent )
{
  if ( !m_Controls )
  {
    m_Controls = std::make_unique<Ui::QmitkSimpleRegistrationViewControls>();
    QmitkDiffusionViewControlsGuard::Install(parent, this);
    m_Controls->setupUi( parent );

    m_Controls->m_FixedImageBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_MovingImageBox->SetDataStorage(this->GetDataStorage());

    mitk::TNodePredicateDataType<mitk::Image>::Pointer isImagePredicate = mitk::TNodePredicateDataType<mitk::Image>::New();
    m_Controls->m_FixedImageBox->SetPredicate(isImagePredicate);
    m_Controls->m_MovingImageBox->SetPredicate(isImagePredicate);

    mitk::TNodePredicateDataType<mitk::FiberBundle>::Pointer isFib = mitk::TNodePredicateDataType<mitk::FiberBundle>::New();
    mitk::TNodePredicateDataType<mitk::MAPRegistrationWrapper>::Pointer isReg = mitk::TNodePredicateDataType<mitk::MAPRegistrationWrapper>::New();

    m_Controls->m_TractoBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_RegObjectBox->SetDataStorage(this->GetDataStorage());
    m_Controls->m_TractoBox->SetPredicate(isFib);
    m_Controls->m_RegObjectBox->SetPredicate(isReg);

    connect(m_Controls->m_FixedImageBox, &QComboBox::currentIndexChanged, this, &QmitkSimpleRegistrationView::FixedImageChanged);
    connect(m_Controls->m_MovingImageBox, &QComboBox::currentIndexChanged, this, &QmitkSimpleRegistrationView::MovingImageChanged);

    connect(m_Controls->m_TractoBox, &QComboBox::currentIndexChanged, this, &QmitkSimpleRegistrationView::TractoChanged);
    connect(m_Controls->m_RegObjectBox, &QComboBox::currentIndexChanged, this, &QmitkSimpleRegistrationView::TractoChanged);

    connect(m_Controls->m_RegistrationStartButton, &QAbstractButton::clicked, this, &QmitkSimpleRegistrationView::StartRegistration);
    connect(m_Controls->m_TractoRegistrationStartButton, &QAbstractButton::clicked, this, &QmitkSimpleRegistrationView::StartTractoRegistration);

    FixedImageChanged();
    MovingImageChanged();
    TractoChanged();
  }
}

void QmitkSimpleRegistrationView::StartTractoRegistration()
{
  mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_Controls->m_TractoBox->GetSelectedNode()->GetData());
  mitk::MAPRegistrationWrapper::Pointer reg = dynamic_cast<mitk::MAPRegistrationWrapper*>(m_Controls->m_RegObjectBox->GetSelectedNode()->GetData());

  // Only affine registrations have a matrix; the box also offers deformable ones made elsewhere.
  mitk::MITKRegistrationHelper::Affine3DTransformType::Pointer affine = mitk::MITKRegistrationHelper::getAffineMatrix(reg, false);
  if (affine.IsNull())
  {
    QMessageBox::warning(nullptr, "Tractogram not transformed", "Only affine registrations can be applied to tractograms.");
    return;
  }

  mitk::FiberBundle::Pointer fib_copy = fib->Clone();
  fib_copy->TransformFibers(affine);
  mitk::DiffusionPropertyHelper::CopyProperties(fib, fib_copy);

  mitk::DataNode::Pointer registration_node = mitk::DataNode::New();
  registration_node->SetData(fib_copy);
  QString name = m_Controls->m_TractoBox->GetSelectedNode()->GetName().c_str();
  registration_node->SetName((name+"_registered").toStdString().c_str());
  GetDataStorage()->Add(registration_node, m_Controls->m_TractoBox->GetSelectedNode());
}

void QmitkSimpleRegistrationView::TractoChanged()
{
  if (m_Controls->m_RegObjectBox->GetSelectedNode().IsNotNull() && m_Controls->m_TractoBox->GetSelectedNode().IsNotNull())
    m_Controls->m_TractoRegistrationStartButton->setEnabled(true);
  else
    m_Controls->m_TractoRegistrationStartButton->setEnabled(false);
}

void QmitkSimpleRegistrationView::FixedImageChanged()
{
  if (m_Controls->m_FixedImageBox->GetSelectedNode().IsNotNull())
  {
    mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(m_Controls->m_FixedImageBox->GetSelectedNode()->GetData());
    int fourth_dim_size = image->GetTimeSteps();
    bool isdiff = mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image);

    if (isdiff)
    {
      m_Controls->m_FixedChannelBox->setEnabled(true);
      if (fourth_dim_size>1)
        m_Controls->m_FixedChannelBox->setMaximum(fourth_dim_size-1);
      else if (isdiff)
        m_Controls->m_FixedChannelBox->setMaximum(mitk::DiffusionPropertyHelper::GetGradientContainer(image)->Size()-1);
    }
    else
    {
      m_Controls->m_FixedChannelBox->setEnabled(false);
    }
  }
  else
  {
    m_Controls->m_FixedChannelBox->setEnabled(false);
  }
  UpdateStartButton();
}

void QmitkSimpleRegistrationView::MovingImageChanged()
{
  if (m_Controls->m_MovingImageBox->GetSelectedNode().IsNotNull())
  {
    mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(m_Controls->m_MovingImageBox->GetSelectedNode()->GetData());
    int fourth_dim_size = image->GetTimeSteps();
    bool isdiff = mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image);

    if (isdiff)
    {
      m_Controls->m_MovingChannelBox->setEnabled(true);
      if (fourth_dim_size>1)
        m_Controls->m_MovingChannelBox->setMaximum(fourth_dim_size-1);
      else if (isdiff)
        m_Controls->m_MovingChannelBox->setMaximum(mitk::DiffusionPropertyHelper::GetGradientContainer(image)->Size()-1);
    }
    else
    {
      m_Controls->m_MovingChannelBox->setEnabled(false);
    }
  }
  else
  {
    m_Controls->m_MovingChannelBox->setEnabled(false);
  }
  UpdateStartButton();
}

void QmitkSimpleRegistrationView::UpdateStartButton()
{
  // 3D images with one channel, and DWIs, of which one volume is registered.
  auto canRegister = [](const mitk::DataNode* node)
  {
    const mitk::Image* image = node != nullptr ? dynamic_cast<const mitk::Image*>(node->GetData()) : nullptr;
    if (image == nullptr)
      return false;
    return mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image) || (image->GetDimension() != 4 && image->GetNumberOfChannels() == 1);
  };
  m_Controls->m_RegistrationStartButton->setEnabled(!m_RegistrationRunning
                                                    && canRegister(m_Controls->m_FixedImageBox->GetSelectedNode())
                                                    && canRegister(m_Controls->m_MovingImageBox->GetSelectedNode()));
}

void QmitkSimpleRegistrationView::OnSelectionChanged(berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer>& )
{
  FixedImageChanged();
  MovingImageChanged();
  TractoChanged();
}


void QmitkSimpleRegistrationView::SetFocus()
{
  m_Controls->m_RegistrationStartButton->setFocus();
  FixedImageChanged();
  MovingImageChanged();
}

