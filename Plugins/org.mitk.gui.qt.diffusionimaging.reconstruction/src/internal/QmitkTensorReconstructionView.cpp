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

#include "QmitkReconstructionResults.h"

#include <QmitkDiffusionViewControlsGuard.h>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkImageVtkMapper2D.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkStatusBar.h>

#include <itkB0ImageExtractionImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

const std::string QmitkTensorReconstructionView::VIEW_ID = "org.mitk.views.tensorreconstruction";

using QmitkReconstructionResults::AddResult;
using QmitkReconstructionResults::ResultNode;

QmitkTensorReconstructionView::~QmitkTensorReconstructionView()
{
  // The generated Ui class owns no widgets (they belong to the parent widget), only itself.
  delete m_Controls;
}

void QmitkTensorReconstructionView::CreateQtPartControl(QWidget *parent)
{
  if (!m_Controls)
  {
    m_Controls = new Ui::QmitkTensorReconstructionViewControls;
    QmitkDiffusionViewControlsGuard::Install(parent, this);
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
  connect(m_Controls->m_StartReconstruction, &QAbstractButton::clicked, this, &QmitkTensorReconstructionView::Reconstruct);
  connect(m_Controls->m_Advanced1, &QAbstractButton::clicked, this, &QmitkTensorReconstructionView::Advanced1CheckboxClicked);
  connect(m_Controls->m_TensorsToDWIButton, &QAbstractButton::clicked, this, &QmitkTensorReconstructionView::TensorsToDwi);
  connect(m_Controls->m_TensorsToOdfButton, &QAbstractButton::clicked, this, &QmitkTensorReconstructionView::TensorsToOdf);
  connect(m_Controls->m_TensorReconstructionThreshold, &QSpinBox::valueChanged, this, &QmitkTensorReconstructionView::PreviewThreshold);

  m_Controls->m_DwiBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_DwiBox->SetPredicate(mitk::NodePredicateIsDWI::New());
  connect(m_Controls->m_DwiBox, &QComboBox::currentIndexChanged, this, &QmitkTensorReconstructionView::UpdateGui);

  auto isDti = mitk::TNodePredicateDataType<mitk::TensorImage>::New();
  m_Controls->m_OdfBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_OdfBox->SetPredicate(isDti);
  connect(m_Controls->m_OdfBox, &QComboBox::currentIndexChanged, this, &QmitkTensorReconstructionView::UpdateGui);

  m_Controls->m_DtiBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_DtiBox->SetPredicate(isDti);
  connect(m_Controls->m_DtiBox, &QComboBox::currentIndexChanged, this, &QmitkTensorReconstructionView::UpdateGui);
}

void QmitkTensorReconstructionView::Advanced1CheckboxClicked()
{
  m_Controls->frame->setVisible(m_Controls->m_Advanced1->isChecked());
}

void QmitkTensorReconstructionView::Deactivated()
{
  // Named: a temporary in the range expression would be released before the loop runs.
  const auto nodes = this->GetDataStorage()->GetAll();
  for (const auto &node : *nodes)
  {
    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dynamic_cast<mitk::Image *>(node->GetData())))
    {
      if (auto overlay = this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", node))
        this->GetDataStorage()->Remove(overlay);
    }
  }
  m_PreviewB0 = nullptr;
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

template <typename Work, typename Done>
void QmitkTensorReconstructionView::RunInBackground(const char *title, Work work, Done done)
{
  QmitkReconstructionResults::RunInBackground(m_Controls->m_DwiBox, m_Busy, [this]() { this->UpdateGui(); }, title, std::move(work), std::move(done));
}

void QmitkTensorReconstructionView::Reconstruct()
{
  auto dwiNode = m_Controls->m_DwiBox->GetSelectedNode();
  mitk::Image::ConstPointer dwi = dwiNode.IsNotNull() ? dynamic_cast<mitk::Image *>(dwiNode->GetData()) : nullptr;
  if (m_Busy || dwi.IsNull())
    return;

  // The combo box lists ITK's linear least squares, then the eigenvalue correction. The corrected
  // fit always zeroes the tensors whose eigenvalues it could not make positive; the check box
  // acts on the linear fit only.
  const bool corrected = m_Controls->m_ReconctructionMethodBox->currentIndex() == 1;
  const auto fit = corrected ? mitk::TensorFit::EigenvalueCorrection : mitk::TensorFit::LinearLeastSquares;
  const double threshold = m_Controls->m_TensorReconstructionThreshold->value();
  const bool zeroNegativeEigenvalues = corrected || m_Controls->m_CheckNegativeEigenvalues->isChecked();
  this->RunInBackground(
    "Reconstruction not possible:",
    [dwi, fit, threshold, zeroNegativeEigenvalues]() { return mitk::ReconstructTensors(dwi, fit, threshold, zeroNegativeEigenvalues); },
    [this, dwiNode, corrected](mitk::TensorImage::Pointer tensors) {
      AddResult(this->GetDataStorage(), ResultNode(tensors, dwiNode, corrected ? "_EigenvalueCorrected_DT" : "_LinearLeastSquares_DT"), dwiNode);
      this->RequestRenderWindowUpdate();
    });
}

void QmitkTensorReconstructionView::TensorsToOdf()
{
  auto tensorNode = m_Controls->m_OdfBox->GetSelectedNode();
  mitk::TensorImage::ConstPointer tensors = tensorNode.IsNotNull() ? dynamic_cast<mitk::TensorImage *>(tensorNode->GetData()) : nullptr;
  if (m_Busy || tensors.IsNull())
    return;

  this->RunInBackground(
    "ODF estimation failed:",
    [tensors]() {
      auto odfs = mitk::convert::GetOdfFromTensorImage(tensors);
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(tensors, odfs);
      return odfs;
    },
    [this, tensorNode](mitk::OdfImage::Pointer odfs) { AddResult(this->GetDataStorage(), ResultNode(odfs, tensorNode, "_Odf"), tensorNode); });
}

void QmitkTensorReconstructionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &)
{
  UpdateGui();
}

void QmitkTensorReconstructionView::UpdateGui()
{
  // One action at a time: the buttons stay disabled while one runs in the background.
  m_Controls->m_StartReconstruction->setEnabled(!m_Busy && m_Controls->m_DwiBox->GetSelectedNode().IsNotNull());
  m_Controls->m_TensorsToDWIButton->setEnabled(!m_Busy && m_Controls->m_DtiBox->GetSelectedNode().IsNotNull());
  m_Controls->m_TensorsToOdfButton->setEnabled(!m_Busy && m_Controls->m_OdfBox->GetSelectedNode().IsNotNull());
}

void QmitkTensorReconstructionView::TensorsToDwi()
{
  auto tensorNode = m_Controls->m_DtiBox->GetSelectedNode();
  mitk::TensorImage::ConstPointer tensors = tensorNode.IsNotNull() ? dynamic_cast<mitk::TensorImage *>(tensorNode->GetData()) : nullptr;
  if (m_Busy || tensors.IsNull())
    return;

  // The combo box lists the sizes of the direction tables.
  const unsigned int numberOfDirections = m_Controls->m_TensorsToDWINumDirsSelect->currentText().toUInt();
  const double bValue = m_Controls->m_TensorsToDWIBValueEdit->value();
  this->RunInBackground(
    "DWI estimation failed:",
    [tensors, numberOfDirections, bValue]() { return mitk::TensorsToDwi(tensors, numberOfDirections, bValue); },
    [this, tensorNode](mitk::Image::Pointer dwi) {
      auto node = ResultNode(dwi, tensorNode, "_DWI");
      mitk::ImageVtkMapper2D::SetDefaultProperties(node);
      AddResult(this->GetDataStorage(), node, tensorNode);
      this->RequestRenderWindowUpdate();
    });
}

void QmitkTensorReconstructionView::PreviewThreshold(int threshold)
{
  auto dwiNode = m_Controls->m_DwiBox->GetSelectedNode();
  const mitk::Image *dwi = dwiNode.IsNotNull() ? dynamic_cast<mitk::Image *>(dwiNode->GetData()) : nullptr;
  if (dwi == nullptr)
    return;

  // The mean b0 image is computed once per DWI (and again when the DWI changes), not on every step
  // of the threshold box. The DWI is read through a const image, which does not block rendering.
  if (m_PreviewB0.IsNull() || m_PreviewDwi.Lock().GetPointer() != dwi || m_PreviewDwiTime != dwi->GetMTime())
  {
    m_PreviewB0 = nullptr;
    using FilterType = itk::B0ImageExtractionImageFilter<short, short>;
    auto filterB0 = FilterType::New();
    try
    {
      mitk::DiffusionPropertyHelper::ImageType::Pointer itkDwi;
      mitk::CastToItkImage(dwi, itkDwi);
      filterB0->SetInput(itkDwi);
      filterB0->SetDirections(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi));
      filterB0->SetReferenceBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
      filterB0->Update();
    }
    catch (const std::exception &e)
    {
      // The preview runs on every step of the threshold box, so it reports in the status bar
      // instead of opening a dialog per key press. A DWI without b0 volume has no preview.
      if (auto overlay = this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", dwiNode))
        this->GetDataStorage()->Remove(overlay);
      mitk::StatusBar::GetInstance()->DisplayText((std::string("No threshold preview: ") + e.what()).c_str());
      mitk::RenderingManager::GetInstance()->RequestUpdateAll();
      return;
    }
    m_PreviewB0 = filterB0->GetOutput();
    m_PreviewB0->DisconnectPipeline();
    m_PreviewDwi = dwi;
    m_PreviewDwiTime = dwi->GetMTime();
  }

  using ImageType = itk::Image<short, 3>;
  auto filterThreshold = itk::BinaryThresholdImageFilter<ImageType, ImageType>::New();
  filterThreshold->SetInput(m_PreviewB0);
  filterThreshold->SetLowerThreshold(threshold);
  filterThreshold->SetInsideValue(0);
  filterThreshold->SetOutsideValue(1); // mark cut off values red
  filterThreshold->Update();

  auto mitkImage = mitk::GrabItkImageMemory(filterThreshold->GetOutput());

  mitk::DataNode::Pointer node = this->GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", dwiNode);
  if (node.IsNull())
  {
    node = mitk::DataNode::New();
    this->GetDataStorage()->Add(node, dwiNode);
    node->SetProperty("name", mitk::StringProperty::New("ThresholdOverlay"));
    node->SetBoolProperty("helper object", true);
  }
  node->SetData(mitkImage);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}
