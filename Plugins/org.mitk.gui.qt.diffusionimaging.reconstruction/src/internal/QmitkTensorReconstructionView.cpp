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

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkImageCast.h>
#include <mitkImageVtkMapper2D.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkStatusBar.h>

#include <itkB0ImageExtractionImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

const std::string QmitkTensorReconstructionView::VIEW_ID = "org.mitk.views.tensorreconstruction";

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
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkTensorReconstructionView::Reconstruct()
{
  QmitkReconstructionResults::Run("Reconstruction not possible:", [this]() {
    auto dwiNode = m_Controls->m_DwiBox->GetSelectedNode();
    if (dwiNode.IsNull())
      return;

    // The combo box lists ITK's linear least squares, then the eigenvalue correction. The corrected
    // fit always zeroes the tensors whose eigenvalues it could not make positive; the check box
    // acts on the linear fit only.
    const bool corrected = m_Controls->m_ReconctructionMethodBox->currentIndex() == 1;
    const auto fit = corrected ? mitk::TensorFit::EigenvalueCorrection : mitk::TensorFit::LinearLeastSquares;
    auto tensors = mitk::ReconstructTensors(dynamic_cast<mitk::Image *>(dwiNode->GetData()),
                                            fit,
                                            m_Controls->m_TensorReconstructionThreshold->value(),
                                            corrected || m_Controls->m_CheckNegativeEigenvalues->isChecked());

    this->GetDataStorage()->Add(ResultNode(tensors, dwiNode, corrected ? "_EigenvalueCorrected_DT" : "_LinearLeastSquares_DT"), dwiNode);
    this->RequestRenderWindowUpdate();
  });
}

void QmitkTensorReconstructionView::TensorsToOdf()
{
  QmitkReconstructionResults::Run("ODF estimation failed:", [this]() {
    auto tensorNode = m_Controls->m_OdfBox->GetSelectedNode();
    auto tensors = tensorNode.IsNotNull() ? dynamic_cast<mitk::TensorImage *>(tensorNode->GetData()) : nullptr;
    if (tensors == nullptr)
      return;

    auto odfs = mitk::convert::GetOdfFromTensorImage(tensors);
    mitk::DiffusionPropertyHelper::CopyDICOMProperties(tensors, odfs);
    this->GetDataStorage()->Add(ResultNode(odfs, tensorNode, "_Odf"), tensorNode);
  });
}

void QmitkTensorReconstructionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &)
{
  UpdateGui();
}

void QmitkTensorReconstructionView::UpdateGui()
{
  m_Controls->m_StartReconstruction->setEnabled(m_Controls->m_DwiBox->GetSelectedNode().IsNotNull());
  m_Controls->m_TensorsToDWIButton->setEnabled(m_Controls->m_DtiBox->GetSelectedNode().IsNotNull());
  m_Controls->m_TensorsToOdfButton->setEnabled(m_Controls->m_OdfBox->GetSelectedNode().IsNotNull());
}

void QmitkTensorReconstructionView::TensorsToDwi()
{
  QmitkReconstructionResults::Run("DWI estimation failed:", [this]() {
    auto tensorNode = m_Controls->m_DtiBox->GetSelectedNode();
    auto tensors = tensorNode.IsNotNull() ? dynamic_cast<mitk::TensorImage *>(tensorNode->GetData()) : nullptr;
    if (tensors == nullptr)
      return;

    // The combo box lists the sizes of the direction tables.
    auto dwi = mitk::TensorsToDwi(tensors,
                                  m_Controls->m_TensorsToDWINumDirsSelect->currentText().toUInt(),
                                  m_Controls->m_TensorsToDWIBValueEdit->value());
    auto node = ResultNode(dwi, tensorNode, "_DWI");
    mitk::ImageVtkMapper2D::SetDefaultProperties(node);
    this->GetDataStorage()->Add(node, tensorNode);
    this->RequestRenderWindowUpdate();
  });
}

void QmitkTensorReconstructionView::PreviewThreshold(int threshold)
{
  auto dwiNode = m_Controls->m_DwiBox->GetSelectedNode();
  auto dwi = dwiNode.IsNotNull() ? dynamic_cast<mitk::Image *>(dwiNode->GetData()) : nullptr;
  if (dwi == nullptr)
    return;

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

  using ImageType = itk::Image<short, 3>;
  auto filterThreshold = itk::BinaryThresholdImageFilter<ImageType, ImageType>::New();
  filterThreshold->SetInput(filterB0->GetOutput());
  filterThreshold->SetLowerThreshold(threshold);
  filterThreshold->SetInsideValue(0);
  filterThreshold->SetOutsideValue(1); // mark cut off values red
  filterThreshold->Update();

  auto mitkImage = mitk::Image::New();
  mitkImage->InitializeByItk(filterThreshold->GetOutput());
  mitkImage->SetVolume(filterThreshold->GetOutput()->GetBufferPointer());

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
