/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkDiffusionQuantificationView.h"

#include "QmitkReconstructionResults.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkNodePredicateOr.h>
#include <mitkOdfImage.h>
#include <mitkShImage.h>
#include <mitkStatusBar.h>
#include <mitkTensorImage.h>

const std::string QmitkDiffusionQuantificationView::VIEW_ID = "org.mitk.views.diffusionquantification";

using QmitkReconstructionResults::AddResult;
using QmitkReconstructionResults::ResultNode;
using QmitkReconstructionResults::ShowAsIndexMap;

namespace
{
  /** Node name suffix and level window (for unscaled values) of each tensor index. */
  struct TensorIndexDisplay
  {
    const char *name;
    double level;
    double window;
  };

  TensorIndexDisplay DisplayOf(mitk::TensorIndex index)
  {
    switch (index)
    {
      case mitk::TensorIndex::FA:
        return {"FA", 0.5, 1.0};
      case mitk::TensorIndex::RA:
        // Basser and Pierpaoli's RA ranges from 0 to sqrt(2).
        return {"RA", 0.7, 1.4};
      case mitk::TensorIndex::AD:
        return {"AD", 0.0015, 0.003};
      case mitk::TensorIndex::RD:
        return {"RD", 0.0015, 0.003};
      case mitk::TensorIndex::CA:
        return {"CA", 0.5, 1.0};
      case mitk::TensorIndex::MD:
        break;
    }
    return {"MD", 0.0015, 0.003};
  }
}

QmitkDiffusionQuantificationView::~QmitkDiffusionQuantificationView()
{
  // The generated Ui class owns no widgets (they belong to the parent widget), only itself.
  delete m_Controls;
}

void QmitkDiffusionQuantificationView::CreateQtPartControl(QWidget *parent)
{
  if (m_Controls)
    return;

  m_Controls = new Ui::QmitkDiffusionQuantificationViewControls;
  m_Controls->setupUi(parent);

  connect(m_Controls->m_GFAButton, &QAbstractButton::clicked, this, [this]() { ComputeGfa(); });
  connect(m_Controls->m_FAButton, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::FA); });
  connect(m_Controls->m_RAButton, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::RA); });
  connect(m_Controls->m_ADButton, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::AD); });
  connect(m_Controls->m_RDButton, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::RD); });
  connect(m_Controls->m_MDButton, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::MD); });
  connect(m_Controls->m_ClusteringAnisotropy, &QAbstractButton::clicked, this, [this]() { ComputeTensorIndex(mitk::TensorIndex::CA); });
  connect(m_Controls->m_MdDwiButton, &QAbstractButton::clicked, this, [this]() { ComputeAdc(false); });
  connect(m_Controls->m_AdcDwiButton, &QAbstractButton::clicked, this, [this]() { ComputeAdc(true); });

  m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
  auto isDti = mitk::TNodePredicateDataType<mitk::TensorImage>::New();
  auto isOdf = mitk::TNodePredicateDataType<mitk::OdfImage>::New();
  auto isSh = mitk::TNodePredicateDataType<mitk::ShImage>::New();
  m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateOr::New(
    isSh, mitk::NodePredicateOr::New(isDti, mitk::NodePredicateOr::New(isOdf, mitk::NodePredicateIsDWI::New()))));

  connect(m_Controls->m_ImageBox, &QComboBox::currentIndexChanged, this, &QmitkDiffusionQuantificationView::UpdateGui);
}

void QmitkDiffusionQuantificationView::SetFocus()
{
  m_Controls->m_ScaleImageValuesBox->setFocus();
}

void QmitkDiffusionQuantificationView::UpdateGui()
{
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  auto data = node.IsNotNull() ? node->GetData() : nullptr;
  const bool odf = dynamic_cast<mitk::OdfImage *>(data) != nullptr || dynamic_cast<mitk::ShImage *>(data) != nullptr;
  const bool tensor = dynamic_cast<mitk::TensorImage *>(data) != nullptr;
  const bool dwi = mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dynamic_cast<mitk::Image *>(data));
  // One computation at a time: the buttons stay disabled while one runs in the background.
  const bool idle = !m_Busy;

  m_Controls->m_GFAButton->setEnabled(idle && odf);
  m_Controls->m_FAButton->setEnabled(idle && tensor);
  m_Controls->m_RAButton->setEnabled(idle && tensor);
  m_Controls->m_ADButton->setEnabled(idle && tensor);
  m_Controls->m_RDButton->setEnabled(idle && tensor);
  m_Controls->m_MDButton->setEnabled(idle && tensor);
  m_Controls->m_ClusteringAnisotropy->setEnabled(idle && tensor);
  m_Controls->m_AdcDwiButton->setEnabled(idle && dwi);
  m_Controls->m_MdDwiButton->setEnabled(idle && dwi);
}

void QmitkDiffusionQuantificationView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &)
{
  UpdateGui();
}

template <typename Work>
void QmitkDiffusionQuantificationView::RunInBackground(Work work, mitk::DataNode::Pointer node, const std::string &suffix, double level, double window)
{
  QmitkReconstructionResults::RunInBackground(
    this, m_Busy, [this]() { this->UpdateGui(); }, "Computation not possible:", std::move(work),
    [this, node, suffix, level, window](mitk::Image::Pointer image) {
      auto result = ResultNode(image, node, suffix);
      ShowAsIndexMap(result, level, window);
      AddResult(this->GetDataStorage(), result, node);
      mitk::StatusBar::GetInstance()->DisplayText("Computation complete.");
      this->RequestRenderWindowUpdate();
    });
}

void QmitkDiffusionQuantificationView::ComputeAdc(bool fit)
{
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  mitk::Image::ConstPointer dwi = node.IsNotNull() ? dynamic_cast<mitk::Image *>(node->GetData()) : nullptr;
  if (m_Busy || dwi.IsNull())
    return;

  const double scale = m_Controls->m_ScaleImageValuesBox->value();
  // The window fits diffusivities in mm^2/s; the image holds them times the scale factor.
  this->RunInBackground([dwi, fit, scale]() { return mitk::AdcImage(dwi, fit, scale); }, node, fit ? "_ADC" : "_MD", 0.0015 * scale, 0.003 * scale);
}

void QmitkDiffusionQuantificationView::ComputeGfa()
{
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  mitk::Image::ConstPointer image = node.IsNotNull() ? dynamic_cast<mitk::Image *>(node->GetData()) : nullptr;
  if (m_Busy || image.IsNull())
    return;

  mitk::StatusBar::GetInstance()->DisplayText(QString("Computing GFA for %1").arg(QString::fromStdString(node->GetName())).toLatin1());
  const double scale = m_Controls->m_ScaleImageValuesBox->value();
  // The window fits GFA values in [0, 1]; the image holds them times the scale factor.
  this->RunInBackground([image, scale]() { return mitk::GfaImage(image, scale); }, node, "_GFA", 0.5 * scale, 1.0 * scale);
}

void QmitkDiffusionQuantificationView::ComputeTensorIndex(mitk::TensorIndex index)
{
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  mitk::TensorImage::ConstPointer tensors = node.IsNotNull() ? dynamic_cast<mitk::TensorImage *>(node->GetData()) : nullptr;
  if (m_Busy || tensors.IsNull())
    return;

  const auto display = DisplayOf(index);
  mitk::StatusBar::GetInstance()->DisplayText(QString("Computing %1 for %2").arg(QString(display.name), QString::fromStdString(node->GetName())).toLatin1());
  const double scale = m_Controls->m_ScaleImageValuesBox->value();
  // The windows fit unscaled values; the image holds them times the scale factor.
  this->RunInBackground([tensors, index, scale]() { return mitk::TensorIndexImage(tensors, index, scale); },
                        node,
                        std::string("_") + display.name,
                        display.level * scale,
                        display.window * scale);
}
