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

  m_Controls->m_GFAButton->setEnabled(odf);
  m_Controls->m_FAButton->setEnabled(tensor);
  m_Controls->m_RAButton->setEnabled(tensor);
  m_Controls->m_ADButton->setEnabled(tensor);
  m_Controls->m_RDButton->setEnabled(tensor);
  m_Controls->m_MDButton->setEnabled(tensor);
  m_Controls->m_ClusteringAnisotropy->setEnabled(tensor);
  m_Controls->m_AdcDwiButton->setEnabled(dwi);
  m_Controls->m_MdDwiButton->setEnabled(dwi);
}

void QmitkDiffusionQuantificationView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &)
{
  UpdateGui();
}

void QmitkDiffusionQuantificationView::ComputeAdc(bool fit)
{
  QmitkReconstructionResults::Run("Computation not possible:", [this, fit]() {
    auto node = m_Controls->m_ImageBox->GetSelectedNode();
    if (node.IsNull())
      return;

    const double scale = m_Controls->m_ScaleImageValuesBox->value();
    auto result = ResultNode(mitk::AdcImage(dynamic_cast<mitk::Image *>(node->GetData()), fit, scale), node, fit ? "_ADC" : "_MD");
    // The window fits diffusivities in mm^2/s; the image holds them times the scale factor.
    ShowAsIndexMap(result, 0.0015 * scale, 0.003 * scale);
    this->GetDataStorage()->Add(result, node);
  });
}

void QmitkDiffusionQuantificationView::ComputeGfa()
{
  QmitkReconstructionResults::Run("Computation not possible:", [this]() {
    auto node = m_Controls->m_ImageBox->GetSelectedNode();
    if (node.IsNull())
      return;

    mitk::StatusBar::GetInstance()->DisplayText(QString("Computing GFA for %1").arg(QString::fromStdString(node->GetName())).toLatin1());
    const double scale = m_Controls->m_ScaleImageValuesBox->value();
    auto result = ResultNode(mitk::GfaImage(dynamic_cast<mitk::Image *>(node->GetData()), scale), node, "_GFA");
    // The window fits GFA values in [0, 1]; the image holds them times the scale factor.
    ShowAsIndexMap(result, 0.5 * scale, 1.0 * scale);
    this->GetDataStorage()->Add(result, node);

    mitk::StatusBar::GetInstance()->DisplayText("Computation complete.");
    this->RequestRenderWindowUpdate();
  });
}

void QmitkDiffusionQuantificationView::ComputeTensorIndex(mitk::TensorIndex index)
{
  QmitkReconstructionResults::Run("Computation not possible:", [this, index]() {
    auto node = m_Controls->m_ImageBox->GetSelectedNode();
    if (node.IsNull())
      return;

    const auto display = DisplayOf(index);
    mitk::StatusBar::GetInstance()->DisplayText(QString("Computing %1 for %2").arg(QString(display.name), QString::fromStdString(node->GetName())).toLatin1());
    const double scale = m_Controls->m_ScaleImageValuesBox->value();
    auto result = ResultNode(mitk::TensorIndexImage(dynamic_cast<mitk::TensorImage *>(node->GetData()), index, scale),
                             node,
                             std::string("_") + display.name);
    // The windows fit unscaled values; the image holds them times the scale factor.
    ShowAsIndexMap(result, display.level * scale, display.window * scale);
    this->GetDataStorage()->Add(result, node);

    mitk::StatusBar::GetInstance()->DisplayText("Computation complete.");
    this->RequestRenderWindowUpdate();
  });
}
