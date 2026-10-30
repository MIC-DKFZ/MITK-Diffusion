/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkQBallReconstructionView.h"

#include "QmitkReconstructionResults.h"

#include <QmitkDiffusionViewControlsGuard.h>

#include <QCheckBox>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkExceptionMacro.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkShImage.h>

#include <array>

const std::string QmitkQBallReconstructionView::VIEW_ID = "org.mitk.views.qballreconstruction";

using QmitkReconstructionResults::AddResult;
using QmitkReconstructionResults::ResultNode;

namespace
{
  /** Method and node name suffix of each entry of m_QBallReconstructionMethodComboBox. */
  struct MethodEntry
  {
    mitk::QballMethod method;
    const char *suffix;
  };
  constexpr std::array<MethodEntry, 6> Methods = {{{mitk::QballMethod::Numerical, "_Numerical_Qball"},
                                                   {mitk::QballMethod::Standard, "_Qball"},
                                                   {mitk::QballMethod::SolidAngle, "_CSA_Qball"},
                                                   {mitk::QballMethod::AdcProfile, "_AdcProfile"},
                                                   {mitk::QballMethod::RawSignal, "_RawSignal"},
                                                   {mitk::QballMethod::MultiShell, "_SH_MultiShell_Qball"}}};
}

QmitkQBallReconstructionView::~QmitkQBallReconstructionView()
{
  // The generated Ui class owns no widgets (they belong to the parent widget), only itself.
  delete m_Controls;
}

void QmitkQBallReconstructionView::CreateQtPartControl(QWidget *parent)
{
  if (!m_Controls)
  {
    m_Controls = new Ui::QmitkQBallReconstructionViewControls;
    QmitkDiffusionViewControlsGuard::Install(parent, this);
    m_Controls->setupUi(parent);
    this->CreateConnections();

    m_Controls->m_QBallReconstructionMaxLLevelComboBox->addItems({"2", "4", "6", "8", "10", "12"});
    m_Controls->m_QBallReconstructionMaxLLevelComboBox->setCurrentIndex(1);
    MethodChosen(m_Controls->m_QBallReconstructionMethodComboBox->currentIndex());
  }
}

void QmitkQBallReconstructionView::CreateConnections()
{
  connect(m_Controls->m_ButtonStandard, &QAbstractButton::clicked, this, &QmitkQBallReconstructionView::Reconstruct);
  connect(m_Controls->m_QBallReconstructionMethodComboBox, &QComboBox::currentIndexChanged, this, &QmitkQBallReconstructionView::MethodChosen);
  connect(m_Controls->m_ConvertButton, &QAbstractButton::clicked, this, &QmitkQBallReconstructionView::ConvertShImage);

  m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateIsDWI::New());

  m_Controls->m_ShImageBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_ShImageBox->SetPredicate(mitk::TNodePredicateDataType<mitk::ShImage>::New());

  connect(m_Controls->m_ImageBox, &QComboBox::currentIndexChanged, this, &QmitkQBallReconstructionView::UpdateGui);
  connect(m_Controls->m_ShImageBox, &QComboBox::currentIndexChanged, this, &QmitkQBallReconstructionView::UpdateGui);

  UpdateGui();
}

template <typename Work, typename Done>
void QmitkQBallReconstructionView::RunInBackground(const char *title, Work work, Done done)
{
  QmitkReconstructionResults::RunInBackground(this, m_Busy, [this]() { this->UpdateGui(); }, title, std::move(work), std::move(done));
}

void QmitkQBallReconstructionView::ConvertShImage()
{
  auto shNode = m_Controls->m_ShImageBox->GetSelectedNode();
  mitk::Image::ConstPointer shImage = shNode.IsNotNull() ? dynamic_cast<mitk::Image *>(shNode->GetData()) : nullptr;
  if (m_Busy || shImage.IsNull())
    return;

  this->RunInBackground(
    "Conversion not possible:",
    [shImage]() {
      auto odfs = mitk::convert::GetOdfFromShImage(shImage);
      mitk::DiffusionPropertyHelper::CopyDICOMProperties(shImage, odfs);
      return odfs;
    },
    [this, shNode](mitk::OdfImage::Pointer odfs) { AddResult(this->GetDataStorage(), ResultNode(odfs, shNode, "_Odf"), shNode); });
}

void QmitkQBallReconstructionView::UpdateGui()
{
  // One action at a time: the buttons stay disabled while one runs in the background.
  m_Controls->m_ButtonStandard->setEnabled(!m_Busy && m_Controls->m_ImageBox->GetSelectedNode().IsNotNull());
  GenerateShellSelectionUI(m_Controls->m_ImageBox->GetSelectedNode());

  m_Controls->m_ConvertButton->setEnabled(!m_Busy && m_Controls->m_ShImageBox->GetSelectedNode().IsNotNull());
}

void QmitkQBallReconstructionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> & /*nodes*/)
{
  UpdateGui();
}

void QmitkQBallReconstructionView::MethodChosen(int method)
{
  m_Controls->m_QBallSelectionBox->setHidden(true);
  m_Controls->m_OutputOdfImage->setHidden(true);

  // The B0 threshold in m_ShFrame applies to every method, SH order and lambda only to the SH
  // methods.
  const bool shMethod = method != 0;
  m_Controls->m_QBallReconstructionMaxLLevelTextLabel_2->setVisible(shMethod);
  m_Controls->m_QBallReconstructionMaxLLevelComboBox->setVisible(shMethod);
  m_Controls->m_QBallReconstructionLambdaTextLabel_2->setVisible(shMethod);
  m_Controls->m_QBallReconstructionLambdaLineEdit->setVisible(shMethod);

  // The cases follow the entries of m_QBallReconstructionMethodComboBox, as Methods does.
  switch (method)
  {
    case 0:
      m_Controls->m_Description->setText("Numerical recon. (Tuch 2004)");
      break;
    case 1:
      m_Controls->m_Description->setText("Spherical harmonics recon. (Descoteaux 2007)");
      m_Controls->m_OutputOdfImage->setHidden(false);
      break;
    case 2:
      m_Controls->m_Description->setText("SH recon. with solid angle consideration (Aganj 2009)");
      m_Controls->m_OutputOdfImage->setHidden(false);
      break;
    case 3:
      m_Controls->m_Description->setText("SH recon. of the plain ADC-profiles");
      m_Controls->m_OutputOdfImage->setHidden(false);
      break;
    case 4:
      m_Controls->m_Description->setText("SH recon. of the raw diffusion signal");
      m_Controls->m_OutputOdfImage->setHidden(false);
      break;
    case 5:
      m_Controls->m_Description->setText("SH recon. of the multi shell diffusion signal (Aganj 2010)");
      m_Controls->m_QBallSelectionBox->setHidden(false);
      m_Controls->m_OutputOdfImage->setHidden(false);
      break;
  }
}

void QmitkQBallReconstructionView::Reconstruct()
{
  auto dwiNode = m_Controls->m_ImageBox->GetSelectedNode();
  mitk::Image::ConstPointer dwi = dwiNode.IsNotNull() ? dynamic_cast<mitk::Image *>(dwiNode->GetData()) : nullptr;
  const int index = m_Controls->m_QBallReconstructionMethodComboBox->currentIndex();
  if (m_Busy || dwi.IsNull() || index < 0 || index >= static_cast<int>(Methods.size()))
    return;
  const auto &entry = Methods[index];

  mitk::QballParameters parameters;
  parameters.method = entry.method;
  parameters.shOrder = m_Controls->m_QBallReconstructionMaxLLevelComboBox->currentText().toUInt();
  parameters.lambda = m_Controls->m_QBallReconstructionLambdaLineEdit->value();
  parameters.b0Threshold = static_cast<short>(m_Controls->m_QBallReconstructionThreasholdEdit->value());
  parameters.sampledOdf = m_Controls->m_OutputOdfImage->isChecked();
  // The b-values of the checked boxes, built for the selected DWI; ReconstructQball rejects a
  // b-value the DWI does not have.
  for (unsigned int i = 0; i < m_ShellCheckBoxes.size(); ++i)
  {
    if (m_ShellCheckBoxes[i]->isChecked())
      parameters.shells.push_back(m_ShellSelectionBValues[i]);
  }
  if (entry.method == mitk::QballMethod::MultiShell && parameters.shells.empty())
  {
    QmitkReconstructionResults::Run("Reconstruction not possible:", []() { mitkThrow() << "Select three equidistant shells."; });
    return;
  }

  const std::string suffix = entry.suffix;
  this->RunInBackground(
    "Reconstruction not possible:",
    [dwi, parameters]() { return mitk::ReconstructQball(dwi, parameters); },
    [this, dwiNode, suffix](const mitk::QballResult &result) {
      if (result.sh.IsNull())
      {
        AddResult(this->GetDataStorage(), ResultNode(result.odf, dwiNode, suffix), dwiNode);
      }
      else
      {
        auto shNode = ResultNode(result.sh, dwiNode, suffix + "_SH");
        AddResult(this->GetDataStorage(), shNode, dwiNode);
        if (result.odf.IsNotNull())
          this->GetDataStorage()->Add(ResultNode(result.odf, dwiNode, suffix + "_Sampled"), shNode);
      }
      this->RequestRenderWindowUpdate();
    });
}

void QmitkQBallReconstructionView::GenerateShellSelectionUI(mitk::DataNode::Pointer node)
{
  // One check box per non-zero shell of the selected DWI. Only the b-values are kept, never the
  // DWI. The boxes (and their check states) stay while the same node with the same shells is
  // selected; the node address is only compared.
  std::vector<unsigned int> bValues;
  auto dwi = node.IsNotNull() ? dynamic_cast<mitk::Image *>(node->GetData()) : nullptr;
  if (dwi != nullptr && mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi))
  {
    for (const auto &shell : mitk::DiffusionPropertyHelper::GetBValueMap(dwi))
    {
      if (shell.first != 0)
        bValues.push_back(shell.first);
    }
  }
  if (node.GetPointer() == m_ShellSelectionNode && bValues == m_ShellSelectionBValues)
    return;

  for (auto box : m_ShellCheckBoxes)
    delete box;
  m_ShellCheckBoxes.clear();
  m_ShellSelectionNode = node.GetPointer();
  m_ShellSelectionBValues = bValues;

  for (auto bValue : bValues)
  {
    auto box = new QCheckBox(QString::number(bValue));
    box->setChecked(true);
    m_Controls->m_QBallSelectionBox->layout()->addWidget(box);
    m_ShellCheckBoxes.push_back(box);
  }
}
