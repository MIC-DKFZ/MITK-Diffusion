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

#include <QCheckBox>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkExceptionMacro.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkShImage.h>

#include <array>

const std::string QmitkQBallReconstructionView::VIEW_ID = "org.mitk.views.qballreconstruction";

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

void QmitkQBallReconstructionView::ConvertShImage()
{
  QmitkReconstructionResults::Run("Conversion not possible:", [this]() {
    auto shNode = m_Controls->m_ShImageBox->GetSelectedNode();
    if (shNode.IsNull())
      return;
    auto shImage = dynamic_cast<mitk::Image *>(shNode->GetData());
    auto odfs = mitk::convert::GetOdfFromShImage(shImage);
    mitk::DiffusionPropertyHelper::CopyDICOMProperties(shImage, odfs);
    this->GetDataStorage()->Add(ResultNode(odfs, shNode, "_Odf"), shNode);
  });
}

void QmitkQBallReconstructionView::UpdateGui()
{
  m_Controls->m_ButtonStandard->setEnabled(m_Controls->m_ImageBox->GetSelectedNode().IsNotNull());
  GenerateShellSelectionUI(m_Controls->m_ImageBox->GetSelectedNode());

  m_Controls->m_ConvertButton->setEnabled(m_Controls->m_ShImageBox->GetSelectedNode().IsNotNull());
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
  QmitkReconstructionResults::Run("Reconstruction not possible:", [this]() {
    auto dwiNode = m_Controls->m_ImageBox->GetSelectedNode();
    const int index = m_Controls->m_QBallReconstructionMethodComboBox->currentIndex();
    if (dwiNode.IsNull() || index < 0 || index >= static_cast<int>(Methods.size()))
      return;
    const auto &entry = Methods[index];

    mitk::QballParameters parameters;
    parameters.method = entry.method;
    parameters.shOrder = m_Controls->m_QBallReconstructionMaxLLevelComboBox->currentText().toUInt();
    parameters.lambda = m_Controls->m_QBallReconstructionLambdaLineEdit->value();
    parameters.b0Threshold = static_cast<short>(m_Controls->m_QBallReconstructionThreasholdEdit->value());
    parameters.sampledOdf = m_Controls->m_OutputOdfImage->isChecked();
    // The shell selection is read from the DWI now; the check boxes belong to the selected DWI.
    for (unsigned int i = 0; i < m_ShellCheckBoxes.size(); ++i)
    {
      if (m_ShellCheckBoxes[i]->isChecked())
        parameters.shells.push_back(m_ShellSelectionBValues[i]);
    }
    if (entry.method == mitk::QballMethod::MultiShell && parameters.shells.empty())
      mitkThrow() << "Select three equidistant shells.";

    const auto result = mitk::ReconstructQball(dynamic_cast<mitk::Image *>(dwiNode->GetData()), parameters);

    const std::string suffix = entry.suffix;
    if (result.sh.IsNull())
    {
      this->GetDataStorage()->Add(ResultNode(result.odf, dwiNode, suffix), dwiNode);
    }
    else
    {
      auto shNode = ResultNode(result.sh, dwiNode, suffix + "_SH");
      this->GetDataStorage()->Add(shNode, dwiNode);
      if (result.odf.IsNotNull())
        this->GetDataStorage()->Add(ResultNode(result.odf, dwiNode, suffix + "_Sampled"), shNode);
    }
    this->RequestRenderWindowUpdate();
  });
}

void QmitkQBallReconstructionView::GenerateShellSelectionUI(mitk::DataNode::Pointer node)
{
  // One check box per non-zero shell of the selected DWI. Nothing keeps a pointer to the DWI: the
  // selection is read from the DWI at reconstruction time. The boxes (and their check states) stay
  // while the same node with the same shells is selected; the node address is only compared.
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
