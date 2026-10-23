/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkOdfMaximaExtractionView.h"

#include "QmitkReconstructionResults.h"

#include <mitkDiffusionPeaks.h>
#include <mitkExceptionMacro.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateIsDWI.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>

#include <algorithm>
#include <initializer_list>

using QmitkReconstructionResults::ResultNode;

QmitkOdfMaximaExtractionView::~QmitkOdfMaximaExtractionView()
{
  // The generated Ui class owns no widgets (they belong to the parent widget), only itself.
  delete m_Controls;
}

void QmitkOdfMaximaExtractionView::CreateQtPartControl(QWidget *parent)
{
  if (m_Controls)
    return;

  m_Controls = new Ui::QmitkOdfMaximaExtractionViewControls;
  m_Controls->setupUi(parent);

  connect(m_Controls->m_StartPeakExtractionButton, &QAbstractButton::clicked, this, &QmitkOdfMaximaExtractionView::StartPeakExtraction);

  m_Controls->m_MaskBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());

  auto isNoDwi = mitk::NodePredicateNot::New(mitk::NodePredicateIsDWI::New());
  auto isNoOdf = mitk::NodePredicateNot::New(mitk::NodePredicateDataType::New("OdfImage"));
  auto isBinary = mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
  m_Controls->m_MaskBox->SetPredicate(mitk::NodePredicateAnd::New(
    mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(isNoOdf, isNoDwi), mitk::NodePredicateDimension::New(3)), isBinary));
  m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateOr::New(mitk::TNodePredicateDataType<mitk::TensorImage>::New(),
                                                                  mitk::TNodePredicateDataType<mitk::ShImage>::New()));
  m_Controls->m_MaskBox->SetZeroEntryText("--");

  connect(m_Controls->m_ImageBox, &QmitkDataStorageComboBox::OnSelectionChanged, this, &QmitkOdfMaximaExtractionView::OnImageSelectionChanged);

  OnImageSelectionChanged();
}

void QmitkOdfMaximaExtractionView::StartPeakExtraction()
{
  QmitkReconstructionResults::Run("Peak extraction not possible:", [this]() {
    auto imageNode = m_Controls->m_ImageBox->GetSelectedNode();
    if (imageNode.IsNull())
      return;
    auto maskNode = m_Controls->m_MaskBox->GetSelectedNode();
    auto mask = maskNode.IsNotNull() ? dynamic_cast<mitk::Image *>(maskNode->GetData()) : nullptr;

    if (auto tensors = dynamic_cast<mitk::TensorImage *>(imageNode->GetData()))
    {
      // The absolute threshold acts as FA threshold; without normalization the vectors are scaled by FA.
      const auto result = mitk::ExtractTensorPeaks(tensors,
                                                   m_Controls->m_AbsoluteThresholdBox->value(),
                                                   m_Controls->m_NormalizationBox->currentIndex() != 0,
                                                   mask);
      this->GetDataStorage()->Add(ResultNode(result.peaks, imageNode, "_PrincipalDirection"), imageNode);
      if (m_Controls->m_OutputNumDirectionsBox->isChecked())
        this->GetDataStorage()->Add(ResultNode(result.numberOfPeaks, imageNode, "_NumDirections"), imageNode);
      return;
    }

    auto shImage = dynamic_cast<mitk::ShImage *>(imageNode->GetData());
    if (shImage == nullptr)
      mitkThrow() << "Select a tensor or an SH image.";

    // The entries of m_NormalizationBox: none, by the largest peak, unit vectors.
    constexpr mitk::PeakNormalization normalizations[] = {
      mitk::PeakNormalization::None, mitk::PeakNormalization::MaxVector, mitk::PeakNormalization::SingleVector};
    mitk::PeakParameters parameters;
    parameters.maxPeaks = m_Controls->m_MaxNumPeaksBox->value();
    parameters.relativeThreshold = m_Controls->m_PeakThresholdBox->value();
    parameters.absoluteThreshold = m_Controls->m_AbsoluteThresholdBox->value();
    parameters.angularThreshold = m_Controls->m_AngularThreshold->value();
    parameters.scaleByGfa = m_Controls->m_ScaleByGfaBox->isChecked();
    parameters.normalization = normalizations[std::clamp(m_Controls->m_NormalizationBox->currentIndex(), 0, 2)];

    const auto result = mitk::ExtractPeaks(shImage, parameters, mask);
    this->GetDataStorage()->Add(ResultNode(result.peaks, imageNode, "_PEAKS"), imageNode);
    if (m_Controls->m_OutputNumDirectionsBox->isChecked())
      this->GetDataStorage()->Add(ResultNode(result.numberOfPeaks, imageNode, "_NUM_DIRECTIONS"), imageNode);
  });
}

void QmitkOdfMaximaExtractionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer> &)
{
  this->OnImageSelectionChanged();
}

void QmitkOdfMaximaExtractionView::OnImageSelectionChanged()
{
  // The combo box admits tensor and SH images only.
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  m_Controls->m_StartPeakExtractionButton->setVisible(node.IsNotNull());

  // For tensors only the threshold (on the FA) and the normalization act; the other parameters
  // belong to the ODF maxima.
  const bool tensors = node.IsNotNull() && dynamic_cast<mitk::TensorImage *>(node->GetData()) != nullptr;
  m_Controls->label_7->setText(tensors ? "FA Threshold:" : "Absolute Threshold:");
  for (QWidget *odfOnly : std::initializer_list<QWidget *>{m_Controls->m_PeakThresholdBox, m_Controls->label_6,
                                                           m_Controls->m_MaxNumPeaksBox, m_Controls->label_5,
                                                           m_Controls->m_AngularThreshold, m_Controls->label_11,
                                                           m_Controls->m_ScaleByGfaBox, m_Controls->label_12})
    odfOnly->setEnabled(!tensors);
}
