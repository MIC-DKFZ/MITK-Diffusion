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

#include <QmitkDiffusionViewControlsGuard.h>

#include <mitkDiffusionMasks.h>
#include <mitkDiffusionPeaks.h>
#include <mitkExceptionMacro.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateOr.h>

#include <algorithm>
#include <initializer_list>

using QmitkReconstructionResults::AddResult;
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
  QmitkDiffusionViewControlsGuard::Install(parent, this);
  m_Controls->setupUi(parent);

  connect(m_Controls->m_StartPeakExtractionButton, &QAbstractButton::clicked, this, &QmitkOdfMaximaExtractionView::StartPeakExtraction);

  m_Controls->m_MaskBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_ImageBox->SetDataStorage(this->GetDataStorage());

  m_Controls->m_MaskBox->SetPredicate(mitk::GetBinaryMaskPredicate());
  m_Controls->m_ImageBox->SetPredicate(mitk::NodePredicateOr::New(mitk::TNodePredicateDataType<mitk::TensorImage>::New(),
                                                                  mitk::TNodePredicateDataType<mitk::ShImage>::New()));
  m_Controls->m_MaskBox->SetZeroEntryText("--");

  connect(m_Controls->m_ImageBox, &QmitkDataStorageComboBox::OnSelectionChanged, this, &QmitkOdfMaximaExtractionView::OnImageSelectionChanged);

  OnImageSelectionChanged();
}

template <typename Work, typename Done>
void QmitkOdfMaximaExtractionView::RunInBackground(const char *title, Work work, Done done)
{
  QmitkReconstructionResults::RunInBackground(this, m_Busy, [this]() { this->OnImageSelectionChanged(); }, title, std::move(work), std::move(done));
}

void QmitkOdfMaximaExtractionView::StartPeakExtraction()
{
  auto imageNode = m_Controls->m_ImageBox->GetSelectedNode();
  if (m_Busy || imageNode.IsNull())
    return;
  auto maskNode = m_Controls->m_MaskBox->GetSelectedNode();
  // A segmentation becomes a mask image in the worker thread, where its errors are reported too.
  mitk::BaseData::ConstPointer maskData = maskNode.IsNotNull() ? maskNode->GetData() : nullptr;
  const bool outputNumberOfPeaks = m_Controls->m_OutputNumDirectionsBox->isChecked();
  const auto done = [this, imageNode, outputNumberOfPeaks](const char *peaksSuffix, const char *numberSuffix) {
    return [this, imageNode, outputNumberOfPeaks, peaksSuffix, numberSuffix](const mitk::PeakResult &result) {
      AddResult(this->GetDataStorage(), ResultNode(result.peaks, imageNode, peaksSuffix), imageNode);
      if (outputNumberOfPeaks)
        AddResult(this->GetDataStorage(), ResultNode(result.numberOfPeaks, imageNode, numberSuffix), imageNode);
    };
  };

  if (mitk::TensorImage::ConstPointer tensors = dynamic_cast<mitk::TensorImage *>(imageNode->GetData()))
  {
    // The absolute threshold acts as FA threshold; without normalization the vectors are scaled by FA.
    const double faThreshold = m_Controls->m_AbsoluteThresholdBox->value();
    const bool normalize = m_Controls->m_NormalizationBox->currentIndex() != 0;
    this->RunInBackground(
      "Peak extraction not possible:",
      [tensors, faThreshold, normalize, maskData]() {
        return mitk::ExtractTensorPeaks(tensors, faThreshold, normalize, mitk::GetMaskImage(maskData));
      },
      done("_PrincipalDirection", "_NumDirections"));
    return;
  }

  mitk::ShImage::ConstPointer shImage = dynamic_cast<mitk::ShImage *>(imageNode->GetData());
  if (shImage.IsNull())
  {
    QmitkReconstructionResults::Run("Peak extraction not possible:", []() { mitkThrow() << "Select a tensor or an SH image."; });
    return;
  }

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

  this->RunInBackground(
    "Peak extraction not possible:",
    [shImage, parameters, maskData]() { return mitk::ExtractPeaks(shImage, parameters, mitk::GetMaskImage(maskData)); },
    done("_PEAKS", "_NUM_DIRECTIONS"));
}

void QmitkOdfMaximaExtractionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer> &)
{
  this->OnImageSelectionChanged();
}

void QmitkOdfMaximaExtractionView::OnImageSelectionChanged()
{
  // The combo box admits tensor and SH images only. One extraction at a time: the button stays
  // disabled while one runs in the background.
  auto node = m_Controls->m_ImageBox->GetSelectedNode();
  m_Controls->m_StartPeakExtractionButton->setVisible(node.IsNotNull());
  m_Controls->m_StartPeakExtractionButton->setEnabled(!m_Busy);

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
