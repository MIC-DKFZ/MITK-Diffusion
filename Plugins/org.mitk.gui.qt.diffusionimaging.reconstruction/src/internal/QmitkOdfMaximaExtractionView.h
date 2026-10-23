/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkOdfMaximaExtractionView_h
#define QmitkOdfMaximaExtractionView_h

#include <QmitkAbstractView.h>

#include "ui_QmitkOdfMaximaExtractionViewControls.h"

/**
 * \brief Peaks of the ODFs of an SH image (local maxima on the sampled ODF, refined by gradient
 * descent) or principal eigenvectors of a tensor image.
 */
class QmitkOdfMaximaExtractionView : public QmitkAbstractView
{
  Q_OBJECT

public:
  QmitkOdfMaximaExtractionView() = default;
  ~QmitkOdfMaximaExtractionView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override {}

protected slots:

  void StartPeakExtraction();
  void OnImageSelectionChanged();

protected:
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer> &nodes) override;

  Ui::QmitkOdfMaximaExtractionViewControls *m_Controls = nullptr;
};

#endif
