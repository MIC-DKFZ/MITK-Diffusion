/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkDiffusionQuantificationView_h
#define QmitkDiffusionQuantificationView_h

#include <QmitkAbstractView.h>

#include <mitkDiffusionIndices.h>

#include <string>

#include "ui_QmitkDiffusionQuantificationViewControls.h"

/**
 * \brief Scalar indices of tensor images (FA, RA, AD, RD, CA, MD), ODF and SH images (GFA) and
 * DWIs (ADC, MD).
 */
class QmitkDiffusionQuantificationView : public QmitkAbstractView
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;

  QmitkDiffusionQuantificationView() = default;
  ~QmitkDiffusionQuantificationView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

protected slots:

  void UpdateGui();

protected:
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer> &nodes) override;

  Ui::QmitkDiffusionQuantificationViewControls *m_Controls = nullptr;

private:
  void ComputeGfa();
  void ComputeTensorIndex(mitk::TensorIndex index);
  void ComputeAdc(bool fit);
};

#endif
