/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionUiTestSelectionView_h
#define QmitkDiffusionUiTestSelectionView_h

#include <QmitkAbstractView.h>

/**
 * \brief Selects data nodes for the diffusion UI tests.
 *
 * The diffusion views take most of their input from the workbench selection, which the Data
 * Manager sets. This view sets it the same way (BlueBerry selection service), without depending
 * on the Data Manager's widgets.
 */
class QmitkDiffusionUiTestSelectionView : public QmitkAbstractView
{
  Q_OBJECT

public:
  static const QString VIEW_ID;

  void Select(const QList<mitk::DataNode::Pointer> &nodes);

protected:
  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;
};

#endif
