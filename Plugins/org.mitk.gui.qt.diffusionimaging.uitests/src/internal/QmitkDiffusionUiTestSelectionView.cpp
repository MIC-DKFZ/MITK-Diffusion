/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestSelectionView.h"

const QString QmitkDiffusionUiTestSelectionView::VIEW_ID = "org.mitk.views.diffusionuitestselection";

void QmitkDiffusionUiTestSelectionView::Select(const QList<mitk::DataNode::Pointer> &nodes)
{
  this->FireNodesSelected(nodes);
}

void QmitkDiffusionUiTestSelectionView::CreateQtPartControl(QWidget *)
{
}

void QmitkDiffusionUiTestSelectionView::SetFocus()
{
}
