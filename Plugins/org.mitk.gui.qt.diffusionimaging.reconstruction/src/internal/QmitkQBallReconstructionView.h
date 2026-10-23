/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkQBallReconstructionView_h
#define QmitkQBallReconstructionView_h

#include <QmitkAbstractView.h>

#include <string>
#include <vector>

#include "ui_QmitkQBallReconstructionViewControls.h"

class QCheckBox;

/**
 * \brief Q-ball reconstruction (numerical, analytical, multi-shell) of a DWI, and conversion of
 * SH images to sampled ODF images.
 */
class QmitkQBallReconstructionView : public QmitkAbstractView
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;

  QmitkQBallReconstructionView() = default;
  ~QmitkQBallReconstructionView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override {}

protected slots:

  void UpdateGui();
  void Reconstruct();
  void ConvertShImage();
  void MethodChosen(int method);

protected:
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer> &nodes) override;

  Ui::QmitkQBallReconstructionViewControls *m_Controls = nullptr;

private:
  void CreateConnections();
  void GenerateShellSelectionUI(mitk::DataNode::Pointer node);
  /** QmitkReconstructionResults::RunInBackground for this view. */
  template <typename Work, typename Done>
  void RunInBackground(const char *title, Work work, Done done);

  /** An action runs in the background. */
  bool m_Busy = false;

  /** Check boxes of the shells of the selected DWI, with their b-values */
  std::vector<QCheckBox *> m_ShellCheckBoxes;
  std::vector<unsigned int> m_ShellSelectionBValues;
  /** The node the check boxes were built for; compared, never dereferenced */
  const mitk::DataNode *m_ShellSelectionNode = nullptr;
};

#endif
