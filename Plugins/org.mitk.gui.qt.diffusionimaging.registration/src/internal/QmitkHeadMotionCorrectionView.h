/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkHeadMotionCorrectionView_h
#define QmitkHeadMotionCorrectionView_h

#include <QmitkAbstractView.h>

#include "ui_QmitkHeadMotionCorrectionViewControls.h"

#include <mitkDataNode.h>
#include <mitkImage.h>

#include <QFutureWatcher>

#include <memory>
#include <string>

/**
 * \brief View for head-motion correction of DWIs (mitk::DWIHeadMotionCorrectionFilter).
 *
 * The correction runs in a worker thread, so the application stays responsive; one run at a time.
 */
class QmitkHeadMotionCorrectionView : public QmitkAbstractView
{
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkHeadMotionCorrectionView();
  ~QmitkHeadMotionCorrectionView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

protected slots:

  void StartCorrection();
  void OnCorrectionFinished();
  void UpdateGUI();             ///< update button activity depending on the selected image

private:

  /** The corrected image, or the reason why there is none. */
  struct Result
  {
    mitk::Image::Pointer image;
    std::string error;
  };

  std::unique_ptr<Ui::QmitkHeadMotionCorrectionViewControls> m_Controls;
  QFutureWatcher<Result> m_Watcher;
  mitk::DataNode::Pointer m_ImageNode;  ///< input of the running correction
};

#endif
