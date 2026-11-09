/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkDenoisingView_h
#define QmitkDenoisingView_h

#include <QmitkAbstractView.h>

#include "ui_QmitkDenoisingViewControls.h"

#include <mitkDataNode.h>
#include <mitkImage.h>

#include <QFutureWatcher>

#include <memory>
#include <string>

/**
 * \class QmitkDenoisingView
 * \brief View for denoising diffusion-weighted images.
 *
 * The denoising (mitk::DiffusionImageDenoising) runs in a worker thread, so the application stays
 * responsive; one run at a time.
 */
class QmitkDenoisingView : public QmitkAbstractView
{
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkDenoisingView();
  ~QmitkDenoisingView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

private slots:

  void StartDenoising();
  void OnDenoisingFinished();
  void UpdateGui();

private:

  /** The denoised image, or the reason why there is none. */
  struct Result
  {
    mitk::Image::Pointer image;
    std::string error;
  };

  std::unique_ptr<Ui::QmitkDenoisingViewControls> m_Controls;
  QFutureWatcher<Result> m_Watcher;
  mitk::DataNode::Pointer m_ImageNode;  ///< input of the running denoising
  QString m_ResultName;
};

#endif
