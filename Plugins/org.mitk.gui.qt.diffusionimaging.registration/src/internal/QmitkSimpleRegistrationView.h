/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkSimpleRegistrationView_h
#define QmitkSimpleRegistrationView_h

#include <QmitkAbstractView.h>
#include "ui_QmitkSimpleRegistrationViewControls.h"

#include <mitkImage.h>
#include <QmitkRegistrationJob.h>

#include <memory>

/**
 * \brief View for rigid and affine registration of images and DWIs, and for applying registrations
 * to tractograms.
 *
 * Of a DWI one volume is registered; the result is applied to all volumes and turns the gradients
 * (mitk::DiffusionImageMapping).
 */
class QmitkSimpleRegistrationView : public QmitkAbstractView
{
  Q_OBJECT

public:

  using ItkFloatImageType = itk::Image<float, 3>;
  static const std::string VIEW_ID;

  QmitkSimpleRegistrationView();
  ~QmitkSimpleRegistrationView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

protected slots:

  void MovingImageChanged();
  void FixedImageChanged();
  void TractoChanged();
  void StartRegistration();
  void StartTractoRegistration();
  void OnRegResultIsAvailable(mitk::MAPRegistrationWrapper::Pointer spResultRegistration, const QmitkRegistrationJob* pRegJob);
  void OnRegistrationError(QString message);
  void OnRegistrationFinished();

protected:

  /// \brief called by QmitkAbstractView when DataManager's selection has changed
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  /** Enables "Start" if both images can be registered and no registration runs. */
  void UpdateStartButton();

  std::unique_ptr<Ui::QmitkSimpleRegistrationViewControls> m_Controls;
  /** Moving image and type of the running registration; only one runs at a time. */
  mitk::DataNode::Pointer  m_MovingImageNode;
  int  m_RegistrationType = 0;
  bool m_RegistrationRunning = false;
};

#endif
