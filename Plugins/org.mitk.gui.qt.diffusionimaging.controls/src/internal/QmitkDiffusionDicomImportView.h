/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkDiffusionDicomImportView_h
#define QmitkDiffusionDicomImportView_h

#include <QmitkAbstractView.h>
#include "ui_QmitkDiffusionDicomImportViewControls.h"

/**
 * \brief Imports diffusion-weighted DICOM series from one or more folders (optionally their
 * subfolders), with mosaic splitting and averaging of duplicate gradients.
 */
class QmitkDiffusionDicomImportView : public QmitkAbstractView
{
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkDiffusionDicomImportView();
  ~QmitkDiffusionDicomImportView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

protected slots:

  void DicomLoadAddFolderNames();
  void DicomLoadDeleteFolderNames();
  void NewDicomLoadStartLoad() ;
  void AverageClicked();
  void Remove();
  void RecursiveSettingsChanged();
  void SetPrefixButtonPushed();
  void ResetPrefixButtonPushed();

protected:

  /** Shows the message in the status bar and logs it. */
  void Status(const QString& status);
  void Error(const QString& status);

  Ui::QmitkDiffusionDicomImportControls *m_Controls;

  QWidget *m_Parent;

};

#endif
