/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _QMITKDIFFUSIONQUANTIFICATIONVIEW_H_INCLUDED
#define _QMITKDIFFUSIONQUANTIFICATIONVIEW_H_INCLUDED

#include <QmitkAbstractView.h>
#include <string>
#include <itkVectorImage.h>

#include "ui_QmitkDiffusionQuantificationViewControls.h"
/*!
 * \ingroup org_mitk_gui_qt_diffusionquantification_internal
 *
 * \brief QmitkDiffusionQuantificationView
 *
 * Document your class here.
 */
class QmitkDiffusionQuantificationView : public QmitkAbstractView
{

  friend struct DqSelListener;

  // this is needed for all Qt objects that should have a MOC object (everything that derives from QObject)
  Q_OBJECT

  public:

  static const std::string VIEW_ID;

  QmitkDiffusionQuantificationView();
  virtual ~QmitkDiffusionQuantificationView();

  typedef itk::VectorImage< short, 3 >                                    ItkDwiType;

  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override;

protected slots:

  void GFA();
  void FA();
  void RA();
  void AD();
  void RD();
  void ClusterAnisotropy();
  void MD();

  void ADC_DWI();
  void MD_DWI();

  void OdfQuantify(int method);
  void OdfQuantification(int method) ;

  void TensorQuantify(int method);
  void TensorQuantification(int method) ;

  void UpdateGui();

protected:

  void DoAdcCalculation(bool fit);

  /// \brief called by QmitkAbstractView when DataManager's selection has changed
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  Ui::QmitkDiffusionQuantificationViewControls* m_Controls;
  static const float m_ScaleDAIValues;
};




#endif // _QMITKDIFFUSIONQUANTIFICATIONVIEW_H_INCLUDED

