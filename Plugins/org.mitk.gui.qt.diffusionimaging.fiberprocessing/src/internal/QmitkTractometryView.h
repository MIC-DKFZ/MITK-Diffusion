/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkTractometryView_h
#define QmitkTractometryView_h

#include <QmitkAbstractView.h>
#include "ui_QmitkTractometryViewControls.h"
#include <itkImage.h>
#include <mitkImage.h>
#include <mitkPixelType.h>
#include <mitkFiberBundle.h>
#include <mitkILifecycleAwarePart.h>
#include <QmitkChartWidget.h>

/*!
\brief Weight fibers by linearly fitting them to the image data.

*/
class QmitkTractometryView : public QmitkAbstractView, public mitk::ILifecycleAwarePart
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkTractometryView();
  virtual ~QmitkTractometryView();

  virtual void CreateQtPartControl(QWidget *parent) override;

  void StaticResamplingTractometry(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector< std::vector< double > >& data, std::string& clipboard_string);

  void NearestCentroidPointTractometry(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector< std::vector< double > >& data, std::string& clipboard_string);

  void AlongTractRadiomicsPreprocessing(mitk::Image::Pointer image, mitk::DataNode::Pointer node, std::vector<std::vector<double> > &data, std::string &clipboard_string);

  virtual void SetFocus() override;

  virtual void Activated() override;
  virtual void Deactivated() override;
  virtual void Visible() override;
  virtual void Hidden() override;

protected slots:

  void UpdateGui();
  void StartTractometry();

protected:

  /// \brief called by QmitkAbstractView when DataManager's selection has changed
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  Ui::QmitkTractometryViewControls* m_Controls;

  std::string RGBToHexString(double *rgb);

  mitk::FiberBundle::Pointer m_ReferenceFib;
  QList<mitk::DataNode::Pointer> m_CurrentSelection;
  unsigned int m_NumSamplingPoints;
  bool  m_Visible;
};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

