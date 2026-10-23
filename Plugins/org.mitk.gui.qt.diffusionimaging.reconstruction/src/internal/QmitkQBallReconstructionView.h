/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _QMITKQBALLRECONSTRUCTIONVIEW_H_INCLUDED
#define _QMITKQBALLRECONSTRUCTIONVIEW_H_INCLUDED

#include <QmitkAbstractView.h>

#include <string>
#include <vector>

#include "ui_QmitkQBallReconstructionViewControls.h"

#include <mitkImage.h>
#include <mitkDiffusionPropertyHelper.h>
#include <itkVectorImage.h>
#include <mitkShImage.h>

typedef short DiffusionPixelType;

class QCheckBox;

/*!
 * \ingroup org_mitk_gui_qt_qballreconstruction_internal
 *
 * \brief QmitkQBallReconstructionView
 *
 * Document your class here.
 */
class QmitkQBallReconstructionView : public QmitkAbstractView
{

  typedef mitk::DiffusionPropertyHelper::GradientDirectionType            GradientDirectionType;
  typedef mitk::DiffusionPropertyHelper::GradientDirectionsContainerType  GradientDirectionContainerType;
  typedef mitk::DiffusionPropertyHelper::BValueMapType                    BValueMapType;
  typedef itk::VectorImage< DiffusionPixelType, 3 >                       ITKDiffusionImageType;

  // this is needed for all Qt objects that should have a MOC object (everything that derives from QObject)
  Q_OBJECT

  public:

  static const std::string VIEW_ID;

  QmitkQBallReconstructionView();
  virtual ~QmitkQBallReconstructionView();

  virtual void CreateQtPartControl(QWidget *parent) override;

  /// \brief Creation of the connections of main and control widget
  virtual void CreateConnections();

  ///
  /// Sets the focus to an internal widget.
  ///
  virtual void SetFocus() override;

protected slots:

  void UpdateGui();
  void ReconstructStandard();
  void ConvertShImage();
  void MethodChoosen(int method);
  void Reconstruct(int method, int normalization);

  void NumericalQBallReconstruction(mitk::DataNode::Pointer node);
  void AnalyticalQBallReconstruction(mitk::DataNode::Pointer node, int normalization);
  void MultiQBallReconstruction(mitk::DataNode::Pointer node);

protected:

  /// \brief called by QmitkAbstractView when DataManager's selection has changed
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  Ui::QmitkQBallReconstructionViewControls* m_Controls;

  template<int L>
  void TemplatedAnalyticalQBallReconstruction(mitk::DataNode* dataNodePointer, double lambda, int normalization);

  template<int L>
  void TemplatedMultiQBallReconstruction(double lambda, mitk::DataNode*);

private:

  /** Check boxes of the shells of the selected DWI, with their b-values */
  std::vector<QCheckBox*> m_ShellCheckBoxes;
  std::vector<unsigned int> m_ShellSelectionBValues;
  /** The node the check boxes were built for; compared, never dereferenced */
  const mitk::DataNode* m_ShellSelectionNode;
  void GenerateShellSelectionUI(mitk::DataNode::Pointer node);
};




#endif // _QMITKQBALLRECONSTRUCTIONVIEW_H_INCLUDED

