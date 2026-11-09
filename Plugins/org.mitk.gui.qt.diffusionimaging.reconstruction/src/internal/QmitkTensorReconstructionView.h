/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkTensorReconstructionView_h
#define QmitkTensorReconstructionView_h

#include <QmitkAbstractView.h>
#include <mitkILifecycleAwarePart.h>
#include <mitkImage.h>
#include <mitkWeakPointer.h>

#include <itkImage.h>

#include <string>

#include "ui_QmitkTensorReconstructionViewControls.h"

/**
 * \brief Tensor reconstruction from a DWI, ODFs and DWIs from tensors, and a preview of the
 * voxels the b0 threshold excludes.
 */
class QmitkTensorReconstructionView : public QmitkAbstractView, public mitk::ILifecycleAwarePart
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;

  QmitkTensorReconstructionView() = default;
  ~QmitkTensorReconstructionView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

  void Activated() override {}
  /** Removes the threshold previews and releases the b0 image they were computed from. */
  void Deactivated() override;
  void Visible() override {}
  void Hidden() override {}

protected slots:

  void UpdateGui();
  void TensorsToOdf();
  void TensorsToDwi();
  void Advanced1CheckboxClicked();
  void Reconstruct();
  /** Shows the voxels whose mean b0 is below the threshold, which get no tensor. */
  void PreviewThreshold(int threshold);

protected:
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer> &nodes) override;

  Ui::QmitkTensorReconstructionViewControls *m_Controls = nullptr;

private:
  void CreateConnections();
  /** QmitkReconstructionResults::RunInBackground for this view. */
  template <typename Work, typename Done>
  void RunInBackground(const char *title, Work work, Done done);

  /** An action runs in the background. */
  bool m_Busy = false;

  /** Mean b0 image of the DWI the threshold preview shows, kept while that DWI is unchanged; each
   * step of the threshold box only thresholds it again. */
  itk::Image<short, 3>::Pointer m_PreviewB0;
  mitk::WeakPointer<const mitk::Image> m_PreviewDwi;
  itk::ModifiedTimeType m_PreviewDwiTime = 0;
};

#endif
