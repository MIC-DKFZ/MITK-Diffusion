/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkPreprocessingView_h
#define QmitkPreprocessingView_h

#include <QmitkAbstractView.h>
#include <mitkILifecycleAwarePart.h>

#include "ui_QmitkPreprocessingViewControls.h"

#include <mitkDataNode.h>
#include <mitkImage.h>

#include <functional>
#include <memory>

/**
 * \brief Views and edits DWIs and other images: gradients, volumes, header, resampling, cropping.
 *
 * The edits themselves are mitk::DiffusionImageEditing functions; each result becomes a new node
 * below the edited one.
 */
class QmitkPreprocessingView : public QmitkAbstractView, public mitk::ILifecycleAwarePart
{
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkPreprocessingView();
  ~QmitkPreprocessingView() override;

  void CreateQtPartControl(QWidget *parent) override;
  void SetFocus() override;

  void Activated() override;
  void Deactivated() override;
  void Visible() override;
  void Hidden() override;

protected slots:

  void AverageGradients();
  void ExtractB0();
  void MergeDwis();
  void DoApplyHeader();
  void DoReduceGradientDirections();
  void DoShowGradientDirections();
  void DoHalfSphereGradientDirections();
  void UpdateDwiBValueMapRounder(int i);
  void DoLengthCorrection();
  void DoResampleImage();
  void DoCropImage();
  void DoUpdateInterpolationGui(int i);
  void DoRemoveGradient();
  void DoExtractGradient();
  void DoFlipAxis();
  void OnImageSelectionChanged();
  void DoFlipGradientDirections();
  void DoClearRotationOfGradients();
  void DoAlignImages();
  void DoSwapGradientDimensions();

private:

  void CreateConnections();

  /** The node chosen in the view's image box and its image, or nullptr. */
  mitk::DataNode::Pointer SelectedNode() const;
  mitk::Image::Pointer SelectedImage() const;

  /**
   * Runs an edit of the selected image, reports exceptions in a message box, and adds the result
   * as a new node named after the selected node plus suffix below it.
   */
  void ApplyEdit(const QString& description, const QString& suffix, const std::function<mitk::Image::Pointer(const mitk::Image*)>& edit, bool reinitializeViews);

  void AddResult(mitk::DataNode* parent, mitk::BaseData* data, const QString& name, bool reinitializeViews);

  void UpdateBValueTableWidget();
  void UpdateGradientDetails();
  void CleanBValueTableWidget();
  void ClearHeaderTables();

  std::unique_ptr<Ui::QmitkPreprocessingViewControls> m_Controls;

  /** The node the controls show; Visible() refreshes them only for another node. */
  mitk::DataNode::Pointer m_ShownNode;
};

#endif
