/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkControlVisualizationPropertiesView_h
#define QmitkControlVisualizationPropertiesView_h

#include <QmitkAbstractView.h>
#include <QmitkSliceNavigationListener.h>
#include <mitkIRenderWindowPartListener.h>

#include "ui_QmitkControlVisualizationPropertiesViewControls.h"

#include <memory>
#include <string>
#include <vector>

class QActionGroup;
class QSlider;

/**
 * \brief Visualization controls for glyph images (ODF, tensor, SH), peak images and fiber
 * bundles, and the thick-slice menu. An edit applies to all selected nodes of the matching kind;
 * the controls show the values of the last selected node of each kind.
 */
class QmitkControlVisualizationPropertiesView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
  // this is needed for all Qt objects that should have a MOC object (everything that derives from QObject)
  Q_OBJECT

  public:

  static const std::string VIEW_ID;

  QmitkControlVisualizationPropertiesView();
  ~QmitkControlVisualizationPropertiesView() override;

  void CreateQtPartControl(QWidget *parent) override;

  void RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart) override;
  void RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart) override;

protected slots:

  void VisibleOdfsON_S();
  void VisibleOdfsON_T();
  void VisibleOdfsON_C();
  void ShowMaxNumberChanged();
  void NormalizationDropdownChanged(int);
  void ScalingFactorChanged(double);
  void AdditionalScaling(int);
  void OnThickSlicesModeSelected( QAction* action );
  void OnTSNumChanged(int num);
  void ResetColoring();
  void FiberSlicingThickness2D();
  void LineWidthChanged();
  void TubeRadiusChanged();
  void RibbonWidthChanged();
  void Toggle3DClipping();
  void FlipPeaks();
  void OnSliceChanged();
  void SetColor();
  void Toggle3DPeaks();

  /** Switches the colourisation mode of the glyphs. */
  void OnColourisationModeChanged();

protected:

  /** Recolours the fibers or peaks of one node when its colour property changes; stops
   *  observing when destroyed. */
  class ColorObserver
  {
  public:
    explicit ColorObserver(mitk::DataNode* node);
    ~ColorObserver();
    ColorObserver(const ColorObserver&) = delete;
    ColorObserver& operator=(const ColorObserver&) = delete;

  private:
    void ColorChanged();

    mitk::DataNode::Pointer m_Node;
    mitk::ColorProperty::Pointer m_Color;
    unsigned long m_Tag = 0;
  };

  void SetFocus() override;

  void NodeRemoved(const mitk::DataNode* node) override;

  /** called by QmitkAbstractView when DataManager's selection has changed */
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  void NodeAdded(const mitk::DataNode *node) override;
  void PlanarFigureFocus(mitk::DataNode* node);

  void Set3DClippingPlane(mitk::DataNode *node, int planeId);
  void SetTs(int m_CurrentThickSlicesMode, int num, mitk::AnatomicalPlane plane);
  void CreateThickSliceMenu();
  void UpdateThickSliceMenu(const QList<mitk::DataNode::Pointer>& nodes);

  /** Shows the 3D clipping plane id (0: none, 1-3: axial, sagittal, coronal) in the radio buttons. */
  void ShowClippingPlaneId(const mitk::DataNode* node);

  /** Calls action for every selected node whose data passes isKind. */
  template <class TPredicate, class TAction>
  void ForSelected(TPredicate isKind, TAction action)
  {
    for (const auto& node : m_SelectedNodes)
      if (node.IsNotNull() && node->GetData() != nullptr && isKind(node->GetData()))
        action(node.GetPointer());
  }

  Ui::QmitkControlVisualizationPropertiesViewControls* m_Controls;

  mitk::DataNode::Pointer m_NodeUsedForOdfVisualization;

  bool m_GlyIsOn_T;
  bool m_GlyIsOn_C;
  bool m_GlyIsOn_S;

  int m_CurrentThickSlicesMode;
  int m_CurrentThickSlicesNum;
  QLabel* m_TSLabel;
  QMenu* m_MyMenu;
  QSlider* m_TSSlider;
  QActionGroup* m_ThickSliceModeGroup;

  std::vector< mitk::DataNode::Pointer > m_SelectedNodes;
  std::vector< std::unique_ptr<ColorObserver> > m_ColorObservers;

  QmitkSliceNavigationListener  m_SliceChangeListener;
};

#endif
