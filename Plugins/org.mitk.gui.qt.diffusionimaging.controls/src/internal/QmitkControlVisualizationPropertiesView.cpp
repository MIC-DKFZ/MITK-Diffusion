/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkControlVisualizationPropertiesView.h"

#include <QmitkRenderWindow.h>
#include <itkCommand.h>
#include <itkFlipPeaksFilter.h>
#include <mitkClippingProperty.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkFiberBundle.h>
#include <mitkImage.h>
#include <mitkImageToItk.h>
#include <mitkImageWriteAccessor.h>
#include <mitkOdfImage.h>
#include <mitkOdfNormalizationMethodProperty.h>
#include <mitkOdfScaleByProperty.h>
#include <mitkPeakImage.h>
#include <mitkPlanarFigure.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkPlaneGeometry.h>
#include <mitkRenderingManager.h>
#include <mitkResliceMethodProperty.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkWorkbenchUtil.h>
#include <usModuleRegistry.h>

#include <QActionGroup>
#include <QColorDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QMenu>
#include <QSignalBlocker>
#include <QSlider>
#include <QWidgetAction>

#include <cstring>

const std::string QmitkControlVisualizationPropertiesView::VIEW_ID = "org.mitk.views.controlvisualizationpropertiesview";

namespace
{
  /** Renderer of the 2D render window that shows the given plane, or nullptr: no editor may be
   *  open, and MITK's MxN editor need not offer every plane. */
  mitk::BaseRenderer* GetRendererOf(mitk::IRenderWindowPart* renderWindowPart, mitk::AnatomicalPlane plane)
  {
    if (renderWindowPart == nullptr)
      return nullptr;
    QmitkRenderWindow* renderWindow = renderWindowPart->GetQmitkRenderWindow(plane);
    return renderWindow != nullptr ? renderWindow->GetRenderer() : nullptr;
  }

  /** ODF, tensor and SH images: the data the ODF mapper draws glyphs for. */
  bool IsGlyphImage(const mitk::BaseData* data)
  {
    return dynamic_cast<const mitk::OdfImage*>(data) || dynamic_cast<const mitk::TensorImage*>(data) || dynamic_cast<const mitk::ShImage*>(data);
  }

  bool IsFiberBundle(const mitk::BaseData* data)
  {
    return dynamic_cast<const mitk::FiberBundle*>(data) != nullptr;
  }

  bool IsPeakImage(const mitk::BaseData* data)
  {
    return dynamic_cast<const mitk::PeakImage*>(data) != nullptr;
  }

  bool IsFiberBundleOrPeakImage(const mitk::BaseData* data)
  {
    return IsFiberBundle(data) || IsPeakImage(data);
  }

  /** Colours the fibers or peaks of the node with the node's colour. */
  void ApplyNodeColor(mitk::DataNode* node)
  {
    float color[3];
    node->GetColor(color);
    if (auto fib = dynamic_cast<mitk::FiberBundle*>(node->GetData()))
      fib->SetFiberColors(color[0]*255, color[1]*255, color[2]*255);
    else if (auto peaks = dynamic_cast<mitk::PeakImage*>(node->GetData()))
      peaks->SetCustomColor(color[0]*255, color[1]*255, color[2]*255);
  }

  /** The 2D plane of a 3D clipping plane id (1: axial, 2: sagittal, 3: coronal). */
  mitk::AnatomicalPlane ClippingPlaneOf(int planeId)
  {
    switch (planeId)
    {
    case 2:
      return mitk::AnatomicalPlane::Sagittal;
    case 3:
      return mitk::AnatomicalPlane::Coronal;
    default:
      return mitk::AnatomicalPlane::Axial;
    }
  }
}

QmitkControlVisualizationPropertiesView::ColorObserver::ColorObserver(mitk::DataNode* node)
  : m_Node(node)
  , m_Color(dynamic_cast<mitk::ColorProperty*>(node->GetProperty("color", nullptr)))
{
  if (m_Color.IsNull())
    return;
  auto command = itk::SimpleMemberCommand<ColorObserver>::New();
  command->SetCallbackFunction(this, &ColorObserver::ColorChanged);
  m_Tag = m_Color->AddObserver(itk::ModifiedEvent(), command);
}

QmitkControlVisualizationPropertiesView::ColorObserver::~ColorObserver()
{
  if (m_Color.IsNotNull())
    m_Color->RemoveObserver(m_Tag);
}

void QmitkControlVisualizationPropertiesView::ColorObserver::ColorChanged()
{
  ApplyNodeColor(m_Node);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

QmitkControlVisualizationPropertiesView::QmitkControlVisualizationPropertiesView()
  : QmitkAbstractView(),
    m_Controls(nullptr),
    m_GlyIsOn_T(false),
    m_GlyIsOn_C(false),
    m_GlyIsOn_S(false),
    m_CurrentThickSlicesMode(1),
    m_CurrentThickSlicesNum(0),
    m_TSLabel(nullptr),
    m_MyMenu(nullptr),
    m_TSSlider(nullptr),
    m_ThickSliceModeGroup(nullptr)
{
}


QmitkControlVisualizationPropertiesView::~QmitkControlVisualizationPropertiesView()
{
  // The observers call back into the nodes' data; the nodes may outlive this view.
  m_ColorObservers.clear();
  // The widgets belong to the parent; only the Ui struct is this view's.
  delete m_Controls;
}

void QmitkControlVisualizationPropertiesView::SetTs(int currentThickSlicesMode, int num, mitk::AnatomicalPlane plane)
{
  mitk::BaseRenderer::Pointer renderer = GetRendererOf(this->GetRenderWindowPart(mitk::WorkbenchUtil::IRenderWindowPartStrategy::OPEN), plane);
  if (renderer.IsNotNull() && renderer->GetCurrentWorldPlaneGeometryNode() != nullptr)
  {
    renderer->GetCurrentWorldPlaneGeometryNode()->SetProperty("reslice.thickslices.num", mitk::IntProperty::New(num));
    if(num>0)
    {
      renderer->GetCurrentWorldPlaneGeometryNode()->SetProperty("reslice.thickslices", mitk::ResliceMethodProperty::New(currentThickSlicesMode));
      renderer->GetCurrentWorldPlaneGeometryNode()->SetProperty("reslice.thickslices.showarea", mitk::BoolProperty::New(true));
    }
    else
    {
      renderer->GetCurrentWorldPlaneGeometryNode()->SetProperty("reslice.thickslices", mitk::ResliceMethodProperty::New(0));
      renderer->GetCurrentWorldPlaneGeometryNode()->SetProperty("reslice.thickslices.showarea", mitk::BoolProperty::New(false));
    }

    renderer->SendUpdateSlice();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}

void QmitkControlVisualizationPropertiesView::OnThickSlicesModeSelected( QAction* action )
{
  m_CurrentThickSlicesMode = action->data().toInt();

  switch( m_CurrentThickSlicesMode )
  {
  case 1:
    this->m_Controls->m_TSMenu->setText("MIP");
    break;

  case 2:
    this->m_Controls->m_TSMenu->setText("SUM");
    break;

  case 3:
    this->m_Controls->m_TSMenu->setText("WEIGH");
    break;

  default:
    return;
  }

  for (auto plane : {mitk::AnatomicalPlane::Axial, mitk::AnatomicalPlane::Sagittal, mitk::AnatomicalPlane::Coronal})
    SetTs(m_CurrentThickSlicesMode, m_CurrentThickSlicesNum, plane);
}


void QmitkControlVisualizationPropertiesView::OnTSNumChanged( int num )
{
  m_CurrentThickSlicesNum = num;
  for (auto plane : {mitk::AnatomicalPlane::Axial, mitk::AnatomicalPlane::Sagittal, mitk::AnatomicalPlane::Coronal})
    SetTs(m_CurrentThickSlicesMode, m_CurrentThickSlicesNum, plane);
  m_TSLabel->setText(QString::number( num*2 + 1 ));
}


void QmitkControlVisualizationPropertiesView::CreateQtPartControl(QWidget *parent)
{
  if (m_Controls)
    return;

  // create GUI widgets
  m_Controls = new Ui::QmitkControlVisualizationPropertiesViewControls;
  m_Controls->setupUi(parent);

  using Self = QmitkControlVisualizationPropertiesView;
  connect(m_Controls->m_VisibleOdfsON_T, &QToolButton::clicked, this, &Self::VisibleOdfsON_T);
  connect(m_Controls->m_VisibleOdfsON_S, &QToolButton::clicked, this, &Self::VisibleOdfsON_S);
  connect(m_Controls->m_VisibleOdfsON_C, &QToolButton::clicked, this, &Self::VisibleOdfsON_C);
  connect(m_Controls->m_ShowMaxNumber, &QSpinBox::editingFinished, this, &Self::ShowMaxNumberChanged);
  connect(m_Controls->m_NormalizationDropdown, &QComboBox::currentIndexChanged, this, &Self::NormalizationDropdownChanged);
  connect(m_Controls->m_ScalingFactor, &QDoubleSpinBox::valueChanged, this, &Self::ScalingFactorChanged);
  connect(m_Controls->m_AdditionalScaling, &QComboBox::currentIndexChanged, this, &Self::AdditionalScaling);
  connect(m_Controls->m_ResetColoring, &QToolButton::clicked, this, &Self::ResetColoring);
  connect(m_Controls->m_ResetColoring2, &QToolButton::clicked, this, &Self::ResetColoring);
  connect(m_Controls->m_FiberClippingBox, &QDoubleSpinBox::editingFinished, this, &Self::FiberSlicingThickness2D);
  connect(m_Controls->m_LineWidth, &QDoubleSpinBox::editingFinished, this, &Self::LineWidthChanged);
  connect(m_Controls->m_TubeWidth, &QDoubleSpinBox::editingFinished, this, &Self::TubeRadiusChanged);
  connect(m_Controls->m_RibbonWidth, &QDoubleSpinBox::editingFinished, this, &Self::RibbonWidthChanged);
  connect(m_Controls->m_OdfColorBox, &QComboBox::currentIndexChanged, this, &Self::OnColourisationModeChanged);
  // switching the radio buttons unchecks one and checks another; the checked one applies
  for (QRadioButton* clip : {m_Controls->m_Clip0, m_Controls->m_Clip1, m_Controls->m_Clip2, m_Controls->m_Clip3})
    connect(clip, &QRadioButton::toggled, this, [this](bool checked) { if (checked) this->Toggle3DClipping(); });
  connect(m_Controls->m_FlipClipBox, &QCheckBox::checkStateChanged, this, &Self::Toggle3DClipping);
  connect(m_Controls->m_Enable3dPeaks, &QCheckBox::checkStateChanged, this, &Self::Toggle3DPeaks);
  connect(m_Controls->m_FlipPeaksButton, &QCommandLinkButton::clicked, this, &Self::FlipPeaks);
  connect(m_Controls->m_SetColor1, &QToolButton::clicked, this, &Self::SetColor);
  connect(m_Controls->m_SetColor2, &QToolButton::clicked, this, &Self::SetColor);

  m_Controls->m_BundleControlsFrame->setVisible(false);
  m_Controls->m_ImageControlsFrame->setVisible(false);
  m_Controls->m_PeakImageFrame->setVisible(false);
  m_Controls->m_3DClippingBox->setVisible(false);
  m_Controls->m_NormalizationFrame->setVisible(false);

  m_MyMenu = new QMenu(parent);
  m_Controls->m_TSMenu->setMenu( m_MyMenu );
  this->CreateThickSliceMenu();

  // The listener follows the editor through RenderWindowPartActivated/Deactivated.
  connect(&m_SliceChangeListener, &QmitkSliceNavigationListener::SliceChanged, this, &Self::OnSliceChanged);
  if (mitk::IRenderWindowPart* renderWindowPart = this->GetRenderWindowPart())
    m_SliceChangeListener.RenderWindowPartActivated(renderWindowPart);
}

void QmitkControlVisualizationPropertiesView::RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart)
{
  m_SliceChangeListener.RenderWindowPartActivated(renderWindowPart);
}

void QmitkControlVisualizationPropertiesView::RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart)
{
  m_SliceChangeListener.RenderWindowPartDeactivated(renderWindowPart);
}

void QmitkControlVisualizationPropertiesView::SetColor()
{
  // One dialog for all selected tractograms and peak images (the buttons are in their frames).
  QColor c = QColorDialog::getColor();
  if (!c.isValid())
    return;

  float rgb[3];
  rgb[0] = static_cast<float>(c.redF());
  rgb[1] = static_cast<float>(c.greenF());
  rgb[2] = static_cast<float>(c.blueF());
  ForSelected(IsFiberBundleOrPeakImage, [&rgb](mitk::DataNode* node)
  {
    // A selected node's ColorObserver recolours it when its colour property changes. Recolouring a
    // tractogram of millions of points takes noticeable time, so it is done here only when setting
    // the colour changed no existing property: same colour (e.g. after colouring by orientation) or
    // a new property.
    const mitk::BaseProperty* colorBefore = node->GetProperty("color", nullptr, false);
    const auto modifiedBefore = colorBefore != nullptr ? colorBefore->GetMTime() : 0;
    node->SetColor(rgb);
    if (colorBefore == nullptr || node->GetProperty("color", nullptr, false) != colorBefore
        || colorBefore->GetMTime() == modifiedBefore)
      ApplyNodeColor(node);
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::SetFocus()
{
  m_Controls->m_TSMenu->setFocus();
}


void QmitkControlVisualizationPropertiesView::NodeRemoved(const mitk::DataNode* node)
{
  // Otherwise the removed ODF node and its image stay in memory until another one is selected.
  if (node == m_NodeUsedForOdfVisualization.GetPointer())
    m_NodeUsedForOdfVisualization = nullptr;
}


// set diffusion image channel to b0 volume
void QmitkControlVisualizationPropertiesView::NodeAdded(const mitk::DataNode *node)
{
  mitk::DataNode* notConst = const_cast<mitk::DataNode*>(node);

  bool isDiffusionImage( mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( dynamic_cast<mitk::Image *>(node->GetData())) );

  if (isDiffusionImage)
  {
    mitk::Image::Pointer dimg = dynamic_cast<mitk::Image*>(notConst->GetData());

    // Show the first b0 volume. MITK's image mapper shows the volume named by
    // "Image.Displayed Component" and sets it to 0 by default. Another value was chosen by the user
    // and restored from a scene file, whose properties MITK's scene reader sets before the node is
    // added; it is kept.
    int displayedComponent = 0;
    if (notConst->GetIntProperty("Image.Displayed Component", displayedComponent) && displayedComponent != 0)
      return;

    int displayChannelPropertyValue = 0;
    const auto& map = mitk::DiffusionPropertyHelper::GetBValueMap(dimg);
    const auto b0 = map.find(0);
    if (b0 != map.end() && !b0->second.empty())
      displayChannelPropertyValue = static_cast<int>(b0->second.front());

    notConst->SetIntProperty("Image.Displayed Component", displayChannelPropertyValue );
  }
}


void QmitkControlVisualizationPropertiesView::ShowClippingPlaneId(const mitk::DataNode* node)
{
  int clippingPlaneId = -1;
  node->GetPropertyValue("3DClippingPlaneId", clippingPlaneId);
  QRadioButton* buttons[] = {m_Controls->m_Clip0, m_Controls->m_Clip1, m_Controls->m_Clip2, m_Controls->m_Clip3};
  buttons[clippingPlaneId>=0 && clippingPlaneId<=3 ? clippingPlaneId : 0]->setChecked(true);
}

/* OnSelectionChanged is registered to SelectionService, therefore no need to
implement SelectionService Listener explicitly */
void QmitkControlVisualizationPropertiesView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& nodes)
{
  m_Controls->m_BundleControlsFrame->setVisible(false);
  m_Controls->m_ImageControlsFrame->setVisible(false);
  m_Controls->m_PeakImageFrame->setVisible(false);
  m_Controls->m_3DClippingBox->setVisible(false);
  m_Controls->m_FlipClipBox->setVisible(false);
  m_Controls->m_Enable3dPeaks->setVisible(false);

  m_Controls->m_NumberGlyphsFrame->setVisible(false);
  m_Controls->m_GlyphFrame->setVisible(false);
  m_Controls->m_TSMenu->setVisible(false);

  // Only the selected tractograms and peak images are recoloured when their colour changes.
  m_ColorObservers.clear();
  m_SelectedNodes.clear();

  // Showing a node's values in the controls must not write them into the other selected nodes
  // through the controls' change signals.
  const QSignalBlocker blockers[] = {
    QSignalBlocker(m_Controls->m_Clip0), QSignalBlocker(m_Controls->m_Clip1), QSignalBlocker(m_Controls->m_Clip2),
    QSignalBlocker(m_Controls->m_Clip3), QSignalBlocker(m_Controls->m_FlipClipBox),
    QSignalBlocker(m_Controls->m_Enable3dPeaks), QSignalBlocker(m_Controls->m_LineWidth),
    QSignalBlocker(m_Controls->m_TubeWidth), QSignalBlocker(m_Controls->m_RibbonWidth),
    QSignalBlocker(m_Controls->m_FiberClippingBox), QSignalBlocker(m_Controls->m_ShowMaxNumber),
    QSignalBlocker(m_Controls->m_NormalizationDropdown), QSignalBlocker(m_Controls->m_ScalingFactor),
    QSignalBlocker(m_Controls->m_AdditionalScaling), QSignalBlocker(m_Controls->m_OdfColorBox)};

  // As in MITK's measurement view, the last selected planar figure gets the focus.
  mitk::DataNode::Pointer lastPlanarFigure;
  for (mitk::DataNode::Pointer node: nodes)
  {
    if(node.IsNull())
      continue;

    mitk::BaseData* nodeData = node->GetData();
    if(nodeData == nullptr)
      continue;

    if (IsPeakImage(nodeData))
    {
      m_Controls->m_PeakImageFrame->setVisible(true);
      m_ColorObservers.push_back(std::make_unique<ColorObserver>(node));
      this->ShowClippingPlaneId(node);

      bool enable3dPeaks = false;
      node->GetBoolProperty("Enable3DPeaks", enable3dPeaks);
      m_Controls->m_Enable3dPeaks->setChecked(enable3dPeaks);
      float lineWidth = static_cast<float>(m_Controls->m_LineWidth->value());
      node->GetFloatProperty("shape.linewidth", lineWidth);
      m_Controls->m_LineWidth->setValue(lineWidth);

      m_Controls->m_Enable3dPeaks->setVisible(true);
      m_Controls->m_3DClippingBox->setVisible(true);
    }
    else if (IsFiberBundle(nodeData))
    {
      this->ShowClippingPlaneId(node);
      m_ColorObservers.push_back(std::make_unique<ColorObserver>(node));

      m_Controls->m_FlipClipBox->setVisible(true);
      m_Controls->m_3DClippingBox->setVisible(true);
      m_Controls->m_BundleControlsFrame->setVisible(true);

      // m_LineWidth belongs to the peak controls; fibers have tube and ribbon width
      float radius = static_cast<float>(m_Controls->m_TubeWidth->value());
      node->GetFloatProperty("shape.tuberadius", radius);
      m_Controls->m_TubeWidth->setValue(radius);

      float ribbonWidth = static_cast<float>(m_Controls->m_RibbonWidth->value());
      node->GetFloatProperty("shape.ribbonwidth", ribbonWidth);
      m_Controls->m_RibbonWidth->setValue(ribbonWidth);

      bool flipClipping = false;
      node->GetBoolProperty("3DClippingPlaneFlip", flipClipping);
      m_Controls->m_FlipClipBox->setChecked(flipClipping);

      float range = static_cast<float>(m_Controls->m_FiberClippingBox->value());
      node->GetFloatProperty("Fiber2DSliceThickness",range);
      mitk::BaseGeometry::Pointer geo = nodeData->GetGeometry();
      mitk::ScalarType max = geo->GetExtentInMM(0);
      max = std::max(max, geo->GetExtentInMM(1));
      max = std::max(max, geo->GetExtentInMM(2));

      m_Controls->m_FiberClippingBox->setMaximum(max);
      m_Controls->m_FiberClippingBox->setValue(range);
    }
    else if(IsGlyphImage(nodeData))
    {
      m_Controls->m_ImageControlsFrame->setVisible(true);
      m_Controls->m_NumberGlyphsFrame->setVisible(true);
      m_Controls->m_GlyphFrame->setVisible(true);
      m_Controls->m_NormalizationFrame->setVisible(true);

      if(m_NodeUsedForOdfVisualization.IsNotNull())
      {
        m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_S", false);
        m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_C", false);
        m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_T", false);
      }
      m_NodeUsedForOdfVisualization = node;
      m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_S", m_GlyIsOn_S);
      m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_C", m_GlyIsOn_C);
      m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_T", m_GlyIsOn_T);

      // tensors have no normalization choice
      const bool isTensorImage = dynamic_cast<mitk::TensorImage*>(nodeData) != nullptr;
      m_Controls->m_NormalizationDropdown->setVisible(!isTensorImage);
      m_Controls->m_NormalizationLabel->setVisible(!isTensorImage);

      int val = m_Controls->m_ShowMaxNumber->value();
      node->GetIntProperty("ShowMaxNumber", val);
      m_Controls->m_ShowMaxNumber->setValue(val);

      // the combo box entries follow the ids of the ODF mapper's enumeration properties; "Global
      // Maximum", which older scenes may hold, draws as "Max" and is shown as such
      if (auto normalization = dynamic_cast<mitk::EnumerationProperty*>(node->GetProperty("Normalization")))
      {
        const auto id = normalization->GetValueAsId();
        m_Controls->m_NormalizationDropdown->setCurrentIndex(id==mitk::ODFN_GLOBAL_MAX ? static_cast<int>(mitk::ODFN_MAX) : static_cast<int>(id));
      }

      float fval = static_cast<float>(m_Controls->m_ScalingFactor->value());
      node->GetFloatProperty("Scaling",fval);
      m_Controls->m_ScalingFactor->setValue(fval);

      if (auto scaleBy = dynamic_cast<mitk::EnumerationProperty*>(node->GetProperty("ScaleBy")))
        m_Controls->m_AdditionalScaling->setCurrentIndex(scaleBy->GetValueAsId());

      bool colourisationModeBit = false;
      node->GetBoolProperty("DiffusionCore.Rendering.OdfVtkMapper.ColourisationModeBit", colourisationModeBit );
      m_Controls->m_OdfColorBox->setCurrentIndex(colourisationModeBit);
    }
    else if(dynamic_cast<mitk::PlanarFigure*>(nodeData))
    {
      lastPlanarFigure = node;
    }
    else if( dynamic_cast<mitk::Image*>(nodeData) )
    {
      m_Controls->m_ImageControlsFrame->setVisible(true);
      m_Controls->m_TSMenu->setVisible(true);
    }

    m_SelectedNodes.push_back(node);
  }

  if (lastPlanarFigure.IsNotNull())
    this->PlanarFigureFocus(lastPlanarFigure);

  if( nodes.empty() ) { return; }

  this->UpdateThickSliceMenu(nodes);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

/** The thick-slice menu is built once; a selection change only updates its state. */
void QmitkControlVisualizationPropertiesView::CreateThickSliceMenu()
{
  m_TSSlider = new QSlider(Qt::Horizontal, m_MyMenu);
  m_TSSlider->setMinimum(0);
  connect( m_TSSlider, &QSlider::valueChanged, this, &QmitkControlVisualizationPropertiesView::OnTSNumChanged );

  QHBoxLayout* tsLayout = new QHBoxLayout;
  tsLayout->setContentsMargins(4,4,4,4);
  tsLayout->addWidget(m_TSSlider);
  tsLayout->addWidget(m_TSLabel = new QLabel(QString::number(1), m_MyMenu));

  QWidget* tsWidget = new QWidget;
  tsWidget->setLayout(tsLayout);

  QWidgetAction* tsSliderAction = new QWidgetAction(m_MyMenu);
  tsSliderAction->setDefaultWidget(tsWidget);
  m_MyMenu->addAction(tsSliderAction);

  m_ThickSliceModeGroup = new QActionGroup(m_MyMenu);
  m_ThickSliceModeGroup->setExclusive(true);
  const std::pair<QString, int> modes[] = {{"MIP (max. intensity proj.)", 1},
                                           {"SUM (sum intensity proj.)", 2},
                                           {"WEIGHTED (gaussian proj.)", 3}};
  for (const auto& mode : modes)
  {
    QAction* action = new QAction(m_MyMenu);
    action->setActionGroup(m_ThickSliceModeGroup);
    action->setText(mode.first);
    action->setCheckable(true);
    action->setData(mode.second);
    m_MyMenu->addAction(action);
  }
  connect( m_ThickSliceModeGroup, &QActionGroup::triggered, this, &QmitkControlVisualizationPropertiesView::OnThickSlicesModeSelected );
}

void QmitkControlVisualizationPropertiesView::UpdateThickSliceMenu(const QList<mitk::DataNode::Pointer>& nodes)
{
  int maxTS = 30;
  for (auto node: nodes)
  {
    mitk::Image* image = node.IsNotNull() ? dynamic_cast<mitk::Image*>(node->GetData()) : nullptr;
    if (image)
    {
      int size = std::max(image->GetDimension(0), std::max(image->GetDimension(1), image->GetDimension(2)));

      if (size>maxTS) { maxTS=size; }
    }
  }
  maxTS /= 2;

  // Thick slices are properties of each 2D renderer's plane geometry node, where SetTs writes
  // them and MITK's render window menu reads them; the axial renderer shows the current state.
  int currentNum = 0;
  int currentMode = 0;
  mitk::BaseRenderer* renderer = GetRendererOf(this->GetRenderWindowPart(), mitk::AnatomicalPlane::Axial);
  mitk::DataNode* planeNode = renderer != nullptr ? renderer->GetCurrentWorldPlaneGeometryNode() : nullptr;
  if (planeNode != nullptr)
  {
    planeNode->GetIntProperty("reslice.thickslices.num", currentNum);
    if (auto mode = dynamic_cast<mitk::ResliceMethodProperty*>(planeNode->GetProperty("reslice.thickslices")))
      currentMode = mode->GetValueAsId();
  }
  if (currentMode == 0)
    currentNum = 0;
  else
    m_CurrentThickSlicesMode = currentMode;

  {
    const QSignalBlocker blocker(m_TSSlider);
    m_TSSlider->setMaximum(maxTS-1);
    m_TSSlider->setValue(currentNum);
  }
  m_CurrentThickSlicesNum = m_TSSlider->value();
  m_TSLabel->setText(QString::number(m_CurrentThickSlicesNum*2+1));

  for (QAction* action : m_ThickSliceModeGroup->actions())
    action->setChecked(action->data().toInt() == m_CurrentThickSlicesMode);
}


void QmitkControlVisualizationPropertiesView::VisibleOdfsON_S()
{
  m_GlyIsOn_S = m_Controls->m_VisibleOdfsON_S->isChecked();
  if (m_NodeUsedForOdfVisualization.IsNull())
  {
    MITK_WARN << "ODF visualization activated but m_NodeUsedForOdfVisualization is nullptr";
    return;
  }
  m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_S", m_GlyIsOn_S);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::FlipPeaks()
{

  for (auto node : m_SelectedNodes)
  {
    auto peaks = dynamic_cast<mitk::PeakImage*>(node->GetData());
    if (peaks==nullptr)
      continue;

    // The flipped peaks are written back into the node's image: one copy of the buffer instead
    // of two into a new image, and the image keeps its properties.
    mitk::PeakImage::ItkPeakImageType::Pointer flipped;
    {
      // The caster reads the image's buffer in place; its read access must end with this scope,
      // before the buffer is written.
      using CasterType = mitk::ImageToItk< mitk::PeakImage::ItkPeakImageType >;
      CasterType::Pointer caster = CasterType::New();
      caster->SetInput(static_cast<const mitk::Image*>(peaks));
      caster->Update();

      itk::FlipPeaksFilter::Pointer flipper = itk::FlipPeaksFilter::New();
      flipper->SetInput(caster->GetOutput());
      flipper->SetFlipX(m_Controls->m_FlipPeaksX->isChecked());
      flipper->SetFlipY(m_Controls->m_FlipPeaksY->isChecked());
      flipper->SetFlipZ(m_Controls->m_FlipPeaksZ->isChecked());
      flipper->Update();
      flipped = flipper->GetOutput();
      flipped->DisconnectPipeline();
    }

    std::size_t numValues = 1;
    for (unsigned int d=0; d<peaks->GetDimension(); ++d)
      numValues *= peaks->GetDimension(d);
    if (flipped->GetPixelContainer()->Size()!=numValues)
    {
      MITK_ERROR << "Flipped peaks of " << node->GetName() << " do not match the image size; image left unchanged.";
      continue;
    }
    {
      mitk::ImageWriteAccessor access(peaks);
      std::memcpy(access.GetData(), flipped->GetBufferPointer(), numValues*sizeof(mitk::PeakImage::PixelType));
    }
    // GetPolyData() rebuilds the peak lines the mappers draw once the image is newer, in its
    // current colour
    peaks->Modified();
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::Toggle3DPeaks()
{
  ForSelected(IsPeakImage, [this](mitk::DataNode* node)
  {
    node->SetBoolProperty( "Enable3DPeaks", m_Controls->m_Enable3dPeaks->isChecked() );
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::Toggle3DClipping()
{
  int planeId = 0;
  if (m_Controls->m_Clip1->isChecked())
    planeId = 1;
  else if (m_Controls->m_Clip2->isChecked())
    planeId = 2;
  else if (m_Controls->m_Clip3->isChecked())
    planeId = 3;

  ForSelected(IsFiberBundleOrPeakImage, [this, planeId](mitk::DataNode* node)
  {
    node->SetBoolProperty( "3DClippingPlaneFlip", m_Controls->m_FlipClipBox->isChecked() );
    node->SetIntProperty( "3DClippingPlaneId", planeId );
    this->Set3DClippingPlane(node, planeId);
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::OnSliceChanged()
{
  mitk::DataStorage::SetOfObjects::ConstPointer nodes = this->GetDataStorage()->GetAll();
  for (const auto& node : *nodes)
  {
    int planeId = 0;
    node->GetIntProperty("3DClippingPlaneId", planeId);
    if (planeId>=1 && planeId<=3)
      this->Set3DClippingPlane(node, planeId);
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::Set3DClippingPlane(mitk::DataNode* node, int planeId)
{
  mitk::IRenderWindowPart* renderWindow = this->GetRenderWindowPart();
  if (renderWindow && node && IsFiberBundleOrPeakImage(node->GetData()))
  {
    mitk::Vector3D planeNormal; planeNormal.Fill(0.0);
    mitk::Point3D planeOrigin; planeOrigin.Fill(0.0);
    if (planeId != 0)
    {
      // the window that shows this plane, whatever its name in the editor
      QmitkRenderWindow* planeWindow = renderWindow->GetQmitkRenderWindow(ClippingPlaneOf(planeId));
      mitk::SliceNavigationController* slicer = planeWindow != nullptr ? planeWindow->GetSliceNavigationController() : nullptr;
      mitk::PlaneGeometry::ConstPointer planeGeo = slicer != nullptr ? slicer->GetCurrentPlaneGeometry() : nullptr;
      if (planeGeo.IsNull())
        return;
      planeOrigin = renderWindow->GetSelectedPosition();
      planeNormal = planeGeo->GetNormal();
    }
    // The mappers read the plane at every render (the fiber mapper without rebuilding its tubes);
    // the callers request the render.
    node->SetProperty( "3DClipping", mitk::ClippingProperty::New( planeOrigin, planeNormal ) );
  }
}

void QmitkControlVisualizationPropertiesView::VisibleOdfsON_T()
{
  m_GlyIsOn_T = m_Controls->m_VisibleOdfsON_T->isChecked();
  if (m_NodeUsedForOdfVisualization.IsNull())
  {
    MITK_WARN << "ODF visualization activated but m_NodeUsedForOdfVisualization is nullptr";
    return;
  }
  m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_T", m_GlyIsOn_T);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::VisibleOdfsON_C()
{
  m_GlyIsOn_C = m_Controls->m_VisibleOdfsON_C->isChecked();
  if (m_NodeUsedForOdfVisualization.IsNull())
  {
    MITK_WARN << "ODF visualization activated but m_NodeUsedForOdfVisualization is nullptr";
    return;
  }
  m_NodeUsedForOdfVisualization->SetBoolProperty("VisibleOdfs_C", m_GlyIsOn_C);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::ShowMaxNumberChanged()
{
  // the spin box minimum is 1
  const int maxNr = m_Controls->m_ShowMaxNumber->value();
  ForSelected(IsGlyphImage, [maxNr](mitk::DataNode* node) { node->SetIntProperty("ShowMaxNumber", maxNr); });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::NormalizationDropdownChanged(int normDropdown)
{
  using PropType = mitk::OdfNormalizationMethodProperty;
  PropType::Pointer normMeth = PropType::New();

  switch(normDropdown)
  {
  case 0:
    normMeth->SetNormalizationToMinMax();
    break;
  case 1:
    normMeth->SetNormalizationToMax();
    break;
  case 2:
    normMeth->SetNormalizationToNone();
    break;
  default:
    normMeth->SetNormalizationToMinMax();
  }

  ForSelected(IsGlyphImage, [&normMeth](mitk::DataNode* node) { node->SetProperty("Normalization", normMeth.GetPointer()); });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::ScalingFactorChanged(double scalingFactor)
{
  ForSelected(IsGlyphImage, [scalingFactor](mitk::DataNode* node) { node->SetFloatProperty("Scaling", scalingFactor); });

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}


void QmitkControlVisualizationPropertiesView::AdditionalScaling(int additionalScaling)
{
  using PropType = mitk::OdfScaleByProperty;
  PropType::Pointer scaleBy = PropType::New();

  switch(additionalScaling)
  {
  case 0:
    scaleBy->SetScaleByNothing();
    break;
  case 1:
    scaleBy->SetScaleByGFA();
    break;
  case 2:
    scaleBy->SetScaleByPrincipalCurvature();
    break;
  default:
    scaleBy->SetScaleByNothing();
  }

  ForSelected(IsGlyphImage, [&scaleBy](mitk::DataNode* node) { node->SetProperty("ScaleBy", scaleBy.GetPointer()); });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::FiberSlicingThickness2D()
{
  const float fibThickness = m_Controls->m_FiberClippingBox->value();
  ForSelected(IsFiberBundle, [fibThickness](mitk::DataNode* node)
  {
    float currentThickness = 0;
    node->GetFloatProperty("Fiber2DSliceThickness", currentThickness);
    if ( fabs(fibThickness-currentThickness) < 0.001 )
      return;

    // the 2D mapper reads the thickness on every render; RequestUpdate2D would give each 2D
    // window a new input copy and a new upload of the whole bundle
    node->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(fibThickness));
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::ResetColoring()
{
  ForSelected(IsFiberBundleOrPeakImage, [](mitk::DataNode* node)
  {
    if (auto fib = dynamic_cast<mitk::FiberBundle*>(node->GetData()))
      fib->ColorFibersByOrientation();
    else if (auto peaks = dynamic_cast<mitk::PeakImage*>(node->GetData()))
      peaks->ColorByOrientation();
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkControlVisualizationPropertiesView::PlanarFigureFocus(mitk::DataNode* node)
{
  mitk::PlanarFigure* _PlanarFigure = dynamic_cast<mitk::PlanarFigure*> (node->GetData());
  auto renderWindowPart = this->GetRenderWindowPart(mitk::WorkbenchUtil::IRenderWindowPartStrategy::OPEN);

  if (_PlanarFigure && _PlanarFigure->GetPlaneGeometry() && renderWindowPart)
  {
    // An editor may not offer all three 2D windows (MxN layouts).
    QmitkRenderWindow* axialRenderWindow = renderWindowPart->GetQmitkRenderWindow(mitk::AnatomicalPlane::Axial);
    QmitkRenderWindow* sagittalRenderWindow = renderWindowPart->GetQmitkRenderWindow(mitk::AnatomicalPlane::Sagittal);
    QmitkRenderWindow* coronalRenderWindow = renderWindowPart->GetQmitkRenderWindow(mitk::AnatomicalPlane::Coronal);

    if (axialRenderWindow && sagittalRenderWindow && coronalRenderWindow)
    {
      QmitkRenderWindow* selectedRenderWindow = nullptr;
      const mitk::PlaneGeometry* _PlaneGeometry = _PlanarFigure->GetPlaneGeometry();

      mitk::VnlVector normal = _PlaneGeometry->GetNormalVnl();
      mitk::VnlVector normal1 = axialRenderWindow->GetRenderer()->GetCurrentWorldPlaneGeometry()->GetNormalVnl();
      mitk::VnlVector normal2 = sagittalRenderWindow->GetRenderer()->GetCurrentWorldPlaneGeometry()->GetNormalVnl();
      mitk::VnlVector normal3 = coronalRenderWindow->GetRenderer()->GetCurrentWorldPlaneGeometry()->GetNormalVnl();

      normal[0]  = fabs(normal[0]);  normal[1]  = fabs(normal[1]);  normal[2]  = fabs(normal[2]);
      normal1[0] = fabs(normal1[0]); normal1[1] = fabs(normal1[1]); normal1[2] = fabs(normal1[2]);
      normal2[0] = fabs(normal2[0]); normal2[1] = fabs(normal2[1]); normal2[2] = fabs(normal2[2]);
      normal3[0] = fabs(normal3[0]); normal3[1] = fabs(normal3[1]); normal3[2] = fabs(normal3[2]);

      double ang1 = angle(normal, normal1);
      double ang2 = angle(normal, normal2);
      double ang3 = angle(normal, normal3);

      if(ang1 < ang2 && ang1 < ang3)
        selectedRenderWindow = axialRenderWindow;
      else if(ang2 < ang3)
        selectedRenderWindow = sagittalRenderWindow;
      else
        selectedRenderWindow = coronalRenderWindow;

      // make node visible
      const mitk::Point3D& centerP = _PlaneGeometry->GetOrigin();
      selectedRenderWindow->GetSliceNavigationController()->ReorientSlices(centerP, _PlaneGeometry->GetNormal());
    }
  }

  // set interactor for new node (if not already set)
  mitk::PlanarFigureInteractor::Pointer figureInteractor
      = dynamic_cast<mitk::PlanarFigureInteractor*>(node->GetDataInteractor().GetPointer());

  if(figureInteractor.IsNull())
  {
    figureInteractor = mitk::PlanarFigureInteractor::New();
    us::Module* planarFigureModule = us::ModuleRegistry::GetModule( "MitkPlanarFigure" );
    figureInteractor->LoadStateMachine("PlanarFigureInteraction.xml", planarFigureModule );
    figureInteractor->SetEventConfig( "PlanarFigureConfig.xml", planarFigureModule );
    figureInteractor->SetDataNode( node );
  }

  node->SetProperty("planarfigure.iseditable",mitk::BoolProperty::New(true));
}


void QmitkControlVisualizationPropertiesView::TubeRadiusChanged()
{
  const float newRadius = m_Controls->m_TubeWidth->value();
  ForSelected(IsFiberBundle, [newRadius](mitk::DataNode* node) { node->SetFloatProperty("shape.tuberadius", newRadius); });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::RibbonWidthChanged()
{
  const float newWidth = m_Controls->m_RibbonWidth->value();
  ForSelected(IsFiberBundle, [newWidth](mitk::DataNode* node) { node->SetFloatProperty("shape.ribbonwidth", newWidth); });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::LineWidthChanged()
{
  // the peak mappers read the line width as float
  const auto newWidth = static_cast<float>(m_Controls->m_LineWidth->value());
  ForSelected(IsPeakImage, [newWidth](mitk::DataNode* node)
  {
    float currentWidth = 0;
    node->GetFloatProperty("shape.linewidth", currentWidth);
    if (currentWidth!=newWidth)
      node->SetFloatProperty("shape.linewidth", newWidth);
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkControlVisualizationPropertiesView::OnColourisationModeChanged()
{
  // an ODF mapper property; only glyph images get it
  const bool colourisationModeBit = m_Controls->m_OdfColorBox->currentIndex();
  ForSelected(IsGlyphImage, [colourisationModeBit](mitk::DataNode* node)
  {
    node->SetProperty( "DiffusionCore.Rendering.OdfVtkMapper.ColourisationModeBit", mitk::BoolProperty::New( colourisationModeBit ) );
  });
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}
