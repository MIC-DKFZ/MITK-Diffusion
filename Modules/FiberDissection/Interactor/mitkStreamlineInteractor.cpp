# include

/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkStreamlineInteractor.h"
#include "mitkStreamlineMapper2D.h"

// MITK includes
#include <mitkInteractionConst.h>
#include <mitkInteractionPositionEvent.h>
#include <mitkInternalEvent.h>
#include <mitkLookupTableProperty.h>
#include <mitkOperationEvent.h>
#include <mitkRotationOperation.h>
#include <mitkScaleOperation.h>
#include <mitkSurface.h>
#include <mitkUndoController.h>
#include <mitkVtkMapper.h>

// VTK includes
#include <vtkCamera.h>
#include <vtkInteractorObserver.h>
#include <vtkInteractorStyle.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVector.h>
#include <vtkVectorOperators.h>

mitk::StreamlineInteractor::StreamlineInteractor()
{
    m_ColorForHighlight[0] = 1.0;
    m_ColorForHighlight[1] = 0.5;
    m_ColorForHighlight[2] = 0.0;
    m_ColorForHighlight[3] = 1.0;

    // TODO if we want to get this configurable, the this is the recipe:
    // - make the 2D mapper add corresponding properties to control "enabled" and "color"
    // - make the interactor evaluate those properties
    // - in an ideal world, modify the state machine on the fly and skip mouse move handling
}

mitk::StreamlineInteractor::~StreamlineInteractor()
{
}

void mitk::StreamlineInteractor::ConnectActionsAndFunctions()
{
  CONNECT_CONDITION("PickedHandle", HasPickedHandle);

  CONNECT_FUNCTION("FeedUndoStack", FeedUndoStack);
}

void mitk::StreamlineInteractor::SetStreamlineNode(DataNode *node)
{
  DataInteractor::SetDataNode(node);

  m_Streamline = dynamic_cast<Streamline *>(node->GetData());

  // setup picking from just this object
  m_Picker.clear();
}

void mitk::StreamlineInteractor::SetManipulatedObjectNode(DataNode *node)
{
  if (node && node->GetData())
  {
    m_ManipulatedObjectGeometry = node->GetData()->GetGeometry();
  }
}

bool mitk::StreamlineInteractor::HasPickedHandle(const InteractionEvent *interactionEvent)
{
  auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
  if (positionEvent == nullptr ||
      m_Streamline.IsNull() ||
      m_ManipulatedObjectGeometry.IsNull() ||
      interactionEvent->GetSender()->GetRenderWindow()->GetNeverRendered())
  {
    return false;
  }

  if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
  {
    m_PickedHandle = PickFrom2D(positionEvent);
  }
  else
  {
    m_PickedHandle = PickFrom3D(positionEvent);
  }

  UpdateHandleHighlight();

  return m_PickedHandle != Streamline::NoHandle;
}

void mitk::StreamlineInteractor::FeedUndoStack(StateMachineAction *, InteractionEvent *)
{
  if (m_UndoEnabled)
  {
    OperationEvent *operationEvent = new OperationEvent(m_ManipulatedObjectGeometry,
                                                        // OperationEvent will destroy operations!
                                                        // --> release() and not get()
                                                        m_FinalDoOperation.release(),
                                                        m_FinalUndoOperation.release(),
                                                        "Direct geometry manipulation");
    mitk::OperationEvent::IncCurrObjectEventId(); // save each modification individually
    m_UndoController->SetOperationEvent(operationEvent);
  }
}

mitk::Streamline::HandleType mitk::StreamlineInteractor::PickFrom2D(const InteractionPositionEvent *positionEvent)
{
  BaseRenderer *renderer = positionEvent->GetSender();

  auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard2D);
  auto Streamline_mapper = dynamic_cast<StreamlineMapper2D *>(mapper);
  auto &picker = m_Picker[renderer];

  if (picker == nullptr)
  {
    picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.005);

    if (Streamline_mapper)
    { // doing this each time is bizarre
      picker->AddPickList(Streamline_mapper->GetVtkProp(renderer));
      picker->PickFromListOn();
    }
  }

  auto displayPosition = positionEvent->GetPointerPositionOnScreen();
  picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

  vtkIdType pickedPointID = picker->GetPointId();
  if (pickedPointID == -1)
  {
    return Streamline::NoHandle;
  }

  vtkPolyData *polydata = Streamline_mapper->GetVtkPolyData(renderer);

  if (polydata && polydata->GetPointData() && polydata->GetPointData()->GetScalars())
  {
    double dataValue = polydata->GetPointData()->GetScalars()->GetTuple1(pickedPointID);
    return m_Streamline->GetHandleFromPointDataValue(dataValue);
  }

  return Streamline::NoHandle;
}

mitk::Streamline::HandleType mitk::StreamlineInteractor::PickFrom3D(const InteractionPositionEvent *positionEvent)
{
  BaseRenderer *renderer = positionEvent->GetSender();
  auto &picker = m_Picker[renderer];
  if (picker == nullptr)
  {
    picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.005);
    auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard3D);
    auto vtk_mapper = dynamic_cast<VtkMapper *>(mapper);
    if (vtk_mapper)
    { // doing this each time is bizarre
      picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
      picker->PickFromListOn();
    }
  }

  auto displayPosition = positionEvent->GetPointerPositionOnScreen();
  picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

  vtkIdType pickedPointID = picker->GetPointId();
  if (pickedPointID == -1)
  {
    return Streamline::NoHandle;
  }

  // _something_ picked
  return m_Streamline->GetHandleFromPointID(pickedPointID);
}

void mitk::StreamlineInteractor::UpdateHandleHighlight()
{
  if (m_HighlightedHandle != m_PickedHandle) {

    auto node = GetDataNode();
    if (node == nullptr) return;

    auto base_prop = node->GetProperty("LookupTable");
    if (base_prop == nullptr) return;

    auto lut_prop = dynamic_cast<LookupTableProperty*>(base_prop);
    if (lut_prop == nullptr) return;

    auto lut = lut_prop->GetLookupTable();
    if (lut == nullptr) return;

    // Table size is expected to constructed as one entry per Streamline-part enum value
    assert(lut->GetVtkLookupTable()->GetNumberOfTableValues() > std::max(m_PickedHandle, m_HighlightedHandle));

    // Reset previously overwritten color
    if (m_HighlightedHandle != Streamline::NoHandle)
    {
        lut->SetTableValue(m_HighlightedHandle, m_ColorReplacedByHighlight);
    }

    // Overwrite currently highlighted color
    if (m_PickedHandle != Streamline::NoHandle)
    {
      lut->GetTableValue(m_PickedHandle, m_ColorReplacedByHighlight);
      lut->SetTableValue(m_PickedHandle, m_ColorForHighlight);
    }

    // Mark node modified to allow repaint
    node->Modified();
    RenderingManager::GetInstance()->RequestUpdateAll(RenderingManager::REQUEST_UPDATE_ALL);

    m_HighlightedHandle = m_PickedHandle;
  }
}
