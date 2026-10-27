/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkStreamlineInteractorBrush.h"
//#include "mitkStreamlineMapper2D.h"

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
#include <vtkPolyLine.h>
#include <vtkFloatArray.h>

mitk::StreamlineInteractorBrush::StreamlineInteractorBrush()
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

mitk::StreamlineInteractorBrush::~StreamlineInteractorBrush()
{
}

void mitk::StreamlineInteractorBrush::ConnectActionsAndFunctions()
{
    CONNECT_CONDITION("isoverstreamline", CheckSelection);
    CONNECT_FUNCTION("selectstreamline", SelectStreamline);

  CONNECT_FUNCTION("addnegstreamline", AddStreamlineNegBundle);
  CONNECT_FUNCTION("addposstreamline", AddStreamlinePosBundle);
  CONNECT_FUNCTION("addtolabelstreamline", AddStreamlinetolabelsBundle);
}

void mitk::StreamlineInteractorBrush::SetNegativeNode(DataNode *node)
{

    m_NegStreamlineNode = node;
    m_NegStreamline= dynamic_cast<mitk::FiberBundle *>(node->GetData());
}

void mitk::StreamlineInteractorBrush::SetToLabelNode(DataNode *node)
{
    m_manStreamlineNode = node;
    DataInteractor::SetDataNode(m_manStreamlineNode);
    m_manStreamline = dynamic_cast<mitk::FiberBundle *>(node->GetData());

}

void mitk::StreamlineInteractorBrush::SetPositiveNode(DataNode *node)
{
    m_PosStreamlineNode = node;
    m_PosStreamline= dynamic_cast<mitk::FiberBundle *>(node->GetData());
}

void mitk::StreamlineInteractorBrush::LabelfromPrediction(bool predlabeling)
{
    m_predlabeling = predlabeling;

}

bool mitk::StreamlineInteractorBrush::CheckSelection(const InteractionEvent *interactionEvent)
{
    /* Save Colorvector here*/

    const auto *positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent != nullptr)
    {
        if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
        {
            BaseRenderer *renderer = positionEvent->GetSender();

            auto &picker = m_Picker[renderer];

            picker = vtkSmartPointer<vtkCellPicker>::New();
            picker->SetTolerance(0.01);
            auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard2D);


              auto vtk_mapper = dynamic_cast<VtkMapper *>(mapper);
              if (vtk_mapper)
              { // doing this each time is bizarre
                picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
                picker->PickFromListOn();
              }

            auto displayPosition = positionEvent->GetPointerPositionOnScreen();
            picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());



            if (picker->GetCellId()==-1 && m_predlabeling==false)
            {

                RenderingManager::GetInstance()->RequestUpdateAll();
                return false;
            }
            else
            {

                return true;
            }
        }
        else {
            BaseRenderer *renderer = positionEvent->GetSender();

            auto &picker = m_Picker[renderer];

            picker = vtkSmartPointer<vtkCellPicker>::New();
            picker->SetTolerance(0.01);
            auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard3D);


              auto vtk_mapper = dynamic_cast<VtkMapper *>(mapper);
              if (vtk_mapper)
              { // doing this each time is bizarre
                picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
                picker->PickFromListOn();
              }

            auto displayPosition = positionEvent->GetPointerPositionOnScreen();
            picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());



            if (picker->GetCellId()==-1 && m_predlabeling==false)
            {
                RenderingManager::GetInstance()->RequestUpdateAll();
                return false;
            }
            else
            {

                return true;
            }

        }


    }
    else {
        RenderingManager::GetInstance()->RequestUpdateAll();
        return false;
    }

}

void mitk::StreamlineInteractorBrush::SelectStreamline(StateMachineAction *, InteractionEvent *interactionEvent)
{
    auto *positionEvent = dynamic_cast<InteractionPositionEvent *>(interactionEvent);

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
        BaseRenderer *renderer = positionEvent->GetSender();

        auto &picker = m_Picker[renderer];

        picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.01);
        auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard2D);


          auto vtk_mapper = dynamic_cast<VtkMapper *>(mapper);
          if (vtk_mapper)
          { // doing this each time is bizarre
            picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
            picker->PickFromListOn();
          }


        auto displayPosition = positionEvent->GetPointerPositionOnScreen();
        picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());
        if (m_predlabeling==false)
            {
        vtkIdType pickedCellID = picker->GetCellId();
        m_manStreamline->SetSingleFiberColor(0.0, 255.0, 0.0, pickedCellID);

        }
    }
    else {
        BaseRenderer *renderer = positionEvent->GetSender();

        auto &picker = m_Picker[renderer];

        picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.01);
        auto mapper = GetDataNode()->GetMapper(BaseRenderer::Standard3D);


          auto vtk_mapper = dynamic_cast<VtkMapper *>(mapper);
          if (vtk_mapper)
          { // doing this each time is bizarre
            picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
            picker->PickFromListOn();
          }
    //        }

        auto displayPosition = positionEvent->GetPointerPositionOnScreen();
        picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());
        if (m_predlabeling==false)
            {
        vtkIdType pickedCellID = picker->GetCellId();
        m_manStreamline->SetSingleFiberColor(0.0, 255.0, 0.0, pickedCellID);
//        m_manStreamline->SetFloatProperty("shape.tuberadius", 1.0);
        }
    }


    RenderingManager::GetInstance()->RequestUpdateAll();
}

vtkIdType mitk::StreamlineInteractorBrush::PickFiber(const InteractionEvent *interactionEvent)
{
    const auto *positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr || GetDataNode() == nullptr)
      return -1;

    BaseRenderer *renderer = positionEvent->GetSender();
    const auto mapperId = renderer->GetMapperID() == BaseRenderer::Standard2D ? BaseRenderer::Standard2D : BaseRenderer::Standard3D;

    auto &picker = m_Picker[renderer];
    picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.01);
    auto vtk_mapper = dynamic_cast<VtkMapper *>(GetDataNode()->GetMapper(mapperId));
    if (vtk_mapper)
    {
      picker->AddPickList(vtk_mapper->GetVtkProp(renderer));
      picker->PickFromListOn();
    }

    auto displayPosition = positionEvent->GetPointerPositionOnScreen();
    picker->Pick(displayPosition[0], displayPosition[1], 0, renderer->GetVtkRenderer());
    return picker->GetCellId();
}

void mitk::StreamlineInteractorBrush::MoveFiber(FiberBundle *source, FiberBundle *target, unsigned int fiberId)
{
    if (source == nullptr || fiberId >= source->GetNumFibers())
      mitkThrow() << "Fiber " << fiberId << " is not part of the bundle to move it from.";

    if (target != nullptr)
    {
      auto movedWeights = vtkSmartPointer<vtkFloatArray>::New();
      auto moved = FiberBundle::New(source->GeneratePolyDataByIds({fiberId}, movedWeights));
      moved->SetFiberWeights(movedWeights);
      auto joined = target->AddBundle(moved);
      target->SetFiberPolyData(joined->GetFiberPolyData());
      target->SetFiberWeights(joined->GetFiberWeights());
    }

    std::vector<unsigned int> remaining;
    remaining.reserve(source->GetNumFibers() - 1);
    for (unsigned int i = 0; i < source->GetNumFibers(); ++i)
    {
      if (i != fiberId)
        remaining.push_back(i);
    }
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    auto polyData = source->GeneratePolyDataByIds(remaining, weights);
    source->SetFiberPolyData(polyData);
    source->SetFiberWeights(weights);
}

namespace
{
  bool IsFiberOf(const mitk::FiberBundle *bundle, vtkIdType cellId)
  {
    return bundle != nullptr && cellId >= 0 && cellId < static_cast<vtkIdType>(bundle->GetNumFibers());
  }
}

void mitk::StreamlineInteractorBrush::AddStreamlinePosBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    DataInteractor::SetDataNode(m_manStreamlineNode);

    const vtkIdType pickedCellID = this->PickFiber(interactionEvent);
    if (IsFiberOf(m_manStreamline, pickedCellID) && m_PosStreamline.IsNotNull())
    {
      MoveFiber(m_manStreamline, m_PosStreamline, static_cast<unsigned int>(pickedCellID));
      m_PosStreamline->SetFiberColors(0, 255, 0);
      if (!m_predlabeling)
        m_manStreamline->SetFiberColors(255, 255, 255);
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}

void mitk::StreamlineInteractorBrush::AddStreamlineNegBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    DataInteractor::SetDataNode(m_manStreamlineNode);

    // On the prediction, the fiber is only excluded from it, not labeled negative.
    const vtkIdType pickedCellID = this->PickFiber(interactionEvent);
    if (IsFiberOf(m_manStreamline, pickedCellID) && (m_predlabeling || m_NegStreamline.IsNotNull()))
    {
      MoveFiber(m_manStreamline, m_predlabeling ? nullptr : m_NegStreamline.GetPointer(), static_cast<unsigned int>(pickedCellID));
      if (!m_predlabeling)
      {
        m_NegStreamline->SetFiberColors(255, 0, 0);
        m_manStreamline->SetFiberColors(255, 255, 255);
      }
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}

void mitk::StreamlineInteractorBrush::AddStreamlinetolabelsBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    // Unlabel: a picked fiber of the visible label bundles goes back to the bundle being labeled.
    struct LabelBundle
    {
      mitk::DataNode *node;
      mitk::FiberBundle *bundle;
      float r, g, b;
    };
    const LabelBundle labelBundles[] = {{m_NegStreamlineNode, m_NegStreamline, 255, 0, 0},
                                        {m_PosStreamlineNode, m_PosStreamline, 0, 255, 0}};

    for (const auto &[node, bundle, r, g, b] : labelBundles)
    {
      if (node == nullptr || !node->IsVisible(nullptr))
        continue;

      DataInteractor::SetDataNode(node);
      const vtkIdType pickedCellID = this->PickFiber(interactionEvent);
      if (!IsFiberOf(bundle, pickedCellID) || m_manStreamline.IsNull())
        continue;

      MoveFiber(bundle, m_manStreamline, static_cast<unsigned int>(pickedCellID));
      bundle->SetFiberColors(r, g, b);
      if (!m_predlabeling)
        m_manStreamline->SetFiberColors(255, 255, 255);
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
}
