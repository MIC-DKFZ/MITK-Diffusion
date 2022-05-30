/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkStreamlineInteractor.h"
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
//  CONNECT_CONDITION("isoverstreamline", HasPickedHandle);
  CONNECT_FUNCTION("addnegstreamline", AddStreamlineNegBundle);
  CONNECT_FUNCTION("addposstreamline", AddStreamlinePosBundle);
//  CONNECT_FUNCTION("FeedUndoStack", FeedUndoStack);
}

void mitk::StreamlineInteractor::SetNegativeNode(DataNode *node)
{

    DataInteractor::SetDataNode(node);
    m_NegStreamline= dynamic_cast<mitk::FiberBundle *>(node->GetData());
    MITK_INFO << "Negative Node added";
}

void mitk::StreamlineInteractor::SetToLabelNode(DataNode *node)
{
    DataInteractor::SetDataNode(node);
    m_manStreamline = dynamic_cast<mitk::FiberBundle *>(node->GetData());
    MITK_INFO << "Label node added";

}

void mitk::StreamlineInteractor::SetPositiveNode(DataNode *node)
{

    DataInteractor::SetDataNode(node);
    m_PosStreamline= dynamic_cast<mitk::FiberBundle *>(node->GetData());
    MITK_INFO << "Positive Node added";
  }

void mitk::StreamlineInteractor::AddStreamlinePosBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    MITK_INFO << "PositiveBundle clicked";

    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {
        MITK_INFO << "no position";
    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);
        MITK_INFO << "2D";
    }
    else
    {
        BaseRenderer *renderer = positionEvent->GetSender();

        auto &picker = m_Picker[renderer];
//        if (picker == nullptr)
//        {



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

        vtkIdType pickedCellID = picker->GetCellId();



        if (picker->GetCellId()==-1)
        {
            MITK_INFO << "Nothing picked";
        }
        else {


            vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

            unsigned int counter = 0;
            for ( int i=0; i<m_PosStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
            {
              vtkCell* cell = m_PosStreamline->GetFiberPolyData()->GetCell(i);
              auto numPoints = cell->GetNumberOfPoints();
              vtkPoints* points = cell->GetPoints();

              vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
              for (unsigned int j=0; j<numPoints; j++)
              {
                double p[3];
                points->GetPoint(j, p);

                vtkIdType id = vNewPoints->InsertNextPoint(p);
                container->GetPointIds()->InsertNextId(id);
              }
          //    weights->InsertValue(counter, fib->GetFiberWeight(i));
              vNewLines->InsertNextCell(container);
              counter++;

            }



            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(pickedCellID);

            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
            for (unsigned int j=0; j<numPoints; j++)
            {
              double p[3];
              points->GetPoint(j, p);

              vtkIdType id = vNewPoints->InsertNextPoint(p);
              container->GetPointIds()->InsertNextId(id);
            }
            vNewLines->InsertNextCell(container);

            vNewPolyData->SetPoints(vNewPoints);
            vNewPolyData->SetLines(vNewLines);

    //        m_PosStreamline = mitk::FiberBundle::New(vNewPolyData);
            m_PosStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
            m_PosStreamline->GetFiberPolyData()->SetLines(vNewLines);
            m_PosStreamline->SetFiberColors(0, 255, 0);

            m_manStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
            m_manStreamline->GetFiberPolyData()->RemoveDeletedCells();
        }
    }
  }

void mitk::StreamlineInteractor::AddStreamlineNegBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    MITK_INFO << "NegativeBundle clicked";
//     auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
//     if (positionEvent == nullptr)
//     {
//       return;
//     }
//    return true;
    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {
        MITK_INFO << "no position";
    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);
        MITK_INFO << "2D";
    }
    else
    {
        BaseRenderer *renderer = positionEvent->GetSender();

        auto &picker = m_Picker[renderer];
//        if (picker == nullptr)
//        {


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

        vtkIdType pickedCellID = picker->GetCellId();


        if (picker->GetCellId()==-1)
        {
            MITK_INFO << "Nothing picked";
        }
        else
        {
            vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

            unsigned int counter = 0;
            for ( int i=0; i<m_NegStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
            {
              vtkCell* cell = m_NegStreamline->GetFiberPolyData()->GetCell(i);
              auto numPoints = cell->GetNumberOfPoints();
              vtkPoints* points = cell->GetPoints();

              vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
              for (unsigned int j=0; j<numPoints; j++)
              {
                double p[3];
                points->GetPoint(j, p);

                vtkIdType id = vNewPoints->InsertNextPoint(p);
                container->GetPointIds()->InsertNextId(id);
              }
          //    weights->InsertValue(counter, fib->GetFiberWeight(i));
              vNewLines->InsertNextCell(container);
              counter++;

            }



            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(pickedCellID);
            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
            for (unsigned int j=0; j<numPoints; j++)
            {
              double p[3];
              points->GetPoint(j, p);

              vtkIdType id = vNewPoints->InsertNextPoint(p);
              container->GetPointIds()->InsertNextId(id);
            }
            vNewLines->InsertNextCell(container);

            vNewPolyData->SetPoints(vNewPoints);
            vNewPolyData->SetLines(vNewLines);

    //        m_NegStreamline = mitk::FiberBundle::New(vNewPolyData);
            m_NegStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
            m_NegStreamline->GetFiberPolyData()->SetLines(vNewLines);
            m_NegStreamline->SetFiberColors(255, 0, 0);

            m_manStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
            m_manStreamline->GetFiberPolyData()->RemoveDeletedCells();

        }









    }



  }
