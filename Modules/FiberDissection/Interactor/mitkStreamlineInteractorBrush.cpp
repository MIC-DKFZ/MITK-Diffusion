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
#include <vtkVectorOperators.h>

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
//  CONNECT_CONDITION("isoverstreamline", HasPickedHandle);
    CONNECT_CONDITION("isoverstreamline", CheckSelection);
    CONNECT_FUNCTION("selectstreamline", SelectStreamline);

  CONNECT_FUNCTION("addnegstreamline", AddStreamlineNegBundle);
  CONNECT_FUNCTION("addposstreamline", AddStreamlinePosBundle);
  CONNECT_FUNCTION("addnegtolabelstreamline", AddNegStreamlinetolabelsBundle);
  CONNECT_FUNCTION("addpostolabelstreamline", AddPosStreamlinetolabelsBundle);
//  CONNECT_FUNCTION("FeedUndoStack", FeedUndoStack);
}

void mitk::StreamlineInteractorBrush::SetNegativeNode(DataNode *node)
{

//    DataInteractor::SetDataNode(node);
    m_NegStreamlineNode = node;
    m_NegStreamline= dynamic_cast<mitk::FiberBundle *>(node->GetData());
}

void mitk::StreamlineInteractorBrush::SetToLabelNode(DataNode *node)
{
    m_manStreamlineNode = node;
    DataInteractor::SetDataNode(m_manStreamlineNode);
    m_manStreamline = dynamic_cast<mitk::FiberBundle *>(node->GetData());
//    m_init = false;

}

void mitk::StreamlineInteractorBrush::SetPositiveNode(DataNode *node)
{

    //    DataInteractor::SetDataNode(node);
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
//                m_manStreamline->SetFiberColors(255, 255, 255);
//                vtkSmartPointer<vtkUnsignedCharArray> FiberColors = m_manStreamline->m_FiberColors;

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
//                m_manStreamline->SetFiberColors(255, 255, 255);
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
//        m_manStreamline->SetFiberColors(255, 255, 255);
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
    //        }

        auto displayPosition = positionEvent->GetPointerPositionOnScreen();
        picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());
        if (m_predlabeling==false)
            {
        vtkIdType pickedCellID = picker->GetCellId();
        m_manStreamline->SetSingleFiberColor(0.0, 255.0, 0.0, pickedCellID);
//        m_manStreamline->SetFloatProperty("shape.tuberadius", 1.0)
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

void mitk::StreamlineInteractorBrush::AddStreamlinePosBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    DataInteractor::SetDataNode(m_manStreamlineNode);

    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {

    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);
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
//        }

      auto displayPosition = positionEvent->GetPointerPositionOnScreen();
      picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

      vtkIdType pickedCellID = picker->GetCellId();



      if (picker->GetCellId()==-1)
      {
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
//            m_PosStreamline->SetFiberWeights(m_PosStreamline->GetFiberWeights());

          m_manStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
          m_manStreamline->GetFiberPolyData()->RemoveDeletedCells();
          if (m_predlabeling==false)
          {
          vtkSmartPointer<vtkPolyData> vNewPolyData2 = vtkSmartPointer<vtkPolyData>::New();
          vtkSmartPointer<vtkCellArray> vNewLines2 = vtkSmartPointer<vtkCellArray>::New();
          vtkSmartPointer<vtkPoints> vNewPoints2 = vtkSmartPointer<vtkPoints>::New();


         counter = 0;
          for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
          {

            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
            for (unsigned int j=0; j<numPoints; j++)
            {
              double p[3];
              points->GetPoint(j, p);

              vtkIdType id = vNewPoints2->InsertNextPoint(p);
              container->GetPointIds()->InsertNextId(id);
            }
        //    weights->InsertValue(counter, fib->GetFiberWeight(i));
            vNewLines2->InsertNextCell(container);
            counter++;

          }
          m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints2);
          m_manStreamline->GetFiberPolyData()->SetLines(vNewLines2);
          m_manStreamline->SetFiberColors(255, 255, 255);
          }





    }
    }
    else
    {
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

        vtkIdType pickedCellID = picker->GetCellId();



        if (picker->GetCellId()==-1)
        {

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
//            m_PosStreamline->SetFiberWeights(m_PosStreamline->GetFiberWeights());


            m_manStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
            m_manStreamline->GetFiberPolyData()->RemoveDeletedCells();

            if (m_predlabeling==false)
            {
            vtkSmartPointer<vtkPolyData> vNewPolyData2 = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkCellArray> vNewLines2 = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> vNewPoints2 = vtkSmartPointer<vtkPoints>::New();


           counter = 0;
            for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
            {

              vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
              auto numPoints = cell->GetNumberOfPoints();
              vtkPoints* points = cell->GetPoints();

              vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
              for (unsigned int j=0; j<numPoints; j++)
              {
                double p[3];
                points->GetPoint(j, p);

                vtkIdType id = vNewPoints2->InsertNextPoint(p);
                container->GetPointIds()->InsertNextId(id);
              }
          //    weights->InsertValue(counter, fib->GetFiberWeight(i));
              vNewLines2->InsertNextCell(container);
              counter++;

            }
            m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints2);
            m_manStreamline->GetFiberPolyData()->SetLines(vNewLines2);
            m_manStreamline->SetFiberColors(255, 255, 255);
            }

        }
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
  }

void mitk::StreamlineInteractorBrush::AddStreamlineNegBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{

    DataInteractor::SetDataNode(m_manStreamlineNode);

    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {

    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);
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
//        }

      auto displayPosition = positionEvent->GetPointerPositionOnScreen();
      picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

      vtkIdType pickedCellID = picker->GetCellId();



      if (picker->GetCellId()==-1)
      {

      }
      else {


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
//            m_NegStreamline->SetFiberWeights(m_NegStreamline->GetFiberWeights());

          m_manStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
          m_manStreamline->GetFiberPolyData()->RemoveDeletedCells();

          if (m_predlabeling==false)
          {
          vtkSmartPointer<vtkPolyData> vNewPolyData2 = vtkSmartPointer<vtkPolyData>::New();
          vtkSmartPointer<vtkCellArray> vNewLines2 = vtkSmartPointer<vtkCellArray>::New();
          vtkSmartPointer<vtkPoints> vNewPoints2 = vtkSmartPointer<vtkPoints>::New();


         counter = 0;
          for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
          {

            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
            for (unsigned int j=0; j<numPoints; j++)
            {
              double p[3];
              points->GetPoint(j, p);

              vtkIdType id = vNewPoints2->InsertNextPoint(p);
              container->GetPointIds()->InsertNextId(id);
            }
        //    weights->InsertValue(counter, fib->GetFiberWeight(i));
            vNewLines2->InsertNextCell(container);
            counter++;

          }
          m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints2);
          m_manStreamline->GetFiberPolyData()->SetLines(vNewLines2);
          m_manStreamline->SetFiberColors(255, 255, 255);
            }





    }
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

            if (m_predlabeling==false)
            {
                vtkSmartPointer<vtkPolyData> vNewPolyData2 = vtkSmartPointer<vtkPolyData>::New();
                vtkSmartPointer<vtkCellArray> vNewLines2 = vtkSmartPointer<vtkCellArray>::New();
                vtkSmartPointer<vtkPoints> vNewPoints2 = vtkSmartPointer<vtkPoints>::New();


                counter = 0;
                for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
                {

                vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
                auto numPoints = cell->GetNumberOfPoints();
                vtkPoints* points = cell->GetPoints();

                vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
                for (unsigned int j=0; j<numPoints; j++)
                {
                  double p[3];
                  points->GetPoint(j, p);

                  vtkIdType id = vNewPoints2->InsertNextPoint(p);
                  container->GetPointIds()->InsertNextId(id);
                }
                //    weights->InsertValue(counter, fib->GetFiberWeight(i));
                vNewLines2->InsertNextCell(container);
                counter++;

                }
                m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints2);
                m_manStreamline->GetFiberPolyData()->SetLines(vNewLines2);
                m_manStreamline->SetFiberColors(255, 255, 255);
            }

        }
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
  }

void mitk::StreamlineInteractorBrush::AddNegStreamlinetolabelsBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    DataInteractor::SetDataNode(m_NegStreamlineNode);


    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {
    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);

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
//        }

      auto displayPosition = positionEvent->GetPointerPositionOnScreen();
      picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

      vtkIdType pickedCellID = picker->GetCellId();



      if (picker->GetCellId()==-1)
      {

      }
      else {


          vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
          vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
          vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

          unsigned int counter = 0;
          for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
          {

            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
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



          vtkCell* cell = m_NegStreamline->GetFiberPolyData()->GetCell(pickedCellID);

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

  //        m_manStreamline = mitk::FiberBundle::New(vNewPolyData);
          m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
          m_manStreamline->GetFiberPolyData()->SetLines(vNewLines);
          m_manStreamline->SetFiberColors(255, 255, 255);
//            m_manStreamline->SetFiberWeights(m_manStreamline->GetFiberWeights());

          m_NegStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
          m_NegStreamline->GetFiberPolyData()->RemoveDeletedCells();





    }
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

        }
        else
        {
            vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

            unsigned int counter = 0;
            for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
            {
              vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
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



            vtkCell* cell = m_NegStreamline->GetFiberPolyData()->GetCell(pickedCellID);
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

    //        m_manStreamline = mitk::FiberBundle::New(vNewPolyData);
            m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
            m_manStreamline->GetFiberPolyData()->SetLines(vNewLines);
            m_manStreamline->SetFiberColors(255, 255, 255);

            m_NegStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
            m_NegStreamline->GetFiberPolyData()->RemoveDeletedCells();

        }
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
  }

void mitk::StreamlineInteractorBrush::AddPosStreamlinetolabelsBundle(StateMachineAction *, InteractionEvent *interactionEvent)
{
    DataInteractor::SetDataNode(m_PosStreamlineNode);


    auto positionEvent = dynamic_cast<const InteractionPositionEvent *>(interactionEvent);
    if (positionEvent == nullptr)
    {

    }

    if (interactionEvent->GetSender()->GetMapperID() == BaseRenderer::Standard2D)
    {
//      m_PickedHandle = PickFrom2D(positionEvent);


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
//        }

      auto displayPosition = positionEvent->GetPointerPositionOnScreen();
      picker->Pick(displayPosition[0], displayPosition[1], 0, positionEvent->GetSender()->GetVtkRenderer());

      vtkIdType pickedCellID = picker->GetCellId();



      if (picker->GetCellId()==-1)
      {

      }
      else {


          vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
          vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
          vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

          unsigned int counter = 0;
          for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
          {

            vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
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



          vtkCell* cell = m_PosStreamline->GetFiberPolyData()->GetCell(pickedCellID);

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

  //        m_manStreamline = mitk::FiberBundle::New(vNewPolyData);
          m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
          m_manStreamline->GetFiberPolyData()->SetLines(vNewLines);
          m_manStreamline->SetFiberColors(255, 255, 255);
//            m_manStreamline->SetFiberWeights(m_manStreamline->GetFiberWeights());

          m_PosStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
          m_PosStreamline->GetFiberPolyData()->RemoveDeletedCells();





    }
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

        }
        else
        {
            vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

            unsigned int counter = 0;
            for ( int i=0; i<m_manStreamline->GetFiberPolyData()->GetNumberOfCells(); i++)
            {
              vtkCell* cell = m_manStreamline->GetFiberPolyData()->GetCell(i);
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



            vtkCell* cell = m_PosStreamline->GetFiberPolyData()->GetCell(pickedCellID);
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

    //        m_manStreamline = mitk::FiberBundle::New(vNewPolyData);
            m_manStreamline->GetFiberPolyData()->SetPoints(vNewPoints);
            m_manStreamline->GetFiberPolyData()->SetLines(vNewLines);
            m_manStreamline->SetFiberColors(255, 255, 255);

            m_PosStreamline->GetFiberPolyData()->DeleteCell(pickedCellID);
            m_PosStreamline->GetFiberPolyData()->RemoveDeletedCells();

        }
    }
    RenderingManager::GetInstance()->RequestUpdateAll();
  }

