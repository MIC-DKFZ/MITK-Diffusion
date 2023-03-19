/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkSphereInteractor.h"
//#include "mitkSphereMapper2D.h"

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
#include <vtkSelectEnclosedPoints.h>


#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "mitkFiberBundle.h"
#include "mitkInteractionPositionEvent.h"

#include <vtkSphereSource.h>

mitk::SphereInteractor::SphereInteractor()
{


    // TODO if we want to get this configurable, the this is the recipe:
    // - make the 2D mapper add corresponding properties to control "enabled" and "color"
    // - make the interactor evaluate those properties
    // - in an ideal world, modify the state machine on the fly and skip mouse move handling
}

const char* mitk::SphereInteractor::DATANODE_PROPERTY_SIZE = "zone.size";
const char* mitk::SphereInteractor::DATANODE_PROPERTY_CREATED = "zone.created";


void mitk::SphereInteractor::UpdateSurface(mitk::DataNode::Pointer dataNode)
{

  if (!dataNode->GetData())
  {
    MITK_WARN("SphereInteractor")("DataInteractor")
      << "Cannot update surface for node as no data is set to the node.";
    return;
  }

  mitk::Point3D origin = dataNode->GetData()->GetGeometry()->GetOrigin();

  float radius;
  if (!dataNode->GetFloatProperty(DATANODE_PROPERTY_SIZE, radius))
  {
    MITK_WARN("SphereInteractor")("DataInteractor")
      << "Cannut update surface for node as no radius is specified in the node properties.";
    return;
  }

  mitk::Surface::Pointer zone = mitk::Surface::New();

  // create a vtk sphere with given radius
  vtkSmartPointer<vtkSphereSource> vtkSphere = vtkSmartPointer<vtkSphereSource>::New();
  vtkSphere->SetRadius(radius);
  vtkSphere->SetCenter(0, 0, 0);
  vtkSphere->SetPhiResolution(20);
  vtkSphere->SetThetaResolution(20);
  vtkSphere->Update();
  zone->SetVtkPolyData(vtkSphere->GetOutput());

  // set vtk sphere and origin to data node (origin must be set
  // again, because of the new sphere set as data)
  dataNode->SetData(zone);
  dataNode->GetData()->GetGeometry()->SetOrigin(origin);

  // update the RenderWindow to show the changed surface
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();

}




mitk::SphereInteractor::~SphereInteractor()
{
}

void mitk::SphereInteractor::ConnectActionsAndFunctions()
{
  CONNECT_FUNCTION("addCenter", AddCenter);
  CONNECT_FUNCTION("changeRadius", ChangeRadius);
  CONNECT_FUNCTION("endCreationStart", EndCreationStart);
//  CONNECT_FUNCTION("endCreation", EndCreation);
  CONNECT_FUNCTION("abortCreation", AbortCreation);
}

void mitk::SphereInteractor::DataNodeChanged()
{
  mitk::DataNode::Pointer dataNode = this->GetDataNode();
  if (dataNode.IsNotNull() && dataNode->GetData() == nullptr)
  {
    dataNode->SetData(mitk::Surface::New());
  }
}

void mitk::SphereInteractor::StartEndNodes(mitk::DataNode::Pointer startDataNode, mitk::DataNode::Pointer endDataNode){
    m_startDataNode = startDataNode;
    m_endDataNode = endDataNode;

    DataInteractor::SetDataNode(startDataNode);
}

void mitk::SphereInteractor::workingBundleNode(mitk::FiberBundle::Pointer workingBundle, mitk::FiberBundle::Pointer reducedBundle){
    m_workingBundle = workingBundle;
    m_reducedFibersBundle = reducedBundle;

}

void mitk::SphereInteractor::AddCenter(mitk::StateMachineAction*, mitk::InteractionEvent* interactionEvent)
{
  // cast InteractionEvent to a position event in order to read out the mouse position
  mitk::InteractionPositionEvent* positionEvent = dynamic_cast<mitk::InteractionPositionEvent*>(interactionEvent);
  mitk::DataNode::Pointer dataNode = this->GetDataNode();
  dataNode->SetBoolProperty(DATANODE_PROPERTY_CREATED, false);

  // make sure that data node contains data
  mitk::BaseData::Pointer dataNodeData = this->GetDataNode()->GetData();
  if (dataNodeData.IsNull())
  {
    dataNodeData = mitk::Surface::New();
    this->GetDataNode()->SetData(dataNodeData);
  }

  // set origin of the data node to the mouse click position
  dataNodeData->GetGeometry()->SetOrigin(positionEvent->GetPositionInWorld());
  MITK_INFO << "Critical Structure added on position " << positionEvent->GetPointerPositionOnScreen() << " (Image Coordinates); " << positionEvent->GetPositionInWorld() << " (World Coordinates)";

  dataNode->SetFloatProperty("opacity", 0.60f);

//  return true;
}


void mitk::SphereInteractor::ChangeRadius(mitk::StateMachineAction*, mitk::InteractionEvent* interactionEvent)
{
  // cast InteractionEvent to a position event in order to read out the mouse position
  mitk::InteractionPositionEvent* positionEvent = dynamic_cast<mitk::InteractionPositionEvent*>(interactionEvent);
  mitk::DataNode::Pointer curNode = this->GetDataNode();
  mitk::Point3D mousePosition = positionEvent->GetPositionInWorld();

  mitk::ScalarType radius = mousePosition.EuclideanDistanceTo(curNode->GetData()->GetGeometry()->GetOrigin());

  curNode->SetFloatProperty(DATANODE_PROPERTY_SIZE, radius);

  mitk::SphereInteractor::UpdateSurface(curNode);

  //return true;
}

void mitk::SphereInteractor::EndCreationStart(mitk::StateMachineAction*, mitk::InteractionEvent* /*interactionEvent*/)
{
   MITK_INFO << "EndCreationStart";

  this->GetDataNode()->SetBoolProperty(DATANODE_PROPERTY_CREATED, true);
   if (this->GetDataNode()->GetName() == m_startDataNode->GetName()) {
     DataInteractor::SetDataNode(m_endDataNode);
   }
  else {
      ExtractFibers();
  }


  //return true;
}

//void mitk::SphereInteractor::EndCreation(mitk::StateMachineAction*, mitk::InteractionEvent* /*interactionEvent*/)
//{
//    MITK_INFO<< "Extract Fibers";
//   this->GetDataNode()->SetBoolProperty(DATANODE_PROPERTY_CREATED, true);
//   ExtractFibers();
//  //return true;
//}

void mitk::SphereInteractor::AbortCreation(mitk::StateMachineAction*, mitk::InteractionEvent*)
{
    MITK_INFO << "Inside Abort Creation";
  this->GetDataNode()->SetData(mitk::Surface::New());

  // update the RenderWindow to remove the surface
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();



  //return true;
}

void mitk::SphereInteractor::ExtractFibers()
{
    MITK_INFO << "Extract Fibers";


    vtkPolyData* polyData = m_workingBundle->GetFiberPolyData();

    vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

    unsigned int counter = 0;
    for (int i = 0; i < polyData->GetNumberOfCells(); i++)
    {
        vtkCell* cell = polyData->GetCell(i);

        // Assume you have a vtkCell object named cell
        vtkPoints* points = cell->GetPoints();

        double startPoint[3];
        double endPoint[3];

        int numPoints = cell->GetNumberOfPoints();

        points->GetPoint(0, startPoint);
        points->GetPoint(numPoints - 1, endPoint);

        mitk::Surface::Pointer startSurfaceData = dynamic_cast<mitk::Surface*>(m_startDataNode->GetData());
        mitk::Surface::Pointer endSurfaceData = dynamic_cast<mitk::Surface*>(m_endDataNode->GetData());



        if ((startSurfaceData->GetGeometry()->IsInside(startPoint) && endSurfaceData->GetGeometry()->IsInside(endPoint)) || ((endSurfaceData->GetGeometry()->IsInside(endPoint) && startSurfaceData->GetGeometry()->IsInside(startPoint))))
        {
            MITK_INFO << "Inside";
            // Create a new vtkPolyLine container to store the new fiber
              vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

              // Loop through each point in the fiber
              for (int j = 0; j < numPoints; j++) {
                  double p[3];
                  points->GetPoint(j, p);

                  // Insert the new point into the vtkPoints object and get its ID
                  vtkIdType id = vNewPoints->InsertNextPoint(p);
                  // Add the ID to the container object to create a new fiber
                  container->GetPointIds()->InsertNextId(id);
              }

              // Insert the fiber weight into the vtkDoubleArray object
//              weights->InsertValue(counter, m_workingBundle->GetFiberWeight(i));
              // Insert the new fiber into the vtkCellArray object
              vNewLines->InsertNextCell(container);

              // Increment the counter variable
              counter++;

        }
    }
        vNewPolyData->SetLines(vNewLines);
        vNewPolyData->SetPoints(vNewPoints);

        // Create a new vtkPolyData object and set it to the new fibers
        m_reducedFibersBundle->GetFiberPolyData()->SetPoints(vNewPoints);
        m_reducedFibersBundle->GetFiberPolyData()->SetLines(vNewLines);


        // Create a new mitk::FiberBundle object and set it to the new fibers
//        m_reducedFibersBundle = mitk::FiberBundle::New(vNewPolyData);
//        m_reducedFibersBundle->SetFiberColors(255, 255, 255);
        m_reducedFibersBundle->SetFiberWeights(weights);
        m_reducedFibersBundle->ColorFibersByOrientation();

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


