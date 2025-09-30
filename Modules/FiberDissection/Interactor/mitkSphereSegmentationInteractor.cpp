/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkSphereSegmentationInteractor.h"
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
#include <mitkSurfaceToImageFilter.h>

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
#include <mitkLabelSetImage.h>

#include <vtkSphereSource.h>

mitk::SphereSegmentationInteractor::SphereSegmentationInteractor()
{
    // TODO if we want to get this configurable, the this is the recipe:
    // - make the 2D mapper add corresponding properties to control "enabled" and "color"
    // - make the interactor evaluate those properties
    // - in an ideal world, modify the state machine on the fly and skip mouse move handling
}



mitk::SphereSegmentationInteractor::~SphereSegmentationInteractor()
{
}



void mitk::SphereSegmentationInteractor::SetImage(mitk::Image::Pointer image)
{
    m_Image = image;
}

void mitk::SphereSegmentationInteractor::StartEndSurfaceNodes(mitk::DataNode::Pointer startSurfaceDataNode, mitk::DataNode::Pointer endSurfaceDataNode)
{
    m_startDataSurfaceNode = startSurfaceDataNode;
    m_endDataSurfaceNode = endSurfaceDataNode;

}


void mitk::SphereSegmentationInteractor::EndCreationStart(mitk::StateMachineAction*, mitk::InteractionEvent* /*interactionEvent*/)
{

  this->GetDataNode()->SetBoolProperty(DATANODE_PROPERTY_CREATED, true);
   if (this->GetDataNode()->GetName() == m_startDataNode->GetName()) {
     DataInteractor::SetDataNode(m_endDataNode);
   }
  else {
      CreateSegmentation();
  }

}

void mitk::SphereSegmentationInteractor::CreateSegmentation()
{
    MITK_INFO << "Start Create Segmentation";

    // Check if the required objects exist
    if (!m_Image || !m_startDataNode || !m_endDataNode) {
      MITK_ERROR << "Required objects are missing.";
      return;
    }

    mitk::Surface::Pointer startSurfaceData = dynamic_cast<mitk::Surface*>(m_startDataNode->GetData());
    mitk::Surface::Pointer endSurfaceData = dynamic_cast<mitk::Surface*>(m_endDataNode->GetData());

    // Check if the required surfaces exist
    if (!startSurfaceData || !endSurfaceData) {
      MITK_ERROR << "Required surface data is missing.";
      return;
    }

    mitk::SurfaceToImageFilter::Pointer startSurfaceToImageFilter = mitk::SurfaceToImageFilter::New();
    startSurfaceToImageFilter->MakeOutputBinaryOn();
    startSurfaceToImageFilter->SetInput(startSurfaceData);
    startSurfaceToImageFilter->SetImage(m_Image);
    startSurfaceToImageFilter->Update();

    mitk::Image::Pointer startImage = startSurfaceToImageFilter->GetOutput();
    if (startImage.IsNull()) {
      MITK_ERROR << "Convert Start Surface to binary image failed";
      return;
    }

    // Create name for the result node
    std::string nameOfResultImage = m_startDataNode->GetName() + "_ROI";
    m_startDataSurfaceNode->SetData(startImage);
    m_startDataSurfaceNode->SetProperty("name", mitk::StringProperty::New(nameOfResultImage));

    mitk::SurfaceToImageFilter::Pointer endSurfaceToImageFilter = mitk::SurfaceToImageFilter::New();
    endSurfaceToImageFilter->MakeOutputBinaryOn();
    endSurfaceToImageFilter->SetInput(endSurfaceData);
    endSurfaceToImageFilter->SetImage(m_Image);
    endSurfaceToImageFilter->Update();

    mitk::Image::Pointer endImage = endSurfaceToImageFilter->GetOutput();
    if (endImage.IsNull()) {
      MITK_ERROR << "Convert End Surface to binary image failed";
      return;
    }

    // Create name for the result node
    nameOfResultImage = m_endDataNode->GetName() + "_ROI";
    m_endDataSurfaceNode->SetData(endImage);
    m_endDataSurfaceNode->SetProperty("name", mitk::StringProperty::New(nameOfResultImage));

    // Hide the original nodes
    m_startDataNode->SetVisibility(false);
    m_endDataNode->SetVisibility(false);


}







