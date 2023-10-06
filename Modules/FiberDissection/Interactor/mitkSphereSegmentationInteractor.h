/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkSphereSegmentationInteractor_h
#define mitkSphereSegmentationInteractor_h


#include "MitkFiberDissectionExports.h"

// MITK includes
#include <mitkSphereInteractor.h>
#include <mitkDataNode.h>
#include <mitkGeometry3D.h>
#include <mitkFiberBundle.h>

// VTK includes
#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>

// System includes
#include <memory>


namespace mitk
{
  class InteractionPositionEvent;

  //! Data interactor to pick Spheres via interaction
  //! with a mitk::Sphere.

  class MITKFIBERDISSECTION_EXPORT SphereSegmentationInteractor : public SphereInteractor
  {
  public:
      mitkClassMacro(SphereSegmentationInteractor, SphereInteractor);
      itkFactorylessNewMacro(Self);
      itkCloneMacro(Self);

      virtual void EndCreationStart(StateMachineAction*, InteractionEvent*) override;
      void StartEndSurfaceNodes(mitk::DataNode::Pointer startSurfaceDataNode, mitk::DataNode::Pointer endSurfaceDataNode);
      void SetImage(mitk::Image::Pointer image);


  protected:

    void CreateSegmentation();


    mitk::DataNode::Pointer m_startDataSurfaceNode;
    mitk::DataNode::Pointer m_endDataSurfaceNode;
    mitk::Image::Pointer m_Image;



  private:


    SphereSegmentationInteractor();
    ~SphereSegmentationInteractor() override;






  };
}
#endif
