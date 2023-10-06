/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkSphereInteractor_h
#define mitkSphereInteractor_h


#include "MitkFiberDissectionExports.h"

// MITK includes
#include <mitkDataInteractor.h>
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

  class MITKFIBERDISSECTION_EXPORT SphereInteractor : public DataInteractor
  {
  public:
    mitkClassMacro(SphereInteractor, DataInteractor);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);



    static const char* DATANODE_PROPERTY_SIZE;
    static const char* DATANODE_PROPERTY_CREATED;


    static void UpdateSurface(itk::SmartPointer<mitk::DataNode>);
    void ExtractFibers();

    void workingBundleNode(mitk::DataNode::Pointer workingBundleNode, mitk::DataNode::Pointer reducedBundleNode);

    mitk::DataNode::Pointer m_startDataNode;
    mitk::DataNode::Pointer m_endDataNode;

    mitk::DataNode::Pointer m_workingBundleNode;
    mitk::DataNode::Pointer m_reducedFibersBundleNode;
    mitk::FiberBundle::Pointer m_workingBundle;
    mitk::FiberBundle::Pointer m_reducedFibersBundle;

    void DataNodeChanged() override;


    void StartEndNodes(mitk::DataNode::Pointer startDataNode, mitk::DataNode::Pointer endDataNode);

    void AddCenter(StateMachineAction*, InteractionEvent*);
    void ChangeRadius(StateMachineAction*, InteractionEvent*);
//    void EndCreation(StateMachineAction*, InteractionEvent*);
    virtual void EndCreationStart(StateMachineAction*, InteractionEvent*);
    void AbortCreation(StateMachineAction*, InteractionEvent*);

  protected:




    SphereInteractor();
    ~SphereInteractor() override;




  private:

    //! Setup the relation between the XML state machine and this object's methods.
    void ConnectActionsAndFunctions() override;






  };
}
#endif
