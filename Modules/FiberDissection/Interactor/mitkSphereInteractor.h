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

    void StartEndNodes(mitk::DataNode::Pointer startDataNode, mitk::DataNode::Pointer endDataNode);
    mitk::DataNode::Pointer m_startDataNode;
    mitk::DataNode::Pointer m_endDataNode;

  protected:

    void DataNodeChanged() override;
    /**
     * \brief Sets origin of the data node to the coordinates of the position event.
     * \return false if interaction event isn't a position event, true otherwise
     */
    void AddCenter(StateMachineAction*, InteractionEvent*);
    void ChangeRadius(StateMachineAction*, InteractionEvent*);
    void EndCreation(StateMachineAction*, InteractionEvent*);
    void EndCreationStart(StateMachineAction*, InteractionEvent*);
    void AbortCreation(StateMachineAction*, InteractionEvent*);
    /**
     * \brief Updates radius attribute according to position event.
     * Calculates distance between the data node origin and the current position
     * event and updates the data node attribue "zone.size" accordingly. The
     * mitk::USZonesInteractor::UpdateSurface() function is called then.
     *
     * \return false if interaction event isn't a position event, true otherwise
     */


  private:
    SphereInteractor();
    ~SphereInteractor() override;

    //! Setup the relation between the XML state machine and this object's methods.
    void ConnectActionsAndFunctions() override;






  };
}
#endif
