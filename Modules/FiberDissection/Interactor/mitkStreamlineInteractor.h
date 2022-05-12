/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkStreamlineInteractor_h
#define mitkStreamlineInteractor_h


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

#include "MitkFiberDissectionExports.h"

namespace mitk
{
  class InteractionPositionEvent;

  //! Data interactor to pick streamlines via interaction
  //! with a mitk::Streamline.
  //!
  //!
  //! To determine what parts of the object are clicked during interaction,
  //! the mappers (2D: custom mapper, 3D: regular surface mapper) are asked
  //! for their VTK objects, picking is performed, and the picked point is
  //! forwarded to the Streamline object for interpretation.
  //!
  //! The interactor fills the undo/redo stack with operations on the modified geometry.
  //!
  //! \sa Streamline
  class MITKFIBERDISSECTION_EXPORT StreamlineInteractor : public DataInteractor
  {
  public:
    mitkClassMacro(StreamlineInteractor, DataInteractor);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

      //! The node holding the Fiberbundle for visual feedback.
      //! This is the node that the interactor is primarily working on
      //! (calls DataInteractor::SetDataNode).

      void SetNegativeNode(DataNode *node);
      void SetToLabelNode(DataNode *node);
      void SetPositiveNode(DataNode *node);

  protected:

       void AddStreamlineNegBundle(StateMachineAction *, InteractionEvent *interactionEvent);

       void AddStreamlinePosBundle(StateMachineAction *, InteractionEvent *interactionEvent);


       std::map<BaseRenderer *, vtkSmartPointer<vtkCellPicker>> m_Picker;

  private:
    StreamlineInteractor();
    ~StreamlineInteractor() override;

    //! Setup the relation between the XML state machine and this object's methods.
    void ConnectActionsAndFunctions() override;

    //! State machine condition: successful Streamline picking
    //! \return true when any part of the Streamline has been picked.
//    bool HasPickedHandle(const InteractionEvent *);

//    void DecideInteraction(StateMachineAction *, InteractionEvent *interactionEvent);

//    //! Pick a Streamline handle from a 2D event (passing by the 2D mapper)
//    Streamline::HandleType PickFrom2D(const InteractionPositionEvent *positionEvent);

//    //! Pick a Streamline handle from a 3D event
//    //! (passing by the general surface mapper and the Streamline object)
//    Streamline::HandleType PickFrom3D(const InteractionPositionEvent *positionEvent);

//    void UpdateHandleHighlight();

    //! the Streamline used for visual feedback and picking
    mitk::FiberBundle::Pointer m_NegStreamline;
    mitk::FiberBundle::Pointer m_PosStreamline;
    mitk::FiberBundle::Pointer m_manStreamline;

    vtkSmartPointer<vtkPolyData>  m_extracted_streamline;

    double m_ColorForHighlight[4];

  };
}
#endif
