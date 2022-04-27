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

// VTK includes
#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>

// System includes
#include <memory>

#include "MitkStreamlineExports.h"

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
  class MITKStreamline_EXPORT StreamlineInteractor : public DataInteractor
  {
  public:
    mitkClassMacro(StreamlineInteractor, DataInteractor);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

      //! The node holding the Streamline for visual feedback.
      //! This is the node that the interactor is primarily working on
      //! (calls DataInteractor::SetDataNode).
      void SetStreamlineNode(DataNode *node);

  private:
    StreamlineInteractor();
    ~StreamlineInteractor() override;

    //! Setup the relation between the XML state machine and this object's methods.
    void ConnectActionsAndFunctions() override;

    //! State machine condition: successful Streamline picking
    //! \return true when any part of the Streamline has been picked.
    bool HasPickedHandle(const InteractionEvent *);

    //! Pick a Streamline handle from a 2D event (passing by the 2D mapper)
    Streamline::HandleType PickFrom2D(const InteractionPositionEvent *positionEvent);

    //! Pick a Streamline handle from a 3D event
    //! (passing by the general surface mapper and the Streamline object)
    Streamline::HandleType PickFrom3D(const InteractionPositionEvent *positionEvent);

    void UpdateHandleHighlight();

    //! the Streamline used for visual feedback and picking
    Streamline::Pointer m_Streamline;

    //! The manipulated object's geometry
    BaseGeometry::Pointer m_ManipulatedObjectGeometry;

    //! For picking on the vtkPolyData representing the Streamline
    std::map<BaseRenderer *, vtkSmartPointer<vtkCellPicker>> m_Picker;

    //! Part of the Streamline that was picked on last check
    Streamline::HandleType m_PickedHandle = Streamline::NoHandle;

    //! Part of the Streamline that is currently highlighted
    Streamline::HandleType m_HighlightedHandle = Streamline::NoHandle;

    //! Color (RGBA) used for highlighting
    double m_ColorForHighlight[4];

    //! Color (RGBA) that has been replaced by m_ColorForHighlight
    double m_ColorReplacedByHighlight[4];

    Point2D m_InitialClickPosition2D; //< Initial screen click position
    double m_InitialClickPosition2DZ; //< Z value of the initial screen click position
    Point3D m_InitialClickPosition3D; //< Initial 3D click position

    Point2D m_InitialStreamlineCenter2D; //< Initial position of the Streamline's center in screen coordinates
    Point3D m_InitialStreamlineCenter3D; //< Initial 3D position of the Streamline's center

    Vector3D m_AxisOfMovement; //< Axis along which we move when translating
    Vector3D m_AxisOfRotation; //< Axis around which we turn when rotating

    std::unique_ptr<Operation> m_FinalDoOperation;   //< Operation for the undo-stack
    std::unique_ptr<Operation> m_FinalUndoOperation; //< Operation for the undo-stack

    //! A copy of the origin geometry, to avoid accumulation of tiny errors
    BaseGeometry::Pointer m_InitialManipulatedObjectGeometry;
  };
}
#endif
