/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleMapper2D_h
#define mitkFiberBundleMapper2D_h

#include <MitkFiberBundleExports.h>
#include <mitkFiberBundle.h>
#include <mitkVtkMapper.h>

#include <vtkSmartPointer.h>

class vtkActor;
class vtkOpenGLPolyDataMapper;

namespace mitk {

/**
 * \brief 2D mapper of fiber bundles: the fibers within Fiber2DSliceThickness of the shown plane,
 * as lines of shape.linewidth with the bundle's colours.
 */
class MITKFIBERBUNDLE_EXPORT FiberBundleMapper2D : public VtkMapper
{

public:
  mitkClassMacro(FiberBundleMapper2D, VtkMapper);
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  mitk::FiberBundle* GetInput();
  void Update(mitk::BaseRenderer * renderer) override;
  static void SetDefaultProperties(DataNode* node, BaseRenderer* renderer = nullptr, bool overwrite = false );
  vtkProp* GetVtkProp(mitk::BaseRenderer* renderer) override;

  class  FBXLocalStorage : public mitk::Mapper::BaseLocalStorage
  {
  public:
    vtkSmartPointer<vtkActor> m_Actor;
    vtkSmartPointer<vtkOpenGLPolyDataMapper> m_Mapper;
    itk::TimeStamp m_LastUpdateTime;
    // line width of the last update; the actor of each render window has its own
    int m_LineWidth = -1;
    FBXLocalStorage();
  };

  /** \brief This member holds all three LocalStorages for the three 2D render windows. */
  mitk::LocalStorageHandler<FBXLocalStorage> m_LocalStorageHandler;

protected:
  FiberBundleMapper2D() = default;
  ~FiberBundleMapper2D() override = default;

  void GenerateDataForRenderer(mitk::BaseRenderer*) override;
};

}

#endif
