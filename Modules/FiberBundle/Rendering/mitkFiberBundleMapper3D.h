/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkFiberBundleMapper3D_h
#define mitkFiberBundleMapper3D_h

#include <MitkFiberBundleExports.h>
#include <mitkFiberBundle.h>
#include <mitkVtkMapper.h>

#include <vtkSmartPointer.h>

class vtkActor;
class vtkLookupTable;
class vtkOpenGLPolyDataMapper;
class vtkPlane;
class vtkPropAssembly;

namespace mitk {

/**
 * \brief 3D mapper of fiber bundles: lines, tubes (shape.tuberadius, in mm) or ribbons
 * (shape.ribbonwidth) with the bundle's colours, optionally clipped at the plane of the
 * 3DClipping property.
 *
 * Tubes are drawn on the GPU from the lines (mitkFiberTubeShader.h); ribbons are built on the
 * CPU with vtkRibbonFilter, whose band orientation follows each whole fiber.
 */
class MITKFIBERBUNDLE_EXPORT FiberBundleMapper3D : public VtkMapper
{
public:

  mitkClassMacro(FiberBundleMapper3D, VtkMapper)
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  vtkProp *GetVtkProp(mitk::BaseRenderer *renderer) override;
  static void SetDefaultProperties(DataNode* node, BaseRenderer* renderer = nullptr, bool overwrite = false );
  void GenerateDataForRenderer(mitk::BaseRenderer* renderer) override;

  class  LocalStorage3D : public mitk::Mapper::BaseLocalStorage
  {
  public:
    vtkSmartPointer<vtkActor> m_FiberActor;
    vtkSmartPointer<vtkOpenGLPolyDataMapper> m_FiberMapper;
    vtkSmartPointer<vtkPropAssembly> m_FiberAssembly;

    // time the input was last built from the bundle
    itk::TimeStamp m_LastUpdateTime;
    // ribbon width of the input; 0: the input is the lines (also for tubes)
    float m_RibbonWidth = 0;
    // The plane of 3DClipping, in the fiber mapper's plane collection while its normal is not
    // zero, and moved in place: a new plane changes the mapper's modified time, and VTK then
    // rebuilds the GPU buffers of the whole bundle (#46).
    vtkSmartPointer<vtkPlane> m_ClippingPlane;
    LocalStorage3D();
  };

  /** \brief This member holds all three LocalStorages for the 3D render window(s). */
  mitk::LocalStorageHandler<LocalStorage3D> m_LocalStorageHandler;


protected:

  FiberBundleMapper3D();
  ~FiberBundleMapper3D() override = default;
  /** Gives the fiber mapper a new input: the lines, or ribbons for ribbonWidth > 0 */
  void InternalGenerateData(mitk::BaseRenderer *renderer, const FiberBundle* fiberBundle, float ribbonWidth);
  /** Applies the plane of the node's 3DClipping property to the fiber mapper when it changed. */
  void UpdateClippingPlane(LocalStorage3D* localStorage);

private:
  vtkSmartPointer<vtkLookupTable> m_lut;
};

}

#endif
