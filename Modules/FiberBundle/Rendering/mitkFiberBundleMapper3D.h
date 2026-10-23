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
class vtkPropAssembly;

namespace mitk {

/**
 * \brief 3D mapper of fiber bundles: lines, tubes (shape.tuberadius) or ribbons
 * (shape.ribbonwidth) with the bundle's colours, optionally clipped at the plane of the
 * 3DClipping property.
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

  /** The node properties that change the geometry the mapper builds. */
  struct ShapeParameters
  {
    float tubeRadius = 0;
    int tubeSides = 0;
    int lineWidth = 0;
    float ribbonWidth = 0;
    bool operator==(const ShapeParameters&) const = default;
  };

  class  LocalStorage3D : public mitk::Mapper::BaseLocalStorage
  {
  public:
    vtkSmartPointer<vtkActor> m_FiberActor;
    vtkSmartPointer<vtkOpenGLPolyDataMapper> m_FiberMapper;
    vtkSmartPointer<vtkPropAssembly> m_FiberAssembly;

    itk::TimeStamp m_LastUpdateTime;
    // shape of the last update; each render window builds its own
    ShapeParameters m_Shape;
    LocalStorage3D();
  };

  /** \brief This member holds all three LocalStorages for the 3D render window(s). */
  mitk::LocalStorageHandler<LocalStorage3D> m_LocalStorageHandler;


protected:

  FiberBundleMapper3D();
  ~FiberBundleMapper3D() override = default;
  void InternalGenerateData(mitk::BaseRenderer *renderer, const FiberBundle* fiberBundle, const ShapeParameters& shape);

private:
  vtkSmartPointer<vtkLookupTable> m_lut;
};

}

#endif
