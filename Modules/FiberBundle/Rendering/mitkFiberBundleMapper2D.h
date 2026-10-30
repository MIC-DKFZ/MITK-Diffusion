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
#include "mitkFiberSlabIndex.h"

#include <vtkSmartPointer.h>

#include <memory>

class vtkActor;
class vtkOpenGLPolyDataMapper;
class vtkPlane;

namespace mitk {

/**
 * \brief 2D mapper of fiber bundles: the fibers within Fiber2DSliceThickness of the shown plane,
 * as lines of shape.linewidth with the bundle's colours.
 *
 * Each window holds only the fibers in a band around its slab (FiberSlabIndex) and takes a new
 * band when the slab leaves it or turns. Two clipping planes cut out the slab, so the image is the
 * same as from the whole bundle.
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
    // The two planes of the slab, kept in the mapper's plane collection and only moved on a slice
    // step: a new collection or input changes the mapper's modified time, and VTK then rebuilds
    // the GPU buffers of the whole input (#46).
    vtkSmartPointer<vtkPlane> m_SlabPlane1;
    vtkSmartPointer<vtkPlane> m_SlabPlane2;
    // time the input was last built from the bundle
    itk::TimeStamp m_LastUpdateTime;
    // the input: the band of points p with m_BandFrom <= m_BandNormal * p <= m_BandTo, or the
    // whole bundle without a band
    bool m_HasBand = false;
    double m_BandNormal[3] = {0, 0, 0};
    double m_BandFrom = 0;
    double m_BandTo = 0;
    FBXLocalStorage();
  };

  /** \brief This member holds all three LocalStorages for the three 2D render windows. */
  mitk::LocalStorageHandler<FBXLocalStorage> m_LocalStorageHandler;

protected:
  FiberBundleMapper2D() = default;
  ~FiberBundleMapper2D() override = default;

  /** Gives the window the band around its slab, or the whole bundle (see UpdateSlab) */
  void GenerateDataForRenderer(mitk::BaseRenderer*) override;
  /** Moves the slab planes to the renderer's current plane and Fiber2DSliceThickness; gives the
   *  slab as the points p with from <= normal * p <= to (unit normal). Returns false without a
   *  plane. */
  bool UpdateSlab(mitk::BaseRenderer* renderer, FBXLocalStorage* localStorage, double normal[3], double& from, double& to);

private:
  // shared by the 2D windows, built from the bundle at m_SlabIndexTime
  std::unique_ptr<FiberSlabIndex> m_SlabIndex;
  itk::TimeStamp m_SlabIndexTime;
};

}

#endif
