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
#include "mitkFiberSlabTiles.h"

#include <vtkSmartPointer.h>

#include <memory>

class vtkCellPicker;
class vtkPolyData;
class vtkUnsignedCharArray;

namespace mitk {

/**
 * \brief 2D mapper of fiber bundles: the fibers within Fiber2DSliceThickness of the shown plane,
 * as lines of shape.linewidth with the bundle's colours.
 *
 * Each window holds only the tiles of fibers near its slab (FiberSlabTiles) and takes new tiles
 * when the slab moves or turns. Two clipping planes cut out the slab, so the image is the same as
 * from the whole bundle. New colours keep the tiles' geometry.
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

  /** The index of the fiber whose cell the picker picked on a fiber mapper's prop, -1 if none.
   *  The cells of the 2D input are pieces of fibers, so the picked cell id is not the fiber
   *  index there; the 3D input holds the bundle's lines, whose cell id is the fiber index. */
  static vtkIdType PickedFiber(vtkCellPicker* picker);

  class  FBXLocalStorage : public mitk::Mapper::BaseLocalStorage
  {
  public:
    FiberSlabTiles m_Tiles;
    // the versions of the mapper's index and colours the tiles were built or coloured from
    unsigned long m_IndexVersion = 0;
    unsigned long m_ColourVersion = 0;
  };

  /** \brief This member holds all three LocalStorages for the three 2D render windows. */
  mitk::LocalStorageHandler<FBXLocalStorage> m_LocalStorageHandler;

protected:
  FiberBundleMapper2D() = default;
  ~FiberBundleMapper2D() override = default;

  /** Moves the slab planes to the renderer's current plane and Fiber2DSliceThickness; gives the
   *  slab as the points p with from <= normal * p <= to (unit normal). Returns false without a
   *  plane. */
  bool UpdateSlab(mitk::BaseRenderer* renderer, FBXLocalStorage* localStorage, double normal[3], double& from, double& to);
  /** Rebuilds the index for new fiber geometry and takes new colours; false if the bundle's
   *  colours do not fit its points */
  bool UpdateFibers(mitk::FiberBundle* fiberBundle);

private:
  // shared by the 2D windows; geometry only, built from m_IndexedFibers at its modified time
  std::unique_ptr<FiberSlabIndex> m_SlabIndex;
  vtkPolyData* m_IndexedFibers = nullptr;
  vtkMTimeType m_IndexedFibersTime = 0;
  unsigned long m_IndexVersion = 0;
  // the bundle's colours (an RGBA tuple per point) as of the bundle's 2D update time m_ColourTime
  vtkSmartPointer<vtkUnsignedCharArray> m_Colours;
  unsigned char m_LowestAlpha = 255;
  itk::ModifiedTimeType m_ColourTime = 0;
  unsigned long m_ColourVersion = 0;
};

}

#endif
