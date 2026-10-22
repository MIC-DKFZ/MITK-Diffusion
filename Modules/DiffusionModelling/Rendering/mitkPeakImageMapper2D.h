/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPeakImageMapper2D_h
#define mitkPeakImageMapper2D_h

#include <MitkDiffusionModellingExports.h>

#include <mitkLocalStorageHandler.h>
#include <mitkPeakImage.h>
#include <mitkVtkMapper.h>

#include <vtkActor.h>
#include <vtkPlane.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

namespace mitk
{
  /**
   * \brief 2D mapper of peak images: the peak lines of the voxel layer the renderer's plane cuts.
   *
   * The peak polydata is in world coordinates; two clipping planes keep a slab of half the
   * smallest voxel spacing on each side of the plane.
   */
  class MITKDIFFUSIONMODELLING_EXPORT PeakImageMapper2D : public VtkMapper
  {
  public:
    mitkClassMacro(PeakImageMapper2D, VtkMapper);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

    PeakImage *GetInput();

    static void SetDefaultProperties(DataNode *node, BaseRenderer *renderer = nullptr, bool overwrite = false);

    /** The peak polydata is already in world coordinates. */
    void UpdateVtkTransform(BaseRenderer *) override {}
    vtkProp *GetVtkProp(BaseRenderer *renderer) override;

    class LocalStorage : public Mapper::BaseLocalStorage
    {
    public:
      LocalStorage();

      vtkSmartPointer<vtkActor> m_Actor;
      vtkSmartPointer<vtkPolyDataMapper> m_Mapper;
      vtkSmartPointer<vtkPlane> m_LowerPlane;
      vtkSmartPointer<vtkPlane> m_UpperPlane;
    };

  protected:
    PeakImageMapper2D() = default;
    ~PeakImageMapper2D() override = default;

    void GenerateDataForRenderer(BaseRenderer *renderer) override;

  private:
    LocalStorageHandler<LocalStorage> m_LocalStorageHandler;
  };
}

#endif
