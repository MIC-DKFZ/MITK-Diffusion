/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPeakImageMapper3D_h
#define mitkPeakImageMapper3D_h

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
   * \brief 3D mapper of peak images, active while the node property "Enable3DPeaks" is on.
   *
   * With a "3DClipping" property of non-zero normal, only the peaks within half the smallest
   * voxel spacing of that plane are drawn.
   */
  class MITKDIFFUSIONMODELLING_EXPORT PeakImageMapper3D : public VtkMapper
  {
  public:
    mitkClassMacro(PeakImageMapper3D, VtkMapper);
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
    PeakImageMapper3D() = default;
    ~PeakImageMapper3D() override = default;

    void GenerateDataForRenderer(BaseRenderer *renderer) override;

  private:
    LocalStorageHandler<LocalStorage> m_LocalStorageHandler;
  };
}

#endif
