/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkPeakImage_h
#define mitkPeakImage_h

#include "mitkImage.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <MitkDiffusionModellingExports.h>

#include <array>
#include <mutex>
#include <optional>

namespace mitk
{

  /**
  * \brief Peak image: a 4D float image whose fourth dimension holds three components (x, y, z in
  * world coordinates) per peak and voxel.
  */
  class MITKDIFFUSIONMODELLING_EXPORT PeakImage : public Image
  {

  public:
    typedef float PixelType;
    typedef itk::Image<PixelType, 4> ItkPeakImageType;

    mitkClassMacro( PeakImage, Image )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    mitkCloneMacro(Self)

    /**
    * One two-point line per peak longer than 0.001, centred on its voxel in world coordinates,
    * with the "FIBER_COLORS" point array. Built on first use and rebuilt when the image or its
    * geometry was modified since; a rebuilt polydata keeps the colouring last chosen.
    */
    vtkSmartPointer< vtkPolyData > GetPolyData() const;

    ItkPeakImageType::Pointer GetItkPeakImage();

    void SetCustomColor(float r, float g, float b);
    void ColorByOrientation();

  protected:
    PeakImage();
    PeakImage(const PeakImage& other);
    ~PeakImage() override;

  private:
    void ConstructPolydata() const;
    /** Records that the polydata is up to date after a colour change, which calls Modified(). */
    void KeepPolyData();

    mutable vtkSmartPointer< vtkPolyData > m_PolyData;
    /** Modification time of this image and its geometry when m_PolyData was built or coloured. */
    mutable itk::ModifiedTimeType m_PolyDataTime = 0;
    mutable std::mutex m_PolyDataMutex;
    /** Colour of SetCustomColor(); none colours by orientation. */
    std::optional<std::array<float, 3>> m_CustomColor;

  };

} // namespace mitk

#endif
