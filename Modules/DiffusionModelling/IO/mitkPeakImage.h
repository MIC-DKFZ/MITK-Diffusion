/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef __mitkPeakImage__h
#define __mitkPeakImage__h

#include "mitkImage.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <MitkDiffusionModellingExports.h>

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
    * One two-point line per peak longer than 0.001, centred on its voxel, with the "FIBER_COLORS"
    * point array; built on first use.
    */
    vtkSmartPointer< vtkPolyData > GetPolyData() const
    {
      if (m_PolyData==nullptr)
        this->ConstructPolydata();
      return m_PolyData;
    }

    ItkPeakImageType::Pointer GetItkPeakImage();

    void SetCustomColor(float r, float g, float b);
    void ColorByOrientation();

  protected:
    PeakImage();
    PeakImage(const PeakImage& other);
    ~PeakImage() override;

    mutable vtkSmartPointer< vtkPolyData > m_PolyData;

  private:
    void ConstructPolydata() const;
    void ColorPolyDataByOrientation() const;

  };

} // namespace mitk

#endif /* __mitkPeakImage__h */
