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

#include <MitkDiffusionCoreExports.h>

namespace mitk
{

  /**
  * \brief this class encapsulates tensor images
  */
  class MITKDIFFUSIONCORE_EXPORT PeakImage : public Image
  {

  public:
    typedef float PixelType;
    typedef itk::Image<PixelType, 4> ItkPeakImageType;

    mitkClassMacro( PeakImage, Image )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    void ConstructPolydata();

    vtkSmartPointer< vtkPolyData > GetPolyData()
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
    ~PeakImage() override;

    mutable vtkSmartPointer< vtkPolyData > m_PolyData;

  };

} // namespace mitk

#endif /* __mitkPeakImage__h */
