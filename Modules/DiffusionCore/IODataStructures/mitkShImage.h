/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical Image Computing.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef __mitkShImage__h
#define __mitkShImage__h

#include "mitkImage.h"
#include "itkVectorImage.h"
#include "mitkImageVtkAccessor.h"

#include <MitkDiffusionCoreExports.h>

namespace mitk
{

  /**
  * \brief this class encapsulates spherical harmonics coefficients
  */
  class MITKDIFFUSIONCORE_EXPORT ShImage : public Image
  {

  public:

    enum SH_CONVENTION {  ///< SH coefficient convention (depends on toolkit)
        FSL,    // used in FSL, Dipy
        MRTRIX  // used in MRtrix, MITK
    };

    typedef itk::Image<float,4 >  ShOnDiskType; // we store the sh images in MRtrix 4D float format and convert on load to 3D multi-component images

    mitkClassMacro( ShImage, Image )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    virtual const vtkImageData* GetNonRgbVtkImageData(int t = 0, int n = 0) const;
    virtual vtkImageData* GetNonRgbVtkImageData(int t = 0, int n = 0);

    const vtkImageData* GetVtkImageData(int t = 0, int n = 0) const override;
    vtkImageData* GetVtkImageData(int t = 0, int n = 0) override;

    virtual void ConstructRgbImage() const;

    unsigned int ShOrder();
    unsigned int NumCoefficients();

    SH_CONVENTION GetShConvention() const;
    void SetShConvention(SH_CONVENTION ShConvention);

  protected:
    ShImage();
    ~ShImage() override;

    template<int num_components>
    void Construct() const;

    void  PrintSelf(std::ostream &os, itk::Indent indent) const override;

    mutable mitk::Image::Pointer m_RgbImage;
    unsigned int m_ShOrder;
    unsigned int m_NumCoefficients;
    SH_CONVENTION m_ShConvention;  // use mrtrix style SH convention
  };

} // namespace mitk

#endif /* __mitkShImage__h */
