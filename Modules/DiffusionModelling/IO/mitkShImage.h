/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkShImage_h
#define mitkShImage_h

#include "mitkRgbCachedImage.h"

#include <MitkDiffusionModellingExports.h>

namespace mitk
{

  /**
  * \brief this class encapsulates spherical harmonics coefficients
  */
  class MITKDIFFUSIONMODELLING_EXPORT ShImage : public RgbCachedImage
  {

  public:

    enum SH_CONVENTION {  ///< SH coefficient convention (depends on toolkit)
        FSL,    // used in FSL, Dipy
        MRTRIX  // used in MRtrix, MITK
    };

    typedef itk::Image<float,4 >  ShOnDiskType; // we store the sh images in MRtrix 4D float format and convert on load to 3D multi-component images

    mitkClassMacro( ShImage, RgbCachedImage )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    mitkCloneMacro(Self)

    /** SH order of the coefficients; throws if their number belongs to no even order. */
    unsigned int ShOrder() const;
    unsigned int NumCoefficients() const;

    SH_CONVENTION GetShConvention() const;
    void SetShConvention(SH_CONVENTION ShConvention);

  protected:
    ShImage();
    ShImage(const ShImage& other);
    ~ShImage() override;

    Image::Pointer CreateRgbImage() const override;

    template<int num_components>
    Image::Pointer Construct() const;

    void  PrintSelf(std::ostream &os, itk::Indent indent) const override;

    SH_CONVENTION m_ShConvention;  // use mrtrix style SH convention
  };

} // namespace mitk

#endif
