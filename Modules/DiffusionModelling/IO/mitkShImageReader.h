/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkShImageReader_h
#define mitkShImageReader_h

#include "mitkCommon.h"
#include "itkVectorContainer.h"
#include "vnl/vnl_vector_fixed.h"
#include "mitkShImage.h"
#include "itkVectorImage.h"
#include <mitkAbstractFileReader.h>
#include <mitkBaseData.h>
#include <mitkMimeType.h>

namespace mitk
{

  /**
   * Reads SH coefficient images in MRtrix's layout: a 4D float image (NIfTI, or NRRD as .shi) with
   * one volume per coefficient, orders 2 to 12. The image holds coefficients in MITK's basis
   * (mitk::sh, mrtrix = true: MRtrix 0.2's) or, for files marked FSL, in MITK's FSL definition.
   * A convention marked by ShImageWriter ("MITK_SH_convention" NRRD key or NIfTI descrip) decides;
   * otherwise the option "SH convention": "auto" reads files whose descrip says
   * "MRtrix version: 3..." (every MRtrix FOD since 3.0) in MRtrix 3's orthonormal basis and all
   * others, which includes files of earlier MITK versions, in MITK's legacy basis. MRtrix 3
   * coefficients are converted on reading (m != 0 times sqrt(2)). "FSL/Dipy" reads coefficients
   * of dipy's default basis (descoteaux07, legacy), given in the frame of FSL bvecs, and converts
   * them into MITK's basis at world directions.
   *
   * SH images in vector layout (a 3D image with one vector of coefficients per voxel, as MITK's
   * generic image writer and older MITK versions wrote them) are not read as SH images; MITK's
   * image reader reads them as plain multi-component images. The two such files in mitk-diffusion-test-data
   * (qBallImage_SHCoeffs.nrrd, StreamlineTractography/qball_image_sh.nii.gz) hold coefficients in
   * the FSL/Dipy basis, not MITK's, so reading them with MITK's default would give wrong ODFs.
   */
  class ShImageReader : public mitk::AbstractFileReader
  {
  public:

    typedef mitk::ShImage OutputType;
    typedef itk::VectorImage<float,3> VectorImageType;

    ShImageReader(const ShImageReader& other);
    ShImageReader();
    ~ShImageReader() override;

    using AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

  protected:


  private:
    ShImageReader* Clone() const override;

    template <int sh_order>
    mitk::Image::Pointer ConvertShImage(ShImage::ShOnDiskType::Pointer img);
  };

} //namespace MITK

#endif
