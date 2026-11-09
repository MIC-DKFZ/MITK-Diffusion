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
#include <itkImage.h>
#include <itkMatrix.h>
#include <itkMetaDataDictionary.h>
#include <itkVector.h>
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
   * generic image writer and older MITK versions wrote them, e.g. qBallImage_SHCoeffs.nrrd and
   * StreamlineTractography/qball_image_sh.nii.gz) are read only when this reader is chosen
   * explicitly (open dialog, PreferenceListReaderOptionsFunctor): its confidence for them is
   * PartiallySupported, so MITK's image reader stays the default and reads them as plain
   * multi-component images. "auto" reads them in MITK's FSL definition (ShImage::FSL), the basis
   * of the Q-ball filters of those versions.
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
    ConfidenceLevel GetConfidenceLevel() const override;

  protected:


  private:
    ShImageReader* Clone() const override;

    template <int sh_order>
    using CoefficientImageType = itk::Image<itk::Vector<float, (sh_order + 1) * (sh_order + 2) / 2>, 3>;

    template <int sh_order>
    mitk::Image::Pointer ConvertShImage(ShImage::ShOnDiskType::Pointer img);

    template <int sh_order>
    mitk::Image::Pointer ConvertVectorShImage(VectorImageType::Pointer img);

    /** Converts the coefficients into MITK's basis as the file and the option say; vectorLayout: read from a 3D vector image. */
    template <int sh_order>
    mitk::Image::Pointer TakeCoefficients(CoefficientImageType<sh_order> *coefficients,
                                          const itk::MetaDataDictionary &dictionary,
                                          const itk::Matrix<double, 3, 3> &direction,
                                          bool vectorLayout);

    mitk::Image::Pointer ReadVectorLayout(const std::string &location);
  };

} //namespace MITK

#endif
