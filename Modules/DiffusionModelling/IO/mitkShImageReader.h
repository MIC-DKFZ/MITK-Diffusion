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

  /** \brief
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
