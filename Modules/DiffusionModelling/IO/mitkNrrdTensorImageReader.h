/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkNrrdTensorImageReader_h
#define mitkNrrdTensorImageReader_h

#include "mitkCommon.h"
#include "itkVectorContainer.h"
#include "vnl/vnl_vector_fixed.h"
#include "mitkTensorImage.h"
#include "itkVectorImage.h"
#include "itkDiffusionTensor3D.h"
#include <mitkAbstractFileReader.h>
#include <mitkBaseData.h>
#include <mitkMimeType.h>

namespace mitk
{

  /** \brief
  */

  class NrrdTensorImageReader : public mitk::AbstractFileReader
  {
  public:

    typedef mitk::TensorImage OutputType;
    typedef itk::Matrix< double, 3, 3 > MeasurementFrameType;

    NrrdTensorImageReader(const NrrdTensorImageReader& other);
    NrrdTensorImageReader();
    ~NrrdTensorImageReader() override;

    using AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

  protected:


  private:
    NrrdTensorImageReader* Clone() const override;
  };

} //namespace MITK

#endif
