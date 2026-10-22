/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageNiftiReader_h
#define __mitkDiffusionImageNiftiReader_h

#include "mitkCommon.h"

// MITK includes
#include "mitkImageSource.h"

#include <mitkDiffusionPropertyHelper.h>

// ITK includes
#include "itkVectorImage.h"
#include "mitkAbstractFileReader.h"



namespace mitk
{

  /** \brief
  */

  class DiffusionImageNiftiReader : public mitk::AbstractFileReader
  {
  public:

    DiffusionImageNiftiReader(const DiffusionImageNiftiReader & other);
    DiffusionImageNiftiReader();
    ~DiffusionImageNiftiReader() override;

    using AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

    using DiffusionPixelType = short;

    using OutputType = mitk::Image;
    using VectorImageType = mitk::DiffusionPropertyHelper::ImageType;
    using GradientDirectionType = mitk::DiffusionPropertyHelper::GradientDirectionType;
    using MeasurementFrameType = mitk::DiffusionPropertyHelper::MeasurementFrameType;
    using GradientDirectionContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

  protected:
    OutputType::Pointer InternalRead();

  private:

    DiffusionImageNiftiReader* Clone() const override;
    us::ServiceRegistration<mitk::IFileReader> m_ServiceReg;
  };

} //namespace MITK

#endif // __mitkDiffusionImageNiftiReader_h
