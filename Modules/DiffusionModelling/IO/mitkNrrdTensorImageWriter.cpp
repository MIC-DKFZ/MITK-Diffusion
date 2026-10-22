/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkNrrdTensorImageWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkImageFileWriter.h"
#include "itkDiffusionTensor3D.h"
#include "mitkImageCast.h"
#include "mitkDiffusionIOMimeTypes.h"
#include <mitkLocaleSwitch.h>


mitk::NrrdTensorImageWriter::NrrdTensorImageWriter()
  : AbstractFileWriter(mitk::TensorImage::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionIOMimeTypes::DTI_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::DTI_MIMETYPE_DESCRIPTION() )
{
  RegisterService();
}

mitk::NrrdTensorImageWriter::NrrdTensorImageWriter(const mitk::NrrdTensorImageWriter& other)
  : AbstractFileWriter(other)
{
}


mitk::NrrdTensorImageWriter::~NrrdTensorImageWriter()
{}


void mitk::NrrdTensorImageWriter::Write()
{
  InputType::ConstPointer input = dynamic_cast<const InputType*>(this->GetInput());
  if (input.IsNull() )
    mitkThrow() << "Sorry, input to NrrdTensorImageWriter is nullptr!";
  if ( this->GetOutputLocation().empty() )
    mitkThrow() << "Sorry, filename has not been set!";
  mitk::LocaleSwitch localeSwitch("C");

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  io->UseCompressionOn();

  typedef itk::ImageFileWriter<TensorImage::ItkTensorImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();

  TensorImage::ItkTensorImageType::Pointer outimage = TensorImage::ItkTensorImageType::New();
  CastToItkImage(input, outimage);

  nrrdWriter->SetInput( outimage );
  nrrdWriter->SetImageIO(io);
  nrrdWriter->SetFileName(this->GetOutputLocation().c_str());
  nrrdWriter->UseCompressionOn();

  nrrdWriter->Update();
}

mitk::NrrdTensorImageWriter* mitk::NrrdTensorImageWriter::Clone() const
{
  return new NrrdTensorImageWriter(*this);
}

mitk::IFileWriter::ConfidenceLevel mitk::NrrdTensorImageWriter::GetConfidenceLevel() const
{
  InputType::ConstPointer input = dynamic_cast<const InputType*>(this->GetInput());
  if (input.IsNull() )
  {
    return Unsupported;
  }
  else
  {
    return Supported;
  }
}
