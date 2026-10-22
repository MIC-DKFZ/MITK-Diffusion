/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkShImageWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkImageFileWriter.h"
#include "mitkImageCast.h"
#include "mitkIOMimeTypes.h"
#include "mitkDiffusionIOMimeTypes.h"
#include <mitkLocaleSwitch.h>
#include <itkShCoefficientImageExporter.h>
#include <itksys/SystemTools.hxx>
#include <itkNrrdImageIO.h>

mitk::ShImageWriter::ShImageWriter()
  : AbstractFileWriter(mitk::ShImage::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionIOMimeTypes::SH_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::SH_MIMETYPE_DESCRIPTION())
{
  RegisterService();
}

mitk::ShImageWriter::ShImageWriter(const mitk::ShImageWriter& other)
  : AbstractFileWriter(other)
{
}

mitk::ShImageWriter::~ShImageWriter()
{}

template <int shOrder>
void mitk::ShImageWriter::WriteShImage(InputType::ConstPointer input)
{
  mitk::LocaleSwitch localeSwitch("C");
  typename itk::ShCoefficientImageExporter< float, shOrder >::InputImageType::Pointer itk_image = itk::ShCoefficientImageExporter< float, shOrder >::InputImageType::New();
  CastToItkImage(input, itk_image);

  typedef itk::ShCoefficientImageExporter< float, shOrder > ExporterType;
  typename ExporterType::Pointer exporter = ExporterType::New();
  exporter->SetInputImage(itk_image);
  exporter->GenerateData();

  std::string ext = itksys::SystemTools::GetFilenameLastExtension(this->GetOutputLocation());

  typedef itk::ImageFileWriter<ShImage::ShOnDiskType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  ShImage::ShOnDiskType::Pointer image = exporter->GetOutputImage();

  // FSL/Dipy coefficients are marked in the file (NRRD key, NIfTI descrip), so that they read back
  // with their convention; files in MITK's own MRtrix convention stay as before.
  if (input->GetShConvention() == ShImage::FSL)
  {
    if (ext==".shi" || ext==".nrrd")
      itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "MITK_SH_convention", "FSL");
    else
      itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "ITK_FileNotes", "MITK_SH_convention=FSL");
  }

  writer->SetInput( image );

  if (ext==".shi")
  {
    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->UseCompressionOn();
    writer->SetImageIO(io);
  }

  writer->SetFileName(this->GetOutputLocation().c_str());
  writer->SetUseCompression(true);

  writer->Update();
}

void mitk::ShImageWriter::Write()
{
  InputType::ConstPointer input = dynamic_cast<const InputType*>(this->GetInput());
  if (input.IsNull())
    mitkThrow() << "Sorry, input to ShImageWriter is nullptr!";
  if ( this->GetOutputLocation().empty() )
    mitkThrow() << "Sorry, filename has not been set!";

  switch (input->GetImageDescriptor()->GetChannelTypeById(0).GetNumberOfComponents())
  {
  case 6:
    WriteShImage<2>(input);
    break;
  case 15:
    WriteShImage<4>(input);
    break;
  case 28:
    WriteShImage<6>(input);
    break;
  case 45:
    WriteShImage<8>(input);
    break;
  case 66:
    WriteShImage<10>(input);
    break;
  case 91:
    WriteShImage<12>(input);
    break;
  default :
    mitkThrow() << "SH order larger 12 not supported";
  }
}

mitk::ShImageWriter* mitk::ShImageWriter::Clone() const
{
  return new ShImageWriter(*this);
}

