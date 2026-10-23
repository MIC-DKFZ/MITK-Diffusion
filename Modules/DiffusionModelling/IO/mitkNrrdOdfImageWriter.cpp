/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkNrrdOdfImageWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkImageFileWriter.h"
#include "mitkImageToItk.h"
#include "mitkIOMimeTypes.h"
#include "mitkDiffusionIOMimeTypes.h"
#include <mitkLocaleSwitch.h>


mitk::NrrdOdfImageWriter::NrrdOdfImageWriter()
  : AbstractFileWriter(mitk::OdfImage::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionIOMimeTypes::ODF_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::ODF_MIMETYPE_DESCRIPTION())
{
  RegisterService();
}

mitk::NrrdOdfImageWriter::NrrdOdfImageWriter(const mitk::NrrdOdfImageWriter& other)
  : AbstractFileWriter(other)
{
}

mitk::NrrdOdfImageWriter::~NrrdOdfImageWriter() = default;


void mitk::NrrdOdfImageWriter::Write()
{
  InputType::ConstPointer input = dynamic_cast<const InputType*>(this->GetInput());
    if (input.IsNull())
        mitkThrow() << "Sorry, input to NrrdOdfImageWriter is nullptr!";
    if ( this->GetOutputLocation().empty() )
        mitkThrow() << "Sorry, filename has not been set!";

    mitk::LocaleSwitch localeSwitch("C");

    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->UseCompressionOn();

    // The file layout is a VectorImage; it is written straight from the ODF buffer, without a copy.
    typedef itk::VectorImage<float, 3> VecImgType;
    typedef itk::ImageFileWriter<VecImgType> WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();

    typedef mitk::ImageToItk<InputType::ItkOdfImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(input);
    caster->Update();
    InputType::ItkOdfImageType::Pointer outimage = caster->GetOutput();

    VecImgType::Pointer vecImg = VecImgType::New();
    vecImg->SetSpacing( outimage->GetSpacing() );   // Set the image spacing
    vecImg->SetOrigin( outimage->GetOrigin() );     // Set the image origin
    vecImg->SetDirection( outimage->GetDirection() );  // Set the image direction
    vecImg->SetLargestPossibleRegion( outimage->GetLargestPossibleRegion());
    vecImg->SetBufferedRegion( outimage->GetLargestPossibleRegion() );
    vecImg->SetVectorLength(ODF_SAMPLING_SIZE);
    vecImg->GetPixelContainer()->SetImportPointer(outimage->GetBufferPointer()->GetDataPointer(),
                                                  outimage->GetLargestPossibleRegion().GetNumberOfPixels() * ODF_SAMPLING_SIZE, false);

    nrrdWriter->SetInput( vecImg );
    nrrdWriter->SetImageIO(io);
    nrrdWriter->SetFileName(this->GetOutputLocation().c_str());
    nrrdWriter->UseCompressionOn();


    nrrdWriter->Update();
}

mitk::NrrdOdfImageWriter* mitk::NrrdOdfImageWriter::Clone() const
{
  return new NrrdOdfImageWriter(*this);
}

