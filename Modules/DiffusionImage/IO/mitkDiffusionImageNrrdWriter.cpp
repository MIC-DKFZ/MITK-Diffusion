/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageNrrdWriter__cpp
#define __mitkDiffusionImageNrrdWriter__cpp

#include "mitkDiffusionImageNrrdWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkImageFileWriter.h"
#include "itksys/SystemTools.hxx"
#include "mitkDiffusionImageMimeTypes.h"
#include "mitkImageCast.h"
#include <mitkLocaleSwitch.h>

#include <format>

#include <iostream>
#include <fstream>


mitk::DiffusionImageNrrdWriter::DiffusionImageNrrdWriter()
  : AbstractFileWriter(mitk::Image::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_NAME() ), mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_DESCRIPTION())
{
  // Above MITK's plain image writers, so that DWIs are saved with their gradients by default;
  // GetConfidenceLevel() rejects images that are not DWIs.
  this->SetRanking(10);
  RegisterService();
}

mitk::DiffusionImageNrrdWriter::DiffusionImageNrrdWriter(const mitk::DiffusionImageNrrdWriter& other)
  : AbstractFileWriter(other)
{
}

mitk::DiffusionImageNrrdWriter::~DiffusionImageNrrdWriter()
{}

void mitk::DiffusionImageNrrdWriter::Write()
{
  mitk::Image::ConstPointer input = dynamic_cast<const mitk::Image *>(this->GetInput());
  if (input.IsNull())
  {
    MITK_ERROR <<"Sorry, input to DiffusionImageNrrdWriter is nullptr!";
    return;
  }
  if ( this->GetOutputLocation().empty() )
  {
    MITK_ERROR << "Sorry, filename has not been set!";
    return ;
  }
  mitk::LocaleSwitch localeSwitch("C");

  VectorImageType::Pointer itkImg;
  mitk::CastToItkImage(input,itkImg);
  auto& dictionary = itkImg->GetMetaDataDictionary();

  // The modified gradients are world vectors. They are written with an identity measurement frame,
  // as the standard NRRD field and as MITK's own key (which earlier MITK versions read); together
  // they mark the file as written by this version. std::format("{}") writes the shortest text that
  // reads back to the same double.
  const std::vector<std::vector<double>> identity = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  itk::EncapsulateMetaData<std::vector<std::vector<double>>>(dictionary, "NRRD_measurement frame", identity);
  itk::EncapsulateMetaData<std::string>(dictionary, "measurement frame", " (1,0,0) (0,1,0) (0,0,1)");
  itk::EncapsulateMetaData<std::string>(dictionary, "modality", "DWMRI");
  itk::EncapsulateMetaData<std::string>(dictionary, "DWMRI_b-value", std::format("{}", mitk::DiffusionPropertyHelper::GetReferenceBValue(input)));
  const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(input);
  for (unsigned int i = 0; i < gradients->Size(); ++i)
  {
    const auto& g = gradients->ElementAt(i);
    itk::EncapsulateMetaData<std::string>(dictionary, std::format("DWMRI_gradient_{:04}", i), std::format("{} {} {}", g[0], g[1], g[2]));
  }

  // default extension is .nrrd
  if (!this->GetRegisteredMimeType().MatchesExtension(this->GetOutputLocation()))
    this->SetOutputLocation(this->GetOutputLocation() + ".nrrd");

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  io->UseCompressionOn();

  typedef itk::ImageFileWriter<VectorImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( itkImg );
  nrrdWriter->SetImageIO(io);
  nrrdWriter->SetFileName(this->GetOutputLocation());
  nrrdWriter->UseCompressionOn();
  nrrdWriter->SetImageIO(io);
  try
  {
    nrrdWriter->Update();
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.GetDescription() << std::endl;
    throw;
  }
}

mitk::DiffusionImageNrrdWriter* mitk::DiffusionImageNrrdWriter::Clone() const
{
  return new DiffusionImageNrrdWriter(*this);
}

mitk::IFileWriter::ConfidenceLevel mitk::DiffusionImageNrrdWriter::GetConfidenceLevel() const
{
  mitk::Image::ConstPointer input = dynamic_cast<const mitk::Image*>(this->GetInput());
  if (input.IsNull() || !mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage( input ) )
  {
    return Unsupported;
  }
  else
  {
    return Supported;
  }
}

#endif //__mitkDiffusionImageNrrdWriter__cpp
