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

#include <iostream>
#include <fstream>


mitk::DiffusionImageNrrdWriter::DiffusionImageNrrdWriter()
  : AbstractFileWriter(mitk::Image::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE() ), mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_DESCRIPTION())
{
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

  VectorImageType::Pointer itkImg;
  mitk::CastToItkImage(input,itkImg);

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

  char keybuffer[512];
  char valbuffer[512];

  //itk::MetaDataDictionary dic = input->GetImage()->GetMetaDataDictionary();

  vnl_matrix_fixed<double,3,3> measurementFrame = mitk::DiffusionPropertyHelper::GetMeasurementFrame(input);
  if (measurementFrame(0,0) || measurementFrame(0,1) || measurementFrame(0,2) ||
    measurementFrame(1,0) || measurementFrame(1,1) || measurementFrame(1,2) ||
    measurementFrame(2,0) || measurementFrame(2,1) || measurementFrame(2,2))
  {
    sprintf( valbuffer, " (%lf,%lf,%lf) (%lf,%lf,%lf) (%lf,%lf,%lf)", measurementFrame(0,0), measurementFrame(0,1), measurementFrame(0,2), measurementFrame(1,0), measurementFrame(1,1), measurementFrame(1,2), measurementFrame(2,0), measurementFrame(2,1), measurementFrame(2,2));
    itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string("measurement frame"),std::string(valbuffer));
  }

  sprintf( valbuffer, "DWMRI");
  itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string("modality"),std::string(valbuffer));

  if (mitk::DiffusionPropertyHelper::GetOriginalGradientContainer(input)->Size())
  {
      MITK_INFO << "Saving original gradient directions";
      sprintf( valbuffer, "%1f", mitk::DiffusionPropertyHelper::GetReferenceBValue(input) );
      itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string("DWMRI_b-value"),std::string(valbuffer));

      for(unsigned int i=0; i<mitk::DiffusionPropertyHelper::GetGradientContainer(input)->Size(); i++)
      {
        sprintf( keybuffer, "DWMRI_gradient_%04d", i );

        sprintf( valbuffer, "%1f %1f %1f", mitk::DiffusionPropertyHelper::GetOriginalGradientContainer(input)->ElementAt(i).get(0),
          mitk::DiffusionPropertyHelper::GetOriginalGradientContainer(input)->ElementAt(i).get(1), mitk::DiffusionPropertyHelper::GetOriginalGradientContainer(input)->ElementAt(i).get(2));

        itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string(keybuffer),std::string(valbuffer));
      }
  }
  else if(mitk::DiffusionPropertyHelper::GetGradientContainer(input)->Size())
  {
      MITK_INFO << "Original gradient directions not found. Saving modified gradient directions";
      sprintf( valbuffer, "%1f", mitk::DiffusionPropertyHelper::GetReferenceBValue(input) );
      itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string("DWMRI_b-value"),std::string(valbuffer));

      for(unsigned int i=0; i<mitk::DiffusionPropertyHelper::GetGradientContainer(input)->Size(); i++)
      {
        sprintf( keybuffer, "DWMRI_gradient_%04d", i );

        sprintf( valbuffer, "%1f %1f %1f", mitk::DiffusionPropertyHelper::GetGradientContainer(input)->ElementAt(i).get(0),
          mitk::DiffusionPropertyHelper::GetGradientContainer(input)->ElementAt(i).get(1), mitk::DiffusionPropertyHelper::GetGradientContainer(input)->ElementAt(i).get(2));

        itk::EncapsulateMetaData<std::string>(itkImg->GetMetaDataDictionary(),std::string(keybuffer),std::string(valbuffer));
      }
  }

  typedef itk::VectorImage<short,3> ImageType;

  std::string ext = this->GetMimeType()->GetExtension(this->GetOutputLocation());
  ext = itksys::SystemTools::LowerCase(ext);

  // default extension is .nrrd
  if( ext == "")
  {
    ext = ".nrrd";
    this->SetOutputLocation(this->GetOutputLocation() + ext);
  }

  if (ext == ".hdwi" || ext == ".nrrd" || ext == ".dwi")
  {
    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->UseCompressionOn();

    typedef itk::ImageFileWriter<ImageType> WriterType;
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
