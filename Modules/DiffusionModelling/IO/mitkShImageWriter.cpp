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
#include <itkImageRegionIteratorWithIndex.h>

#include <cmath>

mitk::ShImageWriter::ShImageWriter()
  : AbstractFileWriter(mitk::ShImage::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionIOMimeTypes::SH_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::SH_MIMETYPE_DESCRIPTION())
{
  RegisterService();
}

mitk::ShImageWriter::ShImageWriter(const mitk::ShImageWriter& other)
  : AbstractFileWriter(other)
{
}

mitk::ShImageWriter::~ShImageWriter() = default;

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

  // The convention is marked in the file (NRRD key, NIfTI descrip), so that the coefficients read
  // back with it. MITK's own basis (MRtrix 0.2's) is written in MRtrix 3's orthonormal basis, which
  // MRtrix 3 and other tools expect: m != 0 divided by sqrt(2); m = 0 sits at index l(l+1)/2.
  std::string convention = "FSL";
  if (input->GetShConvention() == ShImage::MRTRIX)
  {
    convention = "MRtrix3";
    itk::ImageRegionIteratorWithIndex<ShImage::ShOnDiskType> it(image, image->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      const int j = static_cast<int>(it.GetIndex()[3]);
      bool zonal = false;
      for (int l = 0; l * (l + 1) / 2 <= j; l += 2)
        zonal = zonal || l * (l + 1) / 2 == j;
      if (!zonal)
        it.Set(static_cast<float>(it.Get() / std::sqrt(2.0)));
    }
  }
  if (ext==".shi" || ext==".nrrd")
    itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "MITK_SH_convention", convention);
  else
    itk::EncapsulateMetaData<std::string>(image->GetMetaDataDictionary(), "ITK_FileNotes", "MITK_SH_convention=" + convention);

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

