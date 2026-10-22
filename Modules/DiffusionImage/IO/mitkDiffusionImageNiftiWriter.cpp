/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageNiftiWriter__cpp
#define __mitkDiffusionImageNiftiWriter__cpp

#include "mitkDiffusionImageNiftiWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkNiftiImageIO.h"
#include "itkImageFileWriter.h"
#include "itksys/SystemTools.hxx"
#include "mitkDiffusionImageMimeTypes.h"
#include "mitkImageCast.h"
#include <mitkLocaleSwitch.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <vnl/vnl_inverse.h>
#include <iostream>
#include <fstream>

mitk::DiffusionImageNiftiWriter::DiffusionImageNiftiWriter()
  : AbstractFileWriter(mitk::Image::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_NAME() ), mitk::DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_DESCRIPTION())
{
  // Above MITK's plain image writers, so that DWIs are saved with their gradients by default;
  // GetConfidenceLevel() rejects images that are not DWIs.
  this->SetRanking(10);
  RegisterService();
}

mitk::DiffusionImageNiftiWriter::DiffusionImageNiftiWriter(const mitk::DiffusionImageNiftiWriter& other)
  : AbstractFileWriter(other)
{
}

mitk::DiffusionImageNiftiWriter::~DiffusionImageNiftiWriter() = default;

void mitk::DiffusionImageNiftiWriter::Write()
{
  mitk::Image::ConstPointer input = dynamic_cast<const mitk::Image *>(this->GetInput());
  if (input.IsNull())
  {
    MITK_ERROR <<"Sorry, input to DiffusionImageNiftiWriter is nullptr!";
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

  const auto mimeType = this->GetRegisteredMimeType();

  // default extension is .nii.gz
  if (!mimeType.MatchesExtension(this->GetOutputLocation()))
    this->SetOutputLocation(this->GetOutputLocation() + ".nii.gz");

  MITK_INFO << "Writing Nifti-Image";

  using ImageType4D = itk::Image<short,4>;
  ImageType4D::Pointer img4 = mitk::imv::VectorImageToVolumes(itkImg.GetPointer());

  itk::NiftiImageIO::Pointer io4 = itk::NiftiImageIO::New();

  using WriterType4 = itk::ImageFileWriter<ImageType4D>;
  WriterType4::Pointer nrrdWriter4 = WriterType4::New();
  nrrdWriter4->UseInputMetaDataDictionaryOn();
  nrrdWriter4->SetInput( img4 );
  nrrdWriter4->SetFileName(this->GetOutputLocation());
  nrrdWriter4->UseCompressionOn();
  nrrdWriter4->SetImageIO(io4);
  try
  {
    nrrdWriter4->Update();
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.GetDescription() << std::endl;
    throw;
  }


  std::string base_path = itksys::SystemTools::GetFilenamePath(this->GetOutputLocation());
  std::string bvals_file = mimeType.GetFilenameWithoutExtension(this->GetOutputLocation());
  if (!base_path.empty())
      bvals_file = base_path + "/" + bvals_file;
  bvals_file += ".bvals";

  std::string bvecs_file = mimeType.GetFilenameWithoutExtension(this->GetOutputLocation());
  if (!base_path.empty())
      bvecs_file = base_path + "/" + bvecs_file;
  bvecs_file += ".bvecs";

  // FSL bvecs: the modified gradients are world vectors, turned into the voxel frame the reader
  // expects.
  const auto worldFromBvecs = mitk::gradients::FslBvecsMeasurementFrame(itkImg->GetDirection()).transpose();
  const auto bvecsFromWorld = vnl_inverse(worldFromBvecs);
  auto bvecs = GradientDirectionContainerType::New();
  for (const auto& g : mitk::DiffusionPropertyHelper::GetGradientContainer(input)->CastToSTLConstContainer())
    bvecs->push_back(bvecsFromWorld * g);
  mitk::gradients::WriteBvalsBvecs(bvals_file, bvecs_file, bvecs.GetPointer(), mitk::DiffusionPropertyHelper::GetReferenceBValue(input));
}

mitk::DiffusionImageNiftiWriter* mitk::DiffusionImageNiftiWriter::Clone() const
{
  return new DiffusionImageNiftiWriter(*this);
}

mitk::IFileWriter::ConfidenceLevel mitk::DiffusionImageNiftiWriter::GetConfidenceLevel() const
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

#endif //__mitkDiffusionImageNiftiWriter__cpp
