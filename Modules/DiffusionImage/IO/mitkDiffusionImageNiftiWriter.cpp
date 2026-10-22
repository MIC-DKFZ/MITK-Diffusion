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

  using ImageType = itk::VectorImage<short,3>;

  const auto mimeType = this->GetRegisteredMimeType();

  // default extension is .nii.gz
  if (!mimeType.MatchesExtension(this->GetOutputLocation()))
    this->SetOutputLocation(this->GetOutputLocation() + ".nii.gz");

  MITK_INFO << "Writing Nifti-Image";

  using ImageType4D = itk::Image<short,4>;
  ImageType4D::Pointer img4 = ImageType4D::New();

  ImageType::SpacingType spacing = itkImg->GetSpacing();
  ImageType4D::SpacingType spacing4;
  for(int i=0; i<3; i++)
    spacing4[i] = spacing[i];
  spacing4[3] = 1;
  img4->SetSpacing( spacing4 );   // Set the image spacing

  ImageType::PointType origin = itkImg->GetOrigin();
  ImageType4D::PointType origin4;
  for(int i=0; i<3; i++)
    origin4[i] = origin[i];
  origin4[3] = 0;
  img4->SetOrigin( origin4 );     // Set the image origin

  ImageType::DirectionType direction = itkImg->GetDirection();
  ImageType4D::DirectionType direction4;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      direction4[i][j] = direction[i][j];
  for(int i=0; i<4; i++)
    direction4[i][3] = 0;
  for(int i=0; i<4; i++)
    direction4[3][i] = 0;
  direction4[3][3] = 1;
  img4->SetDirection( direction4 );  // Set the image direction

  ImageType::RegionType region = itkImg->GetLargestPossibleRegion();
  ImageType4D::RegionType region4;

  ImageType::RegionType::SizeType size = region.GetSize();
  ImageType4D::RegionType::SizeType size4;

  for(int i=0; i<3; i++)
    size4[i] = size[i];
  size4[3] = itkImg->GetVectorLength();

  ImageType::RegionType::IndexType index = region.GetIndex();
  ImageType4D::RegionType::IndexType index4;
  for(int i=0; i<3; i++)
    index4[i] = index[i];
  index4[3] = 0;

  region4.SetSize(size4);
  region4.SetIndex(index4);
  img4->SetRegions( region4 );

  img4->Allocate();

  itk::ImageRegionIterator<ImageType>   it (itkImg, itkImg->GetLargestPossibleRegion() );
  using VecPixType = ImageType::PixelType;

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    VecPixType vec = it.Get();
    ImageType::IndexType currentIndex = it.GetIndex();
    for(unsigned int ind=0; ind<vec.Size(); ind++)
    {

      for(int i=0; i<3; i++)
        index4[i] = currentIndex[i];
      index4[3] = ind;
      img4->SetPixel(index4, vec[ind]);
    }
  }

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
