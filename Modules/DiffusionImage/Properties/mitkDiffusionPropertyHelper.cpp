/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionPropertyHelper.h"
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <itkImageRegionIterator.h>
#include <mitkProperties.h>
#include <vnl/algo/vnl_matrix_inverse.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>

namespace
{
  template <class TProperty>
  const TProperty* GetDiffusionProperty(const mitk::Image* image, const std::string& name)
  {
    // The property list owns the property, so the raw pointer outlives the temporary smart pointer.
    auto property = dynamic_cast<const TProperty*>(image->GetProperty(name.c_str()).GetPointer());
    if (property == nullptr)
      mitkThrow() << "Image has no property " << name << " of type " << TProperty::GetStaticNameOfClass();
    return property;
  }
}

const std::string mitk::DiffusionPropertyHelper::GRADIENTCONTAINERPROPERTYNAME = "DWMRI.GradientDirections";
const std::string mitk::DiffusionPropertyHelper::ORIGINALGRADIENTCONTAINERPROPERTYNAME = "DWMRI.OriginalGradientDirections";
const std::string mitk::DiffusionPropertyHelper::MEASUREMENTFRAMEPROPERTYNAME = "DWMRI.MeasurementFrame";
const std::string mitk::DiffusionPropertyHelper::REFERENCEBVALUEPROPERTYNAME = "DWMRI.ReferenceBValue";
const std::string mitk::DiffusionPropertyHelper::BVALUEMAPPROPERTYNAME = "DWMRI.BValueMap";
const std::string mitk::DiffusionPropertyHelper::MODALITY = "DWMRI.Modality";
const std::string mitk::DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS = "DWMRI.ApplyMatrixToGradients";
const std::string mitk::DiffusionPropertyHelper::APPLY_MF_TO_GRADIENTS = "DWMRI.ApplyMfToGradients";

mitk::DiffusionPropertyHelper::ImageType::Pointer mitk::DiffusionPropertyHelper::GetItkVectorImage(mitk::Image* image)
{
  ImageType::Pointer vectorImage = ImageType::New();
  mitk::CastToItkImage(image, vectorImage);
  return vectorImage;
}

mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer
mitk::DiffusionPropertyHelper::CalcAveragedDirectionSet(double precision, GradientDirectionsContainerType::ConstPointer directions)
{
  // save old and construct new direction container
  GradientDirectionsContainerType::Pointer newDirections = GradientDirectionsContainerType::New();

  // fill new direction container
  for(GradientDirectionsContainerType::ConstIterator gdcitOld = directions->Begin();
      gdcitOld != directions->End(); ++gdcitOld)
  {
    // already exists?
    bool found = false;
    for(GradientDirectionsContainerType::ConstIterator gdcitNew = newDirections->Begin();
        gdcitNew != newDirections->End(); ++gdcitNew)
    {
      if(AreAlike(gdcitNew.Value(), gdcitOld.Value(), precision))
      {
        found = true;
        break;
      }
    }

    // if not found, add it to new container
    if(!found)
    {
      newDirections->push_back(gdcitOld.Value());
    }
  }

  return newDirections;
}

void mitk::DiffusionPropertyHelper::AverageRedundantGradients(mitk::Image* image, double precision)
{
  GradientDirectionsContainerType::ConstPointer oldDirs = GetOriginalGradientContainer(image);
  GradientDirectionsContainerType::Pointer newDirs = CalcAveragedDirectionSet(precision, oldDirs);

  // if sizes equal, we do not need to do anything in this function
  if(oldDirs->size() == newDirs->size())
    return;

  // new image
  ImageType::Pointer oldImage = ImageType::New();
  mitk::CastToItkImage( image, oldImage);
  ImageType::Pointer newITKImage = ImageType::New();
  newITKImage->SetSpacing( oldImage->GetSpacing() );   // Set the image spacing
  newITKImage->SetOrigin( oldImage->GetOrigin() );     // Set the image origin
  newITKImage->SetDirection( oldImage->GetDirection() );  // Set the image direction
  newITKImage->SetLargestPossibleRegion( oldImage->GetLargestPossibleRegion() );
  newITKImage->SetVectorLength( newDirs->size() );
  newITKImage->SetBufferedRegion( oldImage->GetLargestPossibleRegion() );
  newITKImage->Allocate();

  // average image data that corresponds to identical directions
  itk::ImageRegionIterator< ImageType > newIt(newITKImage, newITKImage->GetLargestPossibleRegion());
  newIt.GoToBegin();
  itk::ImageRegionIterator< ImageType > oldIt(oldImage, oldImage->GetLargestPossibleRegion());
  oldIt.GoToBegin();

  // initial new value of voxel
  ImageType::PixelType newVec;
  newVec.SetSize(newDirs->size());

  // find which gradients should be averaged
  std::vector<std::vector<int> > dirIndices;
  for(GradientDirectionsContainerType::ConstIterator gdcitNew = newDirs->Begin();
      gdcitNew != newDirs->End(); ++gdcitNew)
  {
    dirIndices.push_back(std::vector<int>(0));
    for(GradientDirectionsContainerType::ConstIterator gdcitOld = oldDirs->Begin();
        gdcitOld != oldDirs->End(); ++gdcitOld)
    {
      if(AreAlike(gdcitNew.Value(), gdcitOld.Value(), precision))
      {
        //MITK_INFO << gdcitNew.Value() << "  " << gdcitOld.Value();
        dirIndices[gdcitNew.Index()].push_back(gdcitOld.Index());
      }
    }
  }

  //int ind1 = -1;
  while(!newIt.IsAtEnd())
  {
    // init new vector with zeros
    newVec.Fill(0.0);

    // the old voxel value with duplicates
    ImageType::PixelType oldVec = oldIt.Get();

    for(unsigned int i=0; i<dirIndices.size(); i++)
    {
      // do the averaging
      const int numavg = dirIndices[i].size();
      int sum = 0;
      for(int j=0; j<numavg; j++)
      {
        //MITK_INFO << newVec[i] << " << " << oldVec[dirIndices[i].at(j)];
        sum += oldVec[dirIndices[i].at(j)];
      }
      if(numavg == 0)
      {
        MITK_ERROR << "VectorImage: Error on averaging. Possibly due to corrupted data";
        return;
      }
      newVec[i] = sum / numavg;
    }

    newIt.Set(newVec);

    ++newIt;
    ++oldIt;
  }

  mitk::GrabItkImageMemory( newITKImage, image );

  SetOriginalGradientContainer(image, newDirs);
  InitializeImage(image);
}

void mitk::DiffusionPropertyHelper::GetGradientTransform(const mitk::Image* image, std::optional<MeasurementFrameType>& rotation, std::optional<MeasurementFrameType>& frame)
{
  bool apply_matrix = true;
  image->GetPropertyList()->GetBoolProperty(mitk::DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS.c_str(), apply_matrix);
  bool apply_mf = true;
  image->GetPropertyList()->GetBoolProperty(mitk::DiffusionPropertyHelper::APPLY_MF_TO_GRADIENTS.c_str(), apply_mf);

  rotation.reset();
  if (apply_matrix)
  {
    const mitk::Vector3D s = image->GetGeometry()->GetSpacing();
    MeasurementFrameType imageRotationMatrix;
    for (unsigned int c = 0; c < 3; ++c)
    {
      const mitk::VnlVector column = image->GetGeometry()->GetMatrixColumn(c)/s[c];
      for (unsigned int r = 0; r < 3; ++r)
        imageRotationMatrix[r][c] = column[r];
    }
    rotation = imageRotationMatrix;
  }

  frame.reset();
  if (apply_mf)
    frame = GetMeasurementFrame(image);
}

void mitk::DiffusionPropertyHelper::ApplyMeasurementFrameAndRotationMatrix(mitk::Image* image)
{
  if(  image->GetProperty(mitk::DiffusionPropertyHelper::ORIGINALGRADIENTCONTAINERPROPERTYNAME.c_str()).IsNull() )
    return;

  GradientDirectionsContainerType::ConstPointer originalDirections = GetOriginalGradientContainer(image);
  if( originalDirections.IsNull() || ( originalDirections->size() == 0 ) )
    return;

  std::optional<MeasurementFrameType> rotation, frame;
  GetGradientTransform(image, rotation, frame);
  MITK_INFO << "Modified gradients = original gradients" << (rotation ? " * image rotation" : "") << (frame ? " * measurement frame" : "");

  // Row vectors in the legacy order g * R * MF, one factor after the other as always, so that the
  // results stay the same to the last bit; the NRRD and NIfTI readers switch R off.
  GradientDirectionsContainerType::Pointer directions = GradientDirectionsContainerType::New();
  for (auto vec : originalDirections->CastToSTLConstContainer())
  {
    if (rotation)
      vec = vec * *rotation;
    if (frame)
      vec = vec * *frame;
    directions->push_back(vec);
  }

  SetGradientContainer(image, directions);
}

void mitk::DiffusionPropertyHelper::UnApplyMeasurementFrameAndRotationMatrix(mitk::Image* image)
{
  if(  image->GetProperty(mitk::DiffusionPropertyHelper::GRADIENTCONTAINERPROPERTYNAME.c_str()).IsNull() )
    return;

  GradientDirectionsContainerType::ConstPointer modifiedDirections = GetGradientContainer(image);
  if( modifiedDirections.IsNull() || ( modifiedDirections->size() == 0 ) )
    return;

  std::optional<MeasurementFrameType> rotation, frame;
  GetGradientTransform(image, rotation, frame);
  MITK_INFO << "Original gradients = modified gradients" << (frame ? " * inverse measurement frame" : "") << (rotation ? " * inverse image rotation" : "");
  if (rotation)
    rotation = MeasurementFrameType(vnl_matrix_inverse<double>(rotation->as_matrix()).pinverse());
  if (frame)
    frame = MeasurementFrameType(vnl_matrix_inverse<double>(frame->as_matrix()).pinverse());

  // The inverse of ApplyMeasurementFrameAndRotationMatrix, so the factors come in reverse order.
  GradientDirectionsContainerType::Pointer directions = GradientDirectionsContainerType::New();
  for (auto vec : modifiedDirections->CastToSTLConstContainer())
  {
    if (frame)
      vec = vec * *frame;
    if (rotation)
      vec = vec * *rotation;
    directions->push_back(vec);
  }

  SetOriginalGradientContainer(image, directions);
}

void mitk::DiffusionPropertyHelper::SetApplyMatrixToGradients(mitk::Image* image, bool apply)
{
  image->GetPropertyList()->SetProperty( mitk::DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS.c_str(), mitk::BoolProperty::New(apply) );
}

void mitk::DiffusionPropertyHelper::SetApplyMfToGradients(mitk::Image* image, bool apply)
{
  image->GetPropertyList()->SetProperty( mitk::DiffusionPropertyHelper::APPLY_MF_TO_GRADIENTS.c_str(), mitk::BoolProperty::New(apply) );
}

void mitk::DiffusionPropertyHelper::UpdateBValueMap(mitk::Image* image)
{
  BValueMapType b_ValueMap;

  if(image->GetProperty(GRADIENTCONTAINERPROPERTYNAME.c_str()).IsNotNull())
    b_ValueMap = mitk::BValueMapProperty::CreateBValueMap(GetGradientContainer(image), GetReferenceBValue(image));

  SetBValueMap(image, b_ValueMap);
}

bool mitk::DiffusionPropertyHelper::AreAlike(GradientDirectionType g1, GradientDirectionType g2, double precision)
{
  GradientDirectionType diff = g1 - g2;
  GradientDirectionType diff2 = g1 + g2;
  return diff.two_norm() < precision || diff2.two_norm() < precision;
}

void mitk::DiffusionPropertyHelper::CopyProperties(BaseData *source, BaseData *target, bool ignore_original_gradients)
{
  mitk::PropertyList::Pointer props = source->GetPropertyList()->Clone();
  if (ignore_original_gradients)
    props->RemoveProperty(mitk::DiffusionPropertyHelper::ORIGINALGRADIENTCONTAINERPROPERTYNAME);
  target->SetPropertyList(props);
}

void mitk::DiffusionPropertyHelper::CopyDICOMProperties(const mitk::BaseData* source, mitk::BaseData* target)
{
  mitk::PropertyList::Pointer list = source->GetPropertyList();
  for (auto key : list->GetPropertyKeys())
    if (key.starts_with("DICOM."))
      target->GetPropertyList()->ReplaceProperty(key, list->GetProperty(key));
}

void mitk::DiffusionPropertyHelper::InitializeImage(mitk::Image* image)
{
  if ( image->GetProperty(mitk::DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS.c_str()).IsNull() )
    image->SetProperty( mitk::DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS.c_str(), mitk::BoolProperty::New(true) );

  if ( image->GetProperty(mitk::DiffusionPropertyHelper::APPLY_MF_TO_GRADIENTS.c_str()).IsNull() )
    image->SetProperty( mitk::DiffusionPropertyHelper::APPLY_MF_TO_GRADIENTS.c_str(), mitk::BoolProperty::New(true) );

  if(  image->GetProperty(mitk::DiffusionPropertyHelper::ORIGINALGRADIENTCONTAINERPROPERTYNAME.c_str()).IsNull() )
  {
    // we don't have the original gradient directions. Therefore use the modified directions and rotate them back.
    UnApplyMeasurementFrameAndRotationMatrix(image);
  }
  else
    ApplyMeasurementFrameAndRotationMatrix(image);
  UpdateBValueMap(image);

  // initialize missing properties
  mitk::MeasurementFrameProperty::Pointer mf = dynamic_cast<mitk::MeasurementFrameProperty *>( image->GetProperty(MEASUREMENTFRAMEPROPERTYNAME.c_str()).GetPointer());
  if( mf.IsNull() )
  {
    //no measurement frame present, identity is assumed
    MeasurementFrameType identity;
    identity.set_identity();
    image->GetPropertyList()->ReplaceProperty( mitk::DiffusionPropertyHelper::MEASUREMENTFRAMEPROPERTYNAME.c_str(), mitk::MeasurementFrameProperty::New( identity ));
  }
}

void mitk::DiffusionPropertyHelper::InitializeImage(mitk::Image* image, GradientDirectionsContainerType::ConstPointer originalGradients, float referenceBValue, const MeasurementFrameType& measurementFrame, bool applyImageRotation)
{
  SetOriginalGradientContainer(image, originalGradients);
  SetMeasurementFrame(image, measurementFrame);
  SetReferenceBValue(image, referenceBValue);
  SetApplyMatrixToGradients(image, applyImageRotation);
  InitializeImage(image);
}

bool mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(const mitk::DataNode* node)
{
  if ( node==nullptr )
    return false;
  if ( node->GetData()==nullptr )
    return false;
  return IsDiffusionWeightedImage(dynamic_cast<mitk::Image *>(node->GetData()));
}


bool mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(const mitk::Image * image)
{
  bool isDiffusionWeightedImage( true );

  if( image == nullptr )
  {
    isDiffusionWeightedImage = false;
  }

  if( isDiffusionWeightedImage )
  {
    mitk::FloatProperty::Pointer referenceBValue = dynamic_cast<mitk::FloatProperty *>(image->GetProperty(REFERENCEBVALUEPROPERTYNAME.c_str()).GetPointer());

    if( referenceBValue.IsNull() )
    {
      isDiffusionWeightedImage = false;
    }

  }

  unsigned int gradientDirections( 0 );
  if( isDiffusionWeightedImage )
  {
    mitk::GradientDirectionsProperty::Pointer gradientDirectionsProperty = dynamic_cast<mitk::GradientDirectionsProperty *>(image->GetProperty(GRADIENTCONTAINERPROPERTYNAME.c_str()).GetPointer());

    if( gradientDirectionsProperty.IsNull() )
    {
      isDiffusionWeightedImage = false;
    }
    else
    {
      gradientDirections = gradientDirectionsProperty->GetGradientDirectionsContainer()->size();
    }
  }

  if( isDiffusionWeightedImage )
  {
    unsigned int components = image->GetPixelType().GetNumberOfComponents();

    if( components != gradientDirections )
    {
      isDiffusionWeightedImage = false;
    }
  }

  return isDiffusionWeightedImage;
}

const mitk::DiffusionPropertyHelper::BValueMapType & mitk::DiffusionPropertyHelper::GetBValueMap(const mitk::Image *image)
{
  return GetDiffusionProperty<mitk::BValueMapProperty>(image, BVALUEMAPPROPERTYNAME)->GetBValueMap();
}



float mitk::DiffusionPropertyHelper::GetReferenceBValue(const mitk::Image *image)
{
  return GetDiffusionProperty<mitk::FloatProperty>(image, REFERENCEBVALUEPROPERTYNAME)->GetValue();
}

mitk::DiffusionPropertyHelper::MeasurementFrameType mitk::DiffusionPropertyHelper::GetMeasurementFrame(const mitk::Image *image)
{
  auto mf = dynamic_cast<const mitk::MeasurementFrameProperty *>(image->GetProperty(MEASUREMENTFRAMEPROPERTYNAME.c_str()).GetPointer());
  if (mf == nullptr)
  {
    MeasurementFrameType identity;
    identity.set_identity();
    return identity;
  }
  return mf->GetMeasurementFrame();
}

const mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer mitk::DiffusionPropertyHelper::GetOriginalGradientContainer(const mitk::Image *image)
{
  return GetDiffusionProperty<mitk::GradientDirectionsProperty>(image, ORIGINALGRADIENTCONTAINERPROPERTYNAME)->GetGradientDirectionsContainer();
}

const mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer mitk::DiffusionPropertyHelper::GetGradientContainer(const mitk::Image *image)
{
  return GetDiffusionProperty<mitk::GradientDirectionsProperty>(image, GRADIENTCONTAINERPROPERTYNAME)->GetGradientDirectionsContainer();
}

void mitk::DiffusionPropertyHelper::SetReferenceBValue(mitk::Image* image, float b_value)
{
  image->GetPropertyList()->ReplaceProperty(REFERENCEBVALUEPROPERTYNAME.c_str(), mitk::FloatProperty::New(b_value));
}

void mitk::DiffusionPropertyHelper::SetBValueMap(mitk::Image* image, BValueMapType map)
{
  image->GetPropertyList()->ReplaceProperty( mitk::DiffusionPropertyHelper::BVALUEMAPPROPERTYNAME.c_str(), mitk::BValueMapProperty::New(map));
}

void mitk::DiffusionPropertyHelper::SetOriginalGradientContainer(mitk::Image* image, GradientDirectionsContainerType::ConstPointer g_cont)
{
  image->GetPropertyList()->ReplaceProperty(ORIGINALGRADIENTCONTAINERPROPERTYNAME.c_str(), mitk::GradientDirectionsProperty::New(g_cont));
}

void mitk::DiffusionPropertyHelper::SetOriginalGradientContainer(mitk::Image* image, GradientDirectionsContainerType::Pointer g_cont)
{
  image->GetPropertyList()->ReplaceProperty(ORIGINALGRADIENTCONTAINERPROPERTYNAME.c_str(), mitk::GradientDirectionsProperty::New(g_cont));
}

void mitk::DiffusionPropertyHelper::SetGradientContainer(mitk::Image* image, GradientDirectionsContainerType::ConstPointer g_cont)
{
  image->GetPropertyList()->ReplaceProperty(GRADIENTCONTAINERPROPERTYNAME.c_str(), mitk::GradientDirectionsProperty::New(g_cont));
}

void mitk::DiffusionPropertyHelper::SetGradientContainer(mitk::Image* image, GradientDirectionsContainerType::Pointer g_cont)
{
  image->GetPropertyList()->ReplaceProperty(GRADIENTCONTAINERPROPERTYNAME.c_str(), mitk::GradientDirectionsProperty::New(g_cont));
}

void mitk::DiffusionPropertyHelper::SetMeasurementFrame(mitk::Image* image, MeasurementFrameType mf)
{
  image->GetPropertyList()->ReplaceProperty( MEASUREMENTFRAMEPROPERTYNAME.c_str(), mitk::MeasurementFrameProperty::New( mf ) );
}

void mitk::DiffusionPropertyHelper::SetupProperties()
{
  // No property persistence infos: the plain image writers would store these properties through
  // GetValueAsString(), which is empty for gradients and measurement frame. The DWI writers write
  // them explicitly.
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::BVALUEMAPPROPERTYNAME, "This map stores which b values belong to which gradients.");
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::ORIGINALGRADIENTCONTAINERPROPERTYNAME, "The original gradients used during acquisition. This property may be empty.");
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::REFERENCEBVALUEPROPERTYNAME, "The reference b value the gradients are normalized to.");
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::MEASUREMENTFRAMEPROPERTYNAME, "The measurment frame used during acquisition.");
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::GRADIENTCONTAINERPROPERTYNAME, "The gradients after applying measurement frame and image matrix.");
  mitk::CoreServices::GetPropertyDescriptions()->AddDescription(mitk::DiffusionPropertyHelper::MODALITY, "Defines the modality used for acquisition. DWMRI signifies diffusion weighted images.");
}
