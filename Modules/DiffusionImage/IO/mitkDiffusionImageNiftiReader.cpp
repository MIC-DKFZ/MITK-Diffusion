/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageNiftiReader_cpp
#define __mitkDiffusionImageNiftiReader_cpp

#include "mitkDiffusionImageNiftiReader.h"

#include <iostream>
#include <fstream>


// Diffusion properties
#include <mitkGradientDirectionsProperty.h>
#include <mitkBValueMapProperty.h>
#include <mitkMeasurementFrameProperty.h>
#include <mitkProperties.h>

// ITK includes
#include <itkImageRegionIterator.h>
#include <itkMetaDataObject.h>
#include "itksys/SystemTools.hxx"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include "itkNiftiImageIO.h"
#include <itkBruker2dseqImageIO.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "mitkCustomMimeType.h"
#include "mitkDiffusionImageMimeTypes.h"

#include <mitkITKImageImport.h>
#include <mitkImageWriteAccessor.h>
#include <mitkImageDataItem.h>
#include <mitkLocaleSwitch.h>
#include "mitkIOUtil.h"
#include <mitkDiffusionImageHelperFunctions.h>

namespace mitk
{

DiffusionImageNiftiReader::
DiffusionImageNiftiReader(const DiffusionImageNiftiReader & other)
  : AbstractFileReader(other)
{
}


DiffusionImageNiftiReader* DiffusionImageNiftiReader::Clone() const
{
  return new DiffusionImageNiftiReader(*this);
}


DiffusionImageNiftiReader::
~DiffusionImageNiftiReader()
{}

  DiffusionImageNiftiReader::
DiffusionImageNiftiReader()
  : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_NAME() ), mitk::DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_DESCRIPTION() )

{
  m_ServiceReg = this->RegisterService();
}

std::vector<itk::SmartPointer<mitk::BaseData> >
DiffusionImageNiftiReader::
DoRead()
{
  std::vector<itk::SmartPointer<mitk::BaseData> > result;

  // Since everything is completely read in GenerateOutputInformation() it is stored
  // in a cache variable. A timestamp is associated.
  // If the timestamp of the cache variable is newer than the MTime, we only need to
  // assign the cache variable to the DataObject.
  // Otherwise, the tree must be read again from the file and OuputInformation must
  // be updated!

  if(m_OutputCache.IsNull()) InternalRead();

  result.push_back(m_OutputCache.GetPointer());
  return result;
}


void DiffusionImageNiftiReader::InternalRead()
{
  OutputType::Pointer outputForCache = OutputType::New();
  if ( this->GetInputLocation() == "")
  {
    throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename to be read is empty!");
  }
  else
  {
    try
    {
      mitk::LocaleSwitch localeSwitch("C");

      MITK_INFO << "DiffusionImageNiftiReader: reading image information";
      VectorImageType::Pointer itkVectorImage;

      const auto mimeType = this->GetRegisteredMimeType();
      const bool isNifti = mimeType.MatchesExtension(this->GetInputLocation());

      itk::Bruker2dseqImageIO::Pointer bruker_io = itk::Bruker2dseqImageIO::New();

      typedef itk::Image<DiffusionPixelType,4> ImageType4D;
      ImageType4D::Pointer img4 = ImageType4D::New();
      bool bruker = false;

      if (isNifti)
      {
        // create reader and read file
        itk::NiftiImageIO::Pointer io2 = itk::NiftiImageIO::New();
        if (!io2->CanReadFile(this->GetInputLocation().c_str()))
          mitkThrow() << "ITK cannot read the NIfTI file " << this->GetInputLocation();
        typedef itk::ImageFileReader<ImageType4D> FileReaderType;
        FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetFileName( this->GetInputLocation() );
        reader->SetImageIO(io2);
        reader->Update();
        img4 = reader->GetOutput();
      }
      else if(bruker_io->CanReadFile(this->GetInputLocation().c_str()))
      {
        MITK_INFO << "Reading bruker 2dseq file";
        typedef itk::ImageFileReader<ImageType4D> FileReaderType;
        FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetFileName( this->GetInputLocation() );
        reader->SetImageIO(bruker_io);
        reader->Update();
        img4 = reader->GetOutput();
        bruker = true;
      }

      // convert 4D file to vector image
      itkVectorImage = VectorImageType::New();

      VectorImageType::SpacingType spacing;
      ImageType4D::SpacingType spacing4 = img4->GetSpacing();
      for(int i=0; i<3; i++)
        spacing[i] = spacing4[i];
      itkVectorImage->SetSpacing( spacing );   // Set the image spacing

      VectorImageType::PointType origin;
      ImageType4D::PointType origin4 = img4->GetOrigin();
      for(int i=0; i<3; i++)
        origin[i] = origin4[i];
      itkVectorImage->SetOrigin( origin );     // Set the image origin

      VectorImageType::DirectionType direction;
      if (bruker)
      {
        MITK_INFO << "Setting image matrix to identity (bruker hack!!!)";
        direction.SetIdentity();
      }
      else
      {
        ImageType4D::DirectionType direction4 = img4->GetDirection();
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
            direction[i][j] = direction4[i][j];
      }
      itkVectorImage->SetDirection( direction );  // Set the image direction

      VectorImageType::RegionType region;
      ImageType4D::RegionType region4 = img4->GetLargestPossibleRegion();

      VectorImageType::RegionType::SizeType size;
      ImageType4D::RegionType::SizeType size4 = region4.GetSize();

      for(int i=0; i<3; i++)
        size[i] = size4[i];

      VectorImageType::RegionType::IndexType index;
      ImageType4D::RegionType::IndexType index4 = region4.GetIndex();
      for(int i=0; i<3; i++)
        index[i] = index4[i];

      region.SetSize(size);
      region.SetIndex(index);
      itkVectorImage->SetRegions( region );

      itkVectorImage->SetVectorLength(size4[3]);
      itkVectorImage->Allocate();

      itk::ImageRegionIterator<VectorImageType>   it ( itkVectorImage,  itkVectorImage->GetLargestPossibleRegion() );
      typedef VectorImageType::PixelType VecPixType;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        VecPixType vec = it.Get();
        VectorImageType::IndexType currentIndex = it.GetIndex();
        for(int i=0; i<3; i++)
          index4[i] = currentIndex[i];
        for(unsigned int ind=0; ind<vec.Size(); ind++)
        {
          index4[3] = ind;
          vec[ind] = img4->GetPixel(index4);
        }
        it.Set(vec);
      }

      // Diffusion Image information START
      GradientDirectionContainerType::ConstPointer DiffusionVectors;
      MeasurementFrameType MeasurementFrame;
      MeasurementFrame.set_identity();
      double BValue = -1;
      // Diffusion Image information END

      if (isNifti)
      {
        std::string base_path = itksys::SystemTools::GetFilenamePath(this->GetInputLocation());
        std::string base = mimeType.GetFilenameWithoutExtension(this->GetInputLocation());
        if (!base_path.empty())
        {
          base = base_path + "/" + base;
          base_path += "/";
        }

        // check for possible file names
        std::string bvals_file, bvecs_file;
        if (itksys::SystemTools::FileExists(base+".bvals"))
          bvals_file = base+".bvals";
        else if (itksys::SystemTools::FileExists(base+".bval"))
          bvals_file = base+".bval";
        else if (itksys::SystemTools::FileExists(base_path+"bvals"))
          bvals_file = base_path + "bvals";
        else if (itksys::SystemTools::FileExists(base_path+"bval"))
          bvals_file = base_path + "bval";

        if (itksys::SystemTools::FileExists(base+".bvecs"))
          bvecs_file = base+".bvecs";
        else if (itksys::SystemTools::FileExists(base+".bvec"))
          bvecs_file = base+".bvec";
        else if (itksys::SystemTools::FileExists(base_path+"bvecs"))
          bvecs_file = base_path + "bvecs";
        else if (itksys::SystemTools::FileExists(base_path+"bvec"))
          bvecs_file = base_path + "bvec";

        DiffusionVectors = mitk::gradients::ReadBvalsBvecs(bvals_file, bvecs_file, BValue);
        MeasurementFrame = mitk::gradients::FslBvecsMeasurementFrame(itkVectorImage->GetDirection());
      }
      else
      {
        MITK_INFO << "Parsing bruker method file for gradient information.";
        std::string base_path = itksys::SystemTools::GetFilenamePath(this->GetInputLocation());
        std::string base = mimeType.GetFilenameWithoutExtension(this->GetInputLocation());
        base_path += "/../../";
        std::string method_file = base_path + "method";

        int g_count = 0;
        int num_gradients = 0;
        int num_b = 0;
        std::vector<double> grad_values;
        std::vector<double> b_values;

        std::fstream newfile;
        newfile.open(method_file,ios::in);
        std::string line;
        while(getline(newfile, line))
        {
          std::vector<std::string> result;
          boost::split(result, line, boost::is_any_of("="));

          // check if current section is b-value section
          if (result.size()==2 && result.at(0)=="##$PVM_DwEffBval")
          {
            std::vector<std::string> vec_spec;
            boost::split(vec_spec, result.at(1), boost::is_any_of("( )"), boost::token_compress_on);
            for (auto a : vec_spec)
              if (!a.empty())
              {
                num_b = boost::lexical_cast<int>(a);
                MITK_INFO << "Found " << num_b << " b-values";
                break;
              }
            continue;
          }

          // check if current section is gradient vector value section
          if (result.size()==2 && result.at(0)=="##$PVM_DwGradVec")
          {
            std::vector<std::string> vec_spec;
            boost::split(vec_spec, result.at(1), boost::is_any_of("( ,)"), boost::token_compress_on);
            for (auto a : vec_spec)
              if (!a.empty())
              {
                num_gradients = boost::lexical_cast<int>(a);
                MITK_INFO << "Found " << num_gradients << " gradients";
                break;
              }
            continue;
          }

          // get gradient vector values
          if (num_gradients>0)
          {

            std::vector<std::string> grad_values_line;
            boost::split(grad_values_line, line, boost::is_any_of(" "));
            for (auto a : grad_values_line)
              if (!a.empty())
              {
                grad_values.push_back(boost::lexical_cast<double>(a));
                ++g_count;
                if (g_count%3==0)
                  --num_gradients;
              }
          }

          // get b-values
          if (num_b>0)
          {
            std::vector<std::string> b_values_line;
            boost::split(b_values_line, line, boost::is_any_of(" "));
            for (auto a : b_values_line)
              if (!a.empty())
              {
                b_values.push_back(boost::lexical_cast<double>(a));
                if (b_values.back()>BValue)
                  BValue = b_values.back();
                num_b--;
              }
          }
        }
        newfile.close();

        GradientDirectionContainerType::Pointer temp = GradientDirectionContainerType::New();
        if (!b_values.empty() && grad_values.size()==b_values.size()*3)
        {
//          MITK_INFO << "Switching gradient vector x and y (bruker hack!!!)";
          for(unsigned int i=0; i<b_values.size(); ++i)
          {
            mitk::DiffusionPropertyHelper::GradientDirectionType vec;
            vec[0] = grad_values.at(i*3);
            vec[1] = grad_values.at(i*3+1);
            vec[2] = grad_values.at(i*3+2);

            double b_val = b_values.at(i);

            // Adjust the vector length to encode gradient strength
            if (BValue>0)
            {
              double factor = b_val/BValue;
              if(vec.magnitude() > 0)
              {
                vec.normalize();
                vec[0] = sqrt(factor)*vec[0];
                vec[1] = sqrt(factor)*vec[1];
                vec[2] = sqrt(factor)*vec[2];
              }
            }

            temp->InsertElement(i,vec);
          }
        }
        else
        {
          mitkThrow() << "No valid gradient information found.";
        }
        DiffusionVectors = temp;
      }

      if (DiffusionVectors->Size() != itkVectorImage->GetVectorLength())
        mitkThrow() << "Found " << DiffusionVectors->Size() << " gradients for " << itkVectorImage->GetVectorLength() << " volumes";

      outputForCache = mitk::GrabItkImageMemory( itkVectorImage);

      // create BValueMap
      mitk::BValueMapProperty::BValueMap BValueMap = mitk::BValueMapProperty::CreateBValueMap(DiffusionVectors,BValue);
      mitk::DiffusionPropertyHelper::SetOriginalGradientContainer(outputForCache, DiffusionVectors);
      mitk::DiffusionPropertyHelper::SetMeasurementFrame(outputForCache, MeasurementFrame);
      mitk::DiffusionPropertyHelper::SetBValueMap(outputForCache, BValueMap);
      mitk::DiffusionPropertyHelper::SetReferenceBValue(outputForCache, BValue);
      // The measurement frame already contains the image rotation (FSL bvecs; Bruker images have an
      // identity matrix), so the modified gradients are world vectors.
      mitk::DiffusionPropertyHelper::SetApplyMatrixToGradients(outputForCache, false);
      mitk::DiffusionPropertyHelper::InitializeImage(outputForCache);

      // Since we have already read the tree, we can store it in a cache variable
      // so that it can be assigned to the DataObject in GenerateData();
      m_OutputCache = outputForCache;
      m_CacheTime.Modified();
    }
    catch(std::exception& e)
    {
      MITK_INFO << "Std::Exception while reading file!!";
      MITK_INFO << e.what();
      throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
    }
    catch(...)
    {
      MITK_INFO << "Exception while reading file!!";
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested file!");
    }
  }
}

} //namespace MITK

#endif
