/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkNrrdTensorImageReader.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"

#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include <itkNiftiImageIO.h>

#include "mitkITKImageImport.h"
#include "mitkImageDataItem.h"
#include <mitkLocaleSwitch.h>
#include <mitkIOUtil.h>
#include <itksys/SystemTools.hxx>

#include <fstream>
#include <optional>

namespace
{
  bool HasNrrdMagic(const std::string& file)
  {
    std::ifstream stream(file, std::ios::binary);
    char magic[4] = {};
    stream.read(magic, 4);
    return stream && std::string(magic, 4) == "NRRD";
  }

  /** A copy of a file with a NIfTI extension, which niftilib needs; removed when it goes out of scope. */
  class TemporaryNiftiCopy
  {
  public:
    explicit TemporaryNiftiCopy(const std::string& file)
      : m_Path(mitk::IOUtil::CreateTemporaryFile("temp_dti_XXXXXX.nii.gz"))
    {
      itksys::SystemTools::CopyAFile(file, m_Path);
    }

    ~TemporaryNiftiCopy()
    {
      itksys::SystemTools::RemoveFile(m_Path);
    }

    TemporaryNiftiCopy(const TemporaryNiftiCopy&) = delete;
    TemporaryNiftiCopy& operator=(const TemporaryNiftiCopy&) = delete;

    const std::string& GetPath() const { return m_Path; }

  private:
    std::string m_Path;
  };

  using FrameType = mitk::NrrdTensorImageReader::MeasurementFrameType;

  /**
   * Returns the rotation M from the measurement frame to world coordinates, so that the world
   * tensor is M * T * M^T. Both the standard NRRD field and MITK's own "measurement frame" key hold
   * the frame vectors F as rows, and M = F^T, the same convention as for DWI gradients
   * (modified = MF^T * original). No MITK version wrote a frame into tensor files, so MITK's key
   * needs no legacy direction. No frame or an all-zero frame gives nothing.
   */
  std::optional<FrameType> ReadFrameToWorld(const itk::MetaDataDictionary& dictionary)
  {
    std::optional<FrameType> frame;

    // ITK stores the standard NRRD field (written by Slicer, DTIPrep, ITK) as a matrix, not as a string.
    std::vector<std::vector<double>> rows;
    if (itk::ExposeMetaData(dictionary, "NRRD_measurement frame", rows) && rows.size() == 3 && rows[0].size() == 3 && rows[1].size() == 3 && rows[2].size() == 3)
    {
      FrameType standardFrame;
      for (unsigned int r = 0; r < 3; ++r)
        for (unsigned int c = 0; c < 3; ++c)
          standardFrame(r, c) = rows[r][c];
      frame = standardFrame;
    }

    std::string value;
    if (!frame && itk::ExposeMetaData<std::string>(dictionary, "measurement frame", value))
    {
      FrameType mitkFrame;
      if (sscanf(value.c_str(), " ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ",
                 &mitkFrame(0, 0), &mitkFrame(0, 1), &mitkFrame(0, 2), &mitkFrame(1, 0), &mitkFrame(1, 1), &mitkFrame(1, 2), &mitkFrame(2, 0), &mitkFrame(2, 1), &mitkFrame(2, 2)) != 9)
        mitkThrow() << "Invalid measurement frame " << value;
      frame = mitkFrame;
    }

    if (!frame || frame->GetVnlMatrix().absolute_value_max() <= 1e-9)
      return std::nullopt;
    return FrameType(frame->GetTranspose());
  }
}

namespace mitk
{

  NrrdTensorImageReader::NrrdTensorImageReader(const NrrdTensorImageReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  NrrdTensorImageReader::NrrdTensorImageReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::DTI_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::DTI_MIMETYPE_DESCRIPTION() )
  {
    m_ServiceReg = this->RegisterService();
  }

  NrrdTensorImageReader::~NrrdTensorImageReader()
  {
  }

  std::vector<itk::SmartPointer<BaseData> > NrrdTensorImageReader::DoRead()
  {
    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    std::string location = GetInputLocation();

    if ( location == "")
    {
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename is empty!");
    }
    else
    {
      try
      {
        mitk::LocaleSwitch localeSwitch("C");

        try
        {
          // Tensor files written by MITK are NRRD and go straight to the NRRD branch below; only
          // other files are tried as NIfTI, through a temporary copy.
          if (HasNrrdMagic(location))
            throw itk::ImageFileReaderException(__FILE__, __LINE__, "NRRD file");

          const TemporaryNiftiCopy niftiCopy(location);
          const std::string fname3 = niftiCopy.GetPath();

          typedef itk::VectorImage<float,3> ImageType;
          itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
          typedef itk::ImageFileReader<ImageType> FileReaderType;
          FileReaderType::Pointer reader = FileReaderType::New();
          reader->SetImageIO(io);
          reader->SetFileName(fname3);
          reader->Update();
          ImageType::Pointer img = reader->GetOutput();

          TensorImage::ItkTensorImageType::Pointer vecImg = TensorImage::ItkTensorImageType::New();
          vecImg->SetSpacing( img->GetSpacing() );   // Set the image spacing
          vecImg->SetOrigin( img->GetOrigin() );     // Set the image origin
          vecImg->SetDirection( img->GetDirection() );  // Set the image direction
          vecImg->SetRegions( img->GetLargestPossibleRegion());
          vecImg->Allocate();

          itk::ImageRegionIterator<TensorImage::ItkTensorImageType> ot (vecImg, vecImg->GetLargestPossibleRegion() );
          ot.GoToBegin();

          itk::ImageRegionIterator<ImageType> it (img, img->GetLargestPossibleRegion() );
          it.GoToBegin();

          typedef ImageType::PixelType  VarPixType;
          typedef TensorImage::PixelType FixPixType;
          int numComponents = img->GetNumberOfComponentsPerPixel();

          if (numComponents==6)
          {
            MITK_INFO << "Trying to load dti as 6-comp nifti ...";
            while (!it.IsAtEnd())
            {
              VarPixType vec = it.Get();
              FixPixType fixVec(vec.GetDataPointer());

              TensorImage::PixelType tensor;
              tensor.SetElement(0, vec.GetElement(0));
              tensor.SetElement(1, vec.GetElement(1));
              tensor.SetElement(2, vec.GetElement(2));
              tensor.SetElement(3, vec.GetElement(3));
              tensor.SetElement(4, vec.GetElement(4));
              tensor.SetElement(5, vec.GetElement(5));

              fixVec = tensor;

              ot.Set(fixVec);
              ++ot;
              ++it;
            }
          }
          else if(numComponents==9)
          {
            MITK_INFO << "Trying to load dti as 9-comp nifti ...";
            while (!it.IsAtEnd())
            {
              VarPixType vec = it.Get();
              TensorImage::PixelType tensor;
              tensor.SetElement(0, vec.GetElement(0));
              tensor.SetElement(1, vec.GetElement(1));
              tensor.SetElement(2, vec.GetElement(2));
              tensor.SetElement(3, vec.GetElement(4));
              tensor.SetElement(4, vec.GetElement(5));
              tensor.SetElement(5, vec.GetElement(8));

              FixPixType fixVec(tensor);
              ot.Set(fixVec);
              ++ot;
              ++it;
            }
          }
          else if (numComponents==1)
          {
            MITK_INFO << "Trying to load dti as 4D nifti ...";
            typedef itk::Image<float,4> ImageType;
            typedef itk::ImageFileReader<ImageType> FileReaderType;
            FileReaderType::Pointer reader = FileReaderType::New();
            reader->SetImageIO(io);
            reader->SetFileName(fname3);
            reader->Update();
            ImageType::Pointer img = reader->GetOutput();

            itk::Size<4> size = img->GetLargestPossibleRegion().GetSize();

            while (!ot.IsAtEnd())
            {
              TensorImage::PixelType tensor;
              ImageType::IndexType idx;
              idx[0] = ot.GetIndex()[0]; idx[1] = ot.GetIndex()[1]; idx[2] = ot.GetIndex()[2];

              if (size[3]==6)
              {
                for (unsigned int te=0; te<size[3]; te++)
                {
                  idx[3] = te;
                  tensor.SetElement(te, img->GetPixel(idx));
                }
              }
              else if (size[3]==9)
              {
                idx[3] = 0;
                tensor.SetElement(0, img->GetPixel(idx));
                idx[3] = 1;
                tensor.SetElement(1, img->GetPixel(idx));
                idx[3] = 2;
                tensor.SetElement(2, img->GetPixel(idx));
                idx[3] = 4;
                tensor.SetElement(3, img->GetPixel(idx));
                idx[3] = 5;
                tensor.SetElement(4, img->GetPixel(idx));
                idx[3] = 8;
                tensor.SetElement(5, img->GetPixel(idx));
              }
              else
                throw itk::ImageFileReaderException(__FILE__, __LINE__, "Unknown number of components for DTI file. Should be 6 or 9!");

              FixPixType fixVec(tensor);
              ot.Set(fixVec);
              ++ot;
            }
          }
          OutputType::Pointer resultImage = OutputType::New();
          resultImage->InitializeByItk( vecImg.GetPointer() );
          resultImage->SetVolume( vecImg->GetBufferPointer() );
          result.push_back( resultImage.GetPointer() );
        }
        catch(...)
        {
          MITK_INFO << "Trying to load dti as nrrd ...";

          typedef itk::VectorImage<float,3> ImageType;
          itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
          typedef itk::ImageFileReader<ImageType> FileReaderType;
          FileReaderType::Pointer reader = FileReaderType::New();
          reader->SetImageIO(io);
          reader->SetFileName(location);
          reader->Update();
          ImageType::Pointer img = reader->GetOutput();

          TensorImage::ItkTensorImageType::Pointer vecImg = TensorImage::ItkTensorImageType::New();
          vecImg->SetSpacing( img->GetSpacing() );   // Set the image spacing
          vecImg->SetOrigin( img->GetOrigin() );     // Set the image origin
          vecImg->SetDirection( img->GetDirection() );  // Set the image direction
          vecImg->SetRegions( img->GetLargestPossibleRegion());
          vecImg->Allocate();

          itk::ImageRegionIterator<TensorImage::ItkTensorImageType> ot (vecImg, vecImg->GetLargestPossibleRegion() );
          ot.GoToBegin();

          itk::ImageRegionIterator<ImageType> it (img, img->GetLargestPossibleRegion() );
          it.GoToBegin();

          typedef ImageType::PixelType  VarPixType;
          typedef TensorImage::PixelType FixPixType;
          int numComponents = img->GetNumberOfComponentsPerPixel();

          const std::optional<FrameType> frameToWorld = ReadFrameToWorld(img->GetMetaDataDictionary());

          if (numComponents==6)
          {
            while (!it.IsAtEnd())
            {
              VarPixType vec = it.Get();
              FixPixType fixVec(vec.GetDataPointer());

              if (frameToWorld)
                fixVec = fixVec.Rotate(*frameToWorld);

              ot.Set(fixVec);
              ++ot;
              ++it;
            }
          }
          else if(numComponents==9)
          {
            while (!it.IsAtEnd())
            {
              VarPixType vec = it.Get();
              TensorImage::PixelType tensor;
              tensor.SetElement(0, vec.GetElement(0));
              tensor.SetElement(1, vec.GetElement(1));
              tensor.SetElement(2, vec.GetElement(2));
              tensor.SetElement(3, vec.GetElement(4));
              tensor.SetElement(4, vec.GetElement(5));
              tensor.SetElement(5, vec.GetElement(8));

              if (frameToWorld)
                tensor = tensor.Rotate(*frameToWorld);

              FixPixType fixVec(tensor);
              ot.Set(fixVec);
              ++ot;
              ++it;
            }
          }
          else if (numComponents==1)
          {
            typedef itk::Image<float,4> ImageType;
            itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
            typedef itk::ImageFileReader<ImageType> FileReaderType;
            FileReaderType::Pointer reader = FileReaderType::New();
            reader->SetImageIO(io);
            reader->SetFileName(location);
            reader->Update();
            ImageType::Pointer img = reader->GetOutput();

            itk::Size<4> size = img->GetLargestPossibleRegion().GetSize();

            while (!ot.IsAtEnd())
            {
              TensorImage::PixelType tensor;
              ImageType::IndexType idx;
              idx[0] = ot.GetIndex()[0]; idx[1] = ot.GetIndex()[1]; idx[2] = ot.GetIndex()[2];

              if (size[3]==6)
              {
                for (unsigned int te=0; te<size[3]; te++)
                {
                  idx[3] = te;
                  tensor.SetElement(te, img->GetPixel(idx));
                }
              }
              else if (size[3]==9)
              {
                idx[3] = 0;
                tensor.SetElement(0, img->GetPixel(idx));
                idx[3] = 1;
                tensor.SetElement(1, img->GetPixel(idx));
                idx[3] = 2;
                tensor.SetElement(2, img->GetPixel(idx));
                idx[3] = 4;
                tensor.SetElement(3, img->GetPixel(idx));
                idx[3] = 5;
                tensor.SetElement(4, img->GetPixel(idx));
                idx[3] = 8;
                tensor.SetElement(5, img->GetPixel(idx));
              }
              else
                throw itk::ImageFileReaderException(__FILE__, __LINE__, "Unknown number of komponents for DTI file. Should be 6 or 9!");

              if (frameToWorld)
                tensor = tensor.Rotate(*frameToWorld);
              FixPixType fixVec(tensor);
              ot.Set(fixVec);
              ++ot;
            }
          }
          else
          {
            throw itk::ImageFileReaderException(__FILE__, __LINE__, "Image has wrong number of pixel components!");
          }

          OutputType::Pointer resultImage = OutputType::New();
          resultImage->InitializeByItk( vecImg.GetPointer() );
          resultImage->SetVolume( vecImg->GetBufferPointer() );
          result.push_back( resultImage.GetPointer() );
        }

      }
      catch(std::exception& e)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
      }
      catch(...)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested DTI file!");
      }
    }

    return result;
  }

} //namespace MITK

mitk::NrrdTensorImageReader* mitk::NrrdTensorImageReader::Clone() const
{
  return new NrrdTensorImageReader(*this);
}
