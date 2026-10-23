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

  using TensorImageType = mitk::TensorImage::ItkTensorImageType;
  using TensorType = mitk::TensorImage::PixelType;

  /**
   * The tensor of 6 values (xx, xy, xz, yy, yz, zz) or of a full 3x3 matrix in row-major order,
   * stride values apart.
   */
  TensorType ToTensor(const float *values, unsigned int numValues, std::size_t stride)
  {
    const unsigned int fromSix[6] = {0, 1, 2, 3, 4, 5};
    const unsigned int fromNine[6] = {0, 1, 2, 4, 5, 8};
    const unsigned int *elements = numValues == 9 ? fromNine : fromSix;
    TensorType tensor;
    for (unsigned int i = 0; i < 6; ++i)
      tensor[i] = values[elements[i] * stride];
    return tensor;
  }

  /**
   * Reads a tensor file with a new IO of the given type: a 3D image of 6 or 9 components per voxel,
   * or a 4D image of 6 or 9 volumes. The metadata of the file go to dictionary.
   */
  template <class TImageIO>
  TensorImageType::Pointer ReadTensors(const std::string &file, itk::MetaDataDictionary &dictionary)
  {
    using VectorImageType = itk::VectorImage<float, 3>;
    auto reader = itk::ImageFileReader<VectorImageType>::New();
    reader->SetImageIO(TImageIO::New());
    reader->SetFileName(file);
    reader->Update();
    const VectorImageType::Pointer vectors = reader->GetOutput();
    dictionary = vectors->GetMetaDataDictionary();

    auto tensors = TensorImageType::New();
    tensors->SetSpacing(vectors->GetSpacing());
    tensors->SetOrigin(vectors->GetOrigin());
    tensors->SetDirection(vectors->GetDirection());
    tensors->SetRegions(vectors->GetLargestPossibleRegion());
    tensors->Allocate();
    itk::ImageRegionIterator<TensorImageType> ot(tensors, tensors->GetLargestPossibleRegion());

    const unsigned int numComponents = vectors->GetNumberOfComponentsPerPixel();
    if (numComponents == 6 || numComponents == 9)
    {
      itk::ImageRegionConstIterator<VectorImageType> it(vectors, vectors->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it, ++ot)
        ot.Set(ToTensor(it.Get().GetDataPointer(), numComponents, 1));
      return tensors;
    }
    if (numComponents != 1)
      mitkThrow() << "Unknown number of components for DTI file: " << numComponents << ". Should be 6 or 9!";

    // A 4D file of one volume per tensor element.
    using VolumesType = itk::Image<float, 4>;
    auto volumeReader = itk::ImageFileReader<VolumesType>::New();
    volumeReader->SetImageIO(TImageIO::New());
    volumeReader->SetFileName(file);
    volumeReader->Update();
    const VolumesType::Pointer volumes = volumeReader->GetOutput();
    const auto size = volumes->GetLargestPossibleRegion().GetSize();
    if (size[3] != 6 && size[3] != 9)
      mitkThrow() << "Unknown number of volumes for DTI file: " << size[3] << ". Should be 6 or 9!";
    const std::size_t volumeSize = size[0] * size[1] * size[2];
    if (volumeSize != tensors->GetLargestPossibleRegion().GetNumberOfPixels())
      mitkThrow() << "The volumes of the DTI file do not match its first volume.";

    const float *values = volumes->GetBufferPointer();
    for (std::size_t voxel = 0; !ot.IsAtEnd(); ++ot, ++voxel)
      ot.Set(ToTensor(values + voxel, size[3], volumeSize));
    return tensors;
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
    this->RegisterService();
  }

  NrrdTensorImageReader::~NrrdTensorImageReader() = default;

  std::vector<itk::SmartPointer<BaseData> > NrrdTensorImageReader::DoRead()
  {
    const std::string location = GetInputLocation();
    if (location.empty())
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename is empty!");

    TensorImageType::Pointer tensors;
    try
    {
      mitk::LocaleSwitch localeSwitch("C");
      itk::MetaDataDictionary dictionary;
      if (HasNrrdMagic(location))
      {
        tensors = ReadTensors<itk::NrrdImageIO>(location, dictionary);
        if (const auto frameToWorld = ReadFrameToWorld(dictionary))
        {
          itk::ImageRegionIterator<TensorImageType> it(tensors, tensors->GetLargestPossibleRegion());
          for (; !it.IsAtEnd(); ++it)
            it.Set(it.Get().Rotate(*frameToWorld));
        }
      }
      else
      {
        // Other files are tried as NIfTI, through a copy with a NIfTI extension, which niftilib
        // needs; they reach this reader only when renamed to .dti. Six volumes or components are
        // taken in FSL dtifit's order (xx, xy, xz, yy, yz, zz), in the voxel frame as stored: FSL's
        // flip of the voxel frame (applied to bvecs, see DiffusionImage's NIfTI reader) and the
        // image rotation are not applied. MRtrix's dwi2tensor order (D11, D22, D33, D12, D13, D23,
        // RAS) is not recognized and would be read in the wrong order.
        const TemporaryNiftiCopy niftiCopy(location);
        tensors = ReadTensors<itk::NiftiImageIO>(niftiCopy.GetPath(), dictionary);
      }
    }
    catch(const std::exception& e)
    {
      throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
    }
    catch(...)
    {
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested DTI file!");
    }

    auto resultImage = OutputType::New();
    mitk::GrabItkImageMemory(tensors.GetPointer(), resultImage.GetPointer());
    return { resultImage.GetPointer() };
  }

} //namespace MITK

mitk::NrrdTensorImageReader* mitk::NrrdTensorImageReader::Clone() const
{
  return new NrrdTensorImageReader(*this);
}
