/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkShImageReader.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"

#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkNiftiImageIO.h"
#include "mitkITKImageImport.h"
#include "mitkImageDataItem.h"
#include <mitkLocaleSwitch.h>
#include <itkShCoefficientImageImporter.h>

#include <cmath>
#include <optional>

namespace
{
  // Constant initialization: the reader is constructed while the library is loaded, possibly
  // before dynamically initialized globals of this file.
  constexpr const char *CONVENTION_OPTION = "SH convention";
  constexpr const char *AUTO = "auto (file marker, else MITK legacy)";
  constexpr const char *MRTRIX3 = "MRtrix 3";
  constexpr const char *MITK_LEGACY = "MITK legacy (MRtrix 0.2)";
  constexpr const char *FSL_DIPY = "FSL/Dipy";

  /** How the coefficients in the file are defined. */
  enum class FileBasis
  {
    MitkLegacy, // MITK's internal basis (MRtrix 0.2), taken over unchanged
    Mrtrix3,    // MRtrix 3: orthonormal, m != 0 larger by sqrt(2) than MITK's
    MitkFsl     // MITK's FSL definition (ShImage::FSL), taken over unchanged
  };

  /** Convention marked by ShImageWriter (NRRD key or NIfTI descrip); it overrides the option. */
  std::optional<FileBasis> MitkMarker(const itk::MetaDataDictionary &dictionary)
  {
    std::string value;
    if (itk::ExposeMetaData<std::string>(dictionary, "MITK_SH_convention", value)
        || (itk::ExposeMetaData<std::string>(dictionary, "ITK_FileNotes", value) && value.rfind("MITK_SH_convention=", 0) == 0))
    {
      value = value.substr(value.find('=') + 1);
      if (value == "FSL")
        return FileBasis::MitkFsl;
      if (value == "MRtrix3")
        return FileBasis::Mrtrix3;
    }
    return std::nullopt;
  }

  /** MRtrix writes "MRtrix version: <version>" into the NIfTI descrip; version 3 and later use the orthonormal basis. */
  bool WrittenByMrtrix3(const itk::MetaDataDictionary &dictionary)
  {
    const std::string prefix = "MRtrix version: ";
    std::string notes;
    if (!itk::ExposeMetaData<std::string>(dictionary, "ITK_FileNotes", notes) || notes.rfind(prefix, 0) != 0)
      return false;
    try
    {
      return std::stoi(notes.substr(prefix.size())) >= 3;
    }
    catch (const std::exception &)
    {
      return false;
    }
  }
}

namespace mitk
{
  ShImageReader::ShImageReader(const ShImageReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  ShImageReader::ShImageReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::SH_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::SH_MIMETYPE_DESCRIPTION() )
  {
    Options defaultOptions;
    defaultOptions[CONVENTION_OPTION] = std::string(AUTO);
    defaultOptions[std::string(CONVENTION_OPTION) + ".enum"] = std::vector<std::string>{AUTO, MRTRIX3, MITK_LEGACY, FSL_DIPY};
    this->SetDefaultOptions(defaultOptions);

    this->RegisterService();
  }

  ShImageReader::~ShImageReader() = default;

  template <int shOrder>
  mitk::Image::Pointer ShImageReader::ConvertShImage(ShImage::ShOnDiskType::Pointer img)
  {
    typedef itk::ShCoefficientImageImporter< float, shOrder > ImporterType;
    typename ImporterType::Pointer importer = ImporterType::New();
    importer->SetInputImage(img);
    importer->GenerateData();

    mitk::ShImage::Pointer shImage = mitk::ShImage::New();

    const auto &dictionary = img->GetMetaDataDictionary();
    auto basis = MitkMarker(dictionary);
    if (!basis)
    {
      const std::string option = us::any_cast<std::string>(this->GetOptions()[CONVENTION_OPTION]);
      if (option == MRTRIX3 || (option == AUTO && WrittenByMrtrix3(dictionary)))
        basis = FileBasis::Mrtrix3;
      else if (option == FSL_DIPY)
        basis = FileBasis::MitkFsl;
      else
        basis = FileBasis::MitkLegacy;
    }

    if (*basis == FileBasis::MitkFsl)
      shImage->SetShConvention(mitk::ShImage::SH_CONVENTION::FSL);

    if (*basis == FileBasis::Mrtrix3)
    {
      // MITK's internal basis is MRtrix 0.2's, whose functions with m != 0 are MRtrix 3's divided
      // by sqrt(2); m = 0 sits at index l(l+1)/2 of each even order l.
      auto coefficients = importer->GetCoefficientImage();
      itk::ImageRegionIterator<typename ImporterType::CoefficientImageType> it(coefficients, coefficients->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
      {
        auto pixel = it.Get();
        for (int l = 2; l <= shOrder; l += 2)
          for (int m = -l; m <= l; ++m)
            if (m != 0)
              pixel[l * (l + 1) / 2 + m] = static_cast<float>(pixel[l * (l + 1) / 2 + m] * std::sqrt(2.0));
        it.Set(pixel);
      }
    }

    // The SH image takes over the buffer of the coefficient image.
    mitk::Image::Pointer resultImage = dynamic_cast<mitk::Image*>(shImage.GetPointer());
    mitk::GrabItkImageMemoryChannel(importer->GetCoefficientImage().GetPointer(), nullptr, resultImage.GetPointer(), false);
    return resultImage;
  }

  std::vector<itk::SmartPointer<BaseData> > ShImageReader::DoRead()
  {
    mitk::LocaleSwitch localeSwitch("C");
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
        std::string ext = itksys::SystemTools::GetFilenameLastExtension(location);
        typedef itk::ImageFileReader< ShImage::ShOnDiskType > FileReaderType;
        FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetFileName(location);
        if (ext==".shi")
        {
          itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
          reader->SetImageIO(io);
        }
        reader->Update();
        ShImage::ShOnDiskType::Pointer img = reader->GetOutput();

        switch (img->GetLargestPossibleRegion().GetSize()[3])
        {
        case 6:
          result.push_back( ConvertShImage<2>(img).GetPointer() );
          break;
        case 15:
          result.push_back( ConvertShImage<4>(img).GetPointer() );
          break;
        case 28:
          result.push_back( ConvertShImage<6>(img).GetPointer() );
          break;
        case 45:
          result.push_back( ConvertShImage<8>(img).GetPointer() );
          break;
        case 66:
          result.push_back( ConvertShImage<10>(img).GetPointer() );
          break;
        case 91:
          result.push_back( ConvertShImage<12>(img).GetPointer() );
          break;
        default :
          mitkThrow() << "SH order larger 12 not supported";
        }
      }
      catch(std::exception& e)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
      }
      catch(...)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested SH file!");
      }
    }
    return result;
  }

} //namespace MITK

mitk::ShImageReader* mitk::ShImageReader::Clone() const
{
  return new ShImageReader(*this);
}
