/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkPeakImageReader.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"
#include <mitkPeakImage.h>
#include <mitkIOUtil.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <itkImageFileReader.h>
#include <itksys/SystemTools.hxx>
#include <itkNiftiImageIO.h>
#include <itkNrrdImageIO.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkLocaleSwitch.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMetaDataObject.h>

namespace
{
  // constexpr: the reader is constructed by the module activator during static initialization.
  constexpr const char* ConvertMrtrixOption = "Convert MRtrix peaks (RAS) to MITK (LPS)";

  /**
   * MRtrix (sh2peaks) writes peak directions in scanner RAS coordinates, MITK uses world LPS
   * vectors. MRtrix files are recognized by their NIfTI descrip. The NaN that MRtrix writes for
   * missing peaks is already set to 0 by niftilib.
   */
  void ConvertMrtrixPeaks(mitk::PeakImage::ItkPeakImageType* peaks)
  {
    std::string description;
    if (!itk::ExposeMetaData<std::string>(peaks->GetMetaDataDictionary(), "ITK_FileNotes", description) || !description.starts_with("MRtrix"))
      return;

    itk::ImageRegionIteratorWithIndex<mitk::PeakImage::ItkPeakImageType> it(peaks, peaks->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      if (it.GetIndex()[3] % 3 != 2)
        it.Set(-it.Get());
    }
  }
}

namespace mitk
{

  PeakImageReader::PeakImageReader(const PeakImageReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  PeakImageReader::PeakImageReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::PEAK_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::PEAK_MIMETYPE_DESCRIPTION() )
  {
    Options defaultOptions;
    defaultOptions[ConvertMrtrixOption] = false;
    this->SetDefaultOptions(defaultOptions);

    m_ServiceReg = this->RegisterService();
  }

  PeakImageReader::~PeakImageReader()
  {
  }

  std::vector<itk::SmartPointer<BaseData> > PeakImageReader::DoRead()
  {
    mitk::LocaleSwitch localeSwitch("C");
    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    std::string location = GetInputLocation();
    std::string ext = itksys::SystemTools::GetFilenameLastExtension(location);
    MITK_INFO << "Reading " << location;

    typedef itk::ImageFileReader<PeakImage::ItkPeakImageType> FileReaderType;
    FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName(location);
    if (ext==".peak" || ext==".nrrd")
    {
      itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
      reader->SetImageIO(io);
    }
    else
    {
      itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
      reader->SetImageIO(io);
    }
    reader->Update();
    if (us::any_cast<bool>(this->GetOptions()[ConvertMrtrixOption]))
      ConvertMrtrixPeaks(reader->GetOutput());

    Image::Pointer resultImage = dynamic_cast<Image*>(PeakImage::New().GetPointer());
    mitk::CastToMitkImage(reader->GetOutput(), resultImage);
    resultImage->SetVolume(reader->GetOutput()->GetBufferPointer());

    StringProperty::Pointer nameProp;
    nameProp = StringProperty::New(itksys::SystemTools::GetFilenameWithoutExtension(GetInputLocation()));
    resultImage->SetProperty("name", nameProp);
    dynamic_cast<PeakImage*>(resultImage.GetPointer())->ConstructPolydata();

    result.push_back( resultImage.GetPointer() );

    return result;
  }

} //namespace MITK

mitk::PeakImageReader* mitk::PeakImageReader::Clone() const
{
  return new PeakImageReader(*this);
}
