/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionDICOMFileReader.h"

#include "mitkDiffusionDICOMFileReaderHelper.h"

#include "mitkDiffusionHeaderSiemensDICOMFileReader.h"
#include "mitkDiffusionHeaderSiemensMosaicDICOMFileReader.h"
#include "mitkDiffusionHeaderGEDICOMFileReader.h"
#include "mitkDiffusionHeaderPhilipsDICOMFileReader.h"
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>

#include "mitkStringProperty.h"
#include <mitkImageCast.h>
#include <mitkDICOMTagBasedSorter.h>
#include <mitkDICOMSortByTag.h>
#include <mitkSortByImagePositionPatient.h>
#include <mitkLocaleSwitch.h>

static void PerformHeaderAnalysis( mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType headers )
{
  unsigned int images = headers.size();

  unsigned int unweighted_images = 0;
  unsigned int weighted_images = 0;

  mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType::const_iterator c_iter = headers.begin();
  while( c_iter != headers.end() )
  {
    const mitk::DiffusionImageDICOMHeaderInformation h = *c_iter;
    if( h.baseline ) unweighted_images++;
    if( h.b_value > 0 ) weighted_images++;

    ++c_iter;
  }

  MITK_INFO << "  :: Analyzed volumes " << images  << "\n"
            << "  :: \t"<< unweighted_images << " b = 0" << "\n"
            << "  :: \t"<< weighted_images << " b > 0";
}

/** A DICOM string value without the trailing space or zero byte that pads odd-length values to
 *  even length; gdcm::Scanner returns the values padded, and the series description becomes the
 *  image name. */
static std::string TrimDicomPadding(const char* value)
{
  std::string text(value);
  const auto end = text.find_last_not_of(std::string(" \0", 2));
  return end == std::string::npos ? std::string() : text.substr(0, end + 1);
}

mitk::DiffusionDICOMFileReader::DiffusionDICOMFileReader() = default;

mitk::DiffusionDICOMFileReader::~DiffusionDICOMFileReader() = default;

bool mitk::DiffusionDICOMFileReader
::LoadImages()
{
  unsigned int numberOfOutputs = this->GetNumberOfOutputs();
  bool success = true;

  for(unsigned int o = 0; o < numberOfOutputs; ++o)
  {
    success &= this->LoadSingleOutputImage( this->m_OutputHeaderContainer.at(o),
                                            this->InternalGetOutput(o), this->m_IsMosaicData.at(o) );
  }

  return success;
}

bool mitk::DiffusionDICOMFileReader
::LoadSingleOutputImage( DiffusionHeaderDICOMFileReader::DICOMHeaderListType retrievedHeader,
                         DICOMImageBlockDescriptor& block, bool is_mosaic)
{
  // prepare data reading
  DiffusionDICOMFileReaderHelper helper;
  DiffusionDICOMFileReaderHelper::VolumeFileNamesContainer filenames;

  const DICOMImageFrameList& frames = block.GetImageFrameList();
  int numberOfDWImages = block.GetIntProperty("timesteps", 1);

  int numberOfFramesPerDWImage = frames.size() / numberOfDWImages;

  assert( int( double((double) frames.size() / (double) numberOfDWImages)) == numberOfFramesPerDWImage );
  for( int idx = 0; idx < numberOfDWImages; idx++ )
  {
    std::vector< std::string > FileNamesPerVolume;

    auto timeStepStart = frames.begin() + idx * numberOfFramesPerDWImage;
    auto timeStepEnd   = frames.begin() + (idx+1) * numberOfFramesPerDWImage;
    for (auto frameIter = timeStepStart;
         frameIter != timeStepEnd;
         ++frameIter)
    {
      FileNamesPerVolume.push_back( (*frameIter)->Filename );
    }

    filenames.push_back( FileNamesPerVolume );
  }

  // TODO : only prototyping to test loading of diffusion images
  // we need some solution for the different types
  mitk::Image::Pointer output_image = mitk::Image::New();

  mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer directions =
      mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();


  double max_bvalue = 0;
  for( int idx = 0; idx < numberOfDWImages; idx++ )
  {
    DiffusionImageDICOMHeaderInformation header = retrievedHeader.at(idx);

    if( max_bvalue < header.b_value )
      max_bvalue = header.b_value;
  }

  // normalize the retrieved gradient directions according to the set b-value (maximal one)
  for( int idx = 0; idx < numberOfDWImages; idx++ )
  {
    DiffusionImageDICOMHeaderInformation header = retrievedHeader.at(idx);
    mitk::DiffusionPropertyHelper::GradientDirectionType grad = header.g_vector;

    grad.normalize();
    // A series of b=0 images only has no reference b-value.
    grad *= max_bvalue > 0 ? sqrt( header.b_value / max_bvalue ) : 0.0;

    directions->push_back( grad );
  }

  if( is_mosaic && this->m_ResolveMosaic )
  {
    mitk::DiffusionHeaderSiemensMosaicDICOMFileReader::Pointer mosaic_reader =
        mitk::DiffusionHeaderSiemensMosaicDICOMFileReader::New();

    // retrieve the remaining meta-information needed for mosaic reconstruction
    // it suffices to get it exemplatory from the first file in the file list
    mosaic_reader->RetrieveMosaicInformation( filenames.at(0).at(0) );

    mitk::MosaicDescriptor mdesc = mosaic_reader->GetMosaicDescriptor();
    if (mdesc.nimages == 0)
      mitkThrow() << "No NumberOfImagesInMosaic in " << filenames.at(0).at(0);
    mitk::CastToMitkImage( helper.LoadMosaicToVector<short, 3>( filenames, mdesc ), output_image );

  }
  else
  {
    mitk::CastToMitkImage( helper.LoadToVector<short, 3>( filenames ), output_image );
  }

  // DICOM gradients are the header vectors, rotated by the image matrix if
  // m_ApplyRotationToGradients, without measurement frame; DICOM is not part of the world-vector
  // convention of NRRD and NIfTI, because there is no DICOM test data to check it against.
  mitk::DiffusionPropertyHelper::MeasurementFrameType identity;
  identity.set_identity();
  mitk::DiffusionPropertyHelper::InitializeImage(output_image, directions.GetPointer(), max_bvalue, identity, m_ApplyRotationToGradients);
  output_image->SetProperty("diffusion.dicom.importname", mitk::StringProperty::New( helper.GetOutputName(filenames) ) );

  block.SetMitkImage( (mitk::Image::Pointer) output_image );

  return block.GetMitkImage().IsNotNull();

}

void mitk::DiffusionDICOMFileReader
::AnalyzeInputFiles()
{
  m_Study_names.clear();
  m_Series_names.clear();
  m_OutputHeaderContainer.clear();
  m_OutputReaderContainer.clear();
  m_IsMosaicData.clear();
  m_sop_instance_uids.clear();
  m_frame_of_reference_uids.clear();
  m_series_instance_uids.clear();
  m_study_instance_uids.clear();
  m_patient_names.clear();
  m_patient_ids.clear();
  this->SetGroup3DandT(true);

  Superclass::AnalyzeInputFiles();

  // collect output from superclass
  size_t number_of_outputs = this->GetNumberOfOutputs();

  if(number_of_outputs == 0)
  {
    MITK_ERROR << "Failed to parse input, retrieved 0 outputs from SeriesGDCMReader ";
  }

  MITK_INFO("diffusion.dicomreader") << "Retrieved " << number_of_outputs << " outputs.";
  std::vector< bool > valid_input_list;

  for( unsigned int outputidx = 0; outputidx < this->GetNumberOfOutputs(); outputidx++ )
  {
    DICOMImageBlockDescriptor block_0 = this->GetOutput(outputidx);

    // collect vendor ID from the first output, first image
    StringList inputFilename;
    DICOMImageFrameInfo::Pointer frame_0 = block_0.GetImageFrameList().at(0);
    inputFilename.push_back( frame_0->Filename );


    mitk::DiffusionHeaderDICOMFileReader::Pointer headerReader;

    bool isMosaic = false;
    gdcm::Scanner gdcmScanner;

    gdcm::Tag t_sop_instance_uid(0x0008, 0x0018);
    gdcm::Tag t_frame_of_reference_uid(0x0020, 0x0052);
    gdcm::Tag t_series_instance_uid(0x0020, 0x000E);
    gdcm::Tag t_study_instance_uid(0x0020, 0x000D);
    gdcm::Tag t_patient_name(0x0010, 0x0010);
    gdcm::Tag t_patient_id(0x0010, 0x0020);

    gdcm::Tag t_vendor(0x008, 0x0070);
    gdcm::Tag t_imagetype(0x0008, 0x0008);
    gdcm::Tag t_StudyDescription(0x0008, 0x1030);
    gdcm::Tag t_SeriesDescription(0x0008, 0x103E);

    // add DICOM Tag for vendor
    gdcmScanner.AddTag( t_vendor );
    // add DICOM Tag for image type
    gdcmScanner.AddTag( t_imagetype );
    gdcmScanner.AddTag( t_StudyDescription );
    gdcmScanner.AddTag( t_SeriesDescription );

    gdcmScanner.AddTag( t_sop_instance_uid );
    gdcmScanner.AddTag( t_frame_of_reference_uid );
    gdcmScanner.AddTag( t_series_instance_uid );
    gdcmScanner.AddTag( t_study_instance_uid );
    gdcmScanner.AddTag( t_patient_name );
    gdcmScanner.AddTag( t_patient_id );

    // Per-output DICOM values; stored only once the output turns out to be valid, so that they
    // stay in step with the outputs.
    auto tagValue = [&](const gdcm::Tag& tag) {
      const char* value = gdcmScanner.GetValue( frame_0->Filename.c_str(), tag );
      const std::string text = value != nullptr ? TrimDicomPadding(value) : std::string();
      return text.empty() ? std::string("-") : text;
    };
    std::string sop_instance_uid, frame_of_reference_uid, series_instance_uid, study_instance_uid, patient_name, patient_id;

    if( gdcmScanner.Scan( inputFilename ) )
    {
      sop_instance_uid = tagValue(t_sop_instance_uid);
      frame_of_reference_uid = tagValue(t_frame_of_reference_uid);
      series_instance_uid = tagValue(t_series_instance_uid);
      study_instance_uid = tagValue(t_study_instance_uid);
      patient_name = tagValue(t_patient_name);
      patient_id = tagValue(t_patient_id);

      // retrieve both vendor and image type
      const char* ch_vendor = gdcmScanner.GetValue( frame_0->Filename.c_str(), t_vendor );
      const char* ch_image_type = gdcmScanner.GetValue( frame_0->Filename.c_str(), t_imagetype );

      if( ch_vendor == nullptr || ch_image_type == nullptr )
      {
        MITK_WARN << "Unable to retrieve vendor/image information from " << frame_0->Filename.c_str() << "\n" <<
                     "Output " << outputidx+1 << " is not valid, skipping analysis.";
        valid_input_list.push_back(false);
        continue;

      }

      std::string vendor = std::string( ch_vendor );
      std::string image_type = std::string( ch_image_type );
      MITK_INFO("diffusion.dicomreader") << "Output " << outputidx+1 << "  Got vendor: " << vendor << " image type " << image_type;

      headerReader = CreateHeaderReader(vendor, image_type, isMosaic);
      if( headerReader.IsNull() )
      {
        MITK_ERROR << "No header reader for given vendor. ";
        valid_input_list.push_back(false);
        continue;
      }

    }
    else
    {
      valid_input_list.push_back(false);
      continue;
    }

    // iterate over the threeD+t block; every volume needs its diffusion header
    bool canread = true;

    int numberOfTimesteps = block_0.GetIntProperty("timesteps", 1);
    int framesPerTimestep = block_0.GetImageFrameList().size() / numberOfTimesteps;

    for( int idx = 0; idx < numberOfTimesteps && canread; idx++ )
    {
      int access_idx = idx * framesPerTimestep;
      DICOMImageFrameInfo::Pointer frame = this->GetOutput( outputidx ).GetImageFrameList().at( access_idx );
      canread = headerReader->ReadDiffusionHeader( frame->Filename );
    }

    if( canread )
    {
      // collect the information
      mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType retrievedHeader = headerReader->GetHeaderInformation();

      m_IsMosaicData.push_back(isMosaic);
      m_OutputHeaderContainer.push_back( retrievedHeader );
      m_OutputReaderContainer.push_back( headerReader );
      valid_input_list.push_back(true);

      m_sop_instance_uids.push_back(sop_instance_uid);
      m_frame_of_reference_uids.push_back(frame_of_reference_uid);
      m_series_instance_uids.push_back(series_instance_uid);
      m_study_instance_uids.push_back(study_instance_uid);
      m_patient_names.push_back(patient_name);
      m_patient_ids.push_back(patient_id);

      m_Study_names.push_back(tagValue(t_StudyDescription));
      m_Series_names.push_back(tagValue(t_SeriesDescription));
    }
    else
    {
      valid_input_list.push_back(false);
    }
  }

  // check status of the outputs and remove the non-valid
  std::vector< DICOMImageBlockDescriptor > valid_outputs;
  for ( unsigned int outputidx = 0; outputidx < this->GetNumberOfOutputs(); ++outputidx )
  {
    if( valid_input_list.at(outputidx) )
    {
      valid_outputs.push_back( this->InternalGetOutput( outputidx ) );
    }
  }

  // clear complete list
  this->ClearOutputs();
  this->SetNumberOfOutputs( valid_outputs.size() );

  // insert only the valid ones
  for ( unsigned int outputidx = 0; outputidx < valid_outputs.size(); ++outputidx )
    this->SetOutput( outputidx, valid_outputs.at( outputidx ) );

  for( unsigned int outputidx = 0; outputidx < this->GetNumberOfOutputs(); outputidx++ )
  {
    // TODO : Analyze outputs + header information, i.e. for the loading confidence
    MITK_INFO("diffusion.dicomreader") << "---- DICOM Analysis Report ----  :: Output " << outputidx+1 << " of " << this->GetNumberOfOutputs();

    try{
      PerformHeaderAnalysis( this->m_OutputHeaderContainer.at( outputidx) );
    }
    catch( const std::exception& se)
    {
      MITK_ERROR << "STD Exception " << se.what();
    }

    MITK_INFO("diffusion.dicomreader") << "===========================================";
  }
}


mitk::DiffusionHeaderDICOMFileReader::Pointer mitk::DiffusionDICOMFileReader
::CreateHeaderReader(const std::string& vendor, const std::string& imageType, bool& isMosaic)
{
  isMosaic = false;
  auto contains = [&vendor](const char* text) { return vendor.find(text) != std::string::npos; };

  if (contains("SIEMENS") || contains("Siemens HealthCare GmbH") || contains("Siemens Healthineers"))
  {
    if (imageType.find("MOSAIC") != std::string::npos)
    {
      isMosaic = true;
      return mitk::DiffusionHeaderSiemensMosaicDICOMFileReader::New().GetPointer();
    }
    return mitk::DiffusionHeaderSiemensDICOMFileReader::New().GetPointer();
  }
  if (contains("GE"))
    return mitk::DiffusionHeaderGEDICOMFileReader::New().GetPointer();
  if (contains("Philips"))
    return mitk::DiffusionHeaderPhilipsDICOMFileReader::New().GetPointer();
  return nullptr;
}

std::vector<mitk::Image::Pointer> mitk::DiffusionDICOMFileReader
::LoadDiffusionDicomSeries(const StringList& files, bool resolveMosaic, bool applyRotationToGradients, const std::string& fallbackName)
{
  // The header readers parse decimal numbers with the C library.
  mitk::LocaleSwitch localeSwitch("C");

  // In addition to the series split of ClassicDICOMSeriesReader, split by geometry, so that every
  // output is one 3D+t block.
  mitk::DICOMTagBasedSorter::Pointer tagSorter = mitk::DICOMTagBasedSorter::New();
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0028, 0x0010) ); // Number of Rows
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0028, 0x0011) ); // Number of Columns
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0028, 0x0030) ); // Pixel Spacing
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0018, 0x1164) ); // Imager Pixel Spacing
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0020, 0x0037) ); // Image Orientation (Patient)
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0020, 0x000e) ); // Series Instance UID
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0018, 0x0050) ); // Slice Thickness
  tagSorter->AddDistinguishingTag( mitk::DICOMTag(0x0028, 0x0008) ); // Number of Frames

  auto tag_sop_instance_uid = mitk::DICOMSortByTag::New( mitk::DICOMTag(0x0008, 0x0018), nullptr ); // SOP instance UID (last resort, not really meaningful but decides clearly)
  auto tag_trigger_time = mitk::DICOMSortByTag::New( mitk::DICOMTag(0x0018, 0x1060), tag_sop_instance_uid.GetPointer() ); // trigger time
  auto tag_aqcuisition_time = mitk::DICOMSortByTag::New( mitk::DICOMTag(0x0008, 0x0032), tag_trigger_time.GetPointer()); // aqcuisition time
  auto tag_aqcuisition_number = mitk::DICOMSortByTag::New( mitk::DICOMTag(0x0020, 0x0012), tag_aqcuisition_time.GetPointer()); // aqcuisition number
  mitk::DICOMSortCriterion::ConstPointer sorting = mitk::SortByImagePositionPatient::New(tag_aqcuisition_number.GetPointer()).GetPointer();
  tagSorter->SetSortCriterion( sorting );

  Pointer reader = Self::New();
  reader->SetResolveMosaic( resolveMosaic );
  reader->SetApplyRotationToGradients( applyRotationToGradients );
  reader->AddSortingElement( tagSorter );
  reader->SetInputFiles( files );
  try
  {
    reader->AnalyzeInputFiles();
  }
  catch( const itk::ExceptionObject &e)
  {
    MITK_ERROR << "Failed to analyze data. " << e.what();
  }
  catch( const std::exception &se)
  {
    MITK_ERROR << "Std Exception " << se.what();
  }

  reader->LoadImages();

  std::vector<Image::Pointer> images;
  for( unsigned int o = 0; o < reader->GetNumberOfOutputs(); o++ )
  {
    Image::Pointer image = reader->GetOutput(o).GetMitkImage();

    std::string name = reader->GetSeriesName(o);
    if (name == "-")
      name = reader->GetStudyName(o);
    if (name == "-")
    {
      name = fallbackName;
      if (name.empty())
        image->GetPropertyList()->GetStringProperty("diffusion.dicom.importname", name);
    }
    image->SetProperty("name", StringProperty::New(name));

    const std::pair<const char*, const std::vector<std::string>*> dicomValues[] = {
      {"DICOM.patient_id", &reader->m_patient_ids},
      {"DICOM.patient_name", &reader->m_patient_names},
      {"DICOM.study_instance_uid", &reader->m_study_instance_uids},
      {"DICOM.series_instance_uid", &reader->m_series_instance_uids},
      {"DICOM.sop_instance_uid", &reader->m_sop_instance_uids},
      {"DICOM.frame_of_reference_uid", &reader->m_frame_of_reference_uids}};
    for (const auto& [key, values] : dicomValues)
    {
      if (values->size() > o)
        image->GetPropertyList()->SetStringProperty(key, values->at(o).c_str());
    }

    images.push_back(image);
  }
  return images;
}

bool mitk::DiffusionDICOMFileReader
::CanHandleFile(const std::string & /* filename */)
{
  //FIXME :
  return true;
}

