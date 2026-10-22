/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageNrrdReader_cpp
#define __mitkDiffusionImageNrrdReader_cpp

#include "mitkDiffusionImageNrrdReader.h"

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
#include "itkNrrdImageIO.h"

#include "mitkCustomMimeType.h"
#include "mitkDiffusionImageMimeTypes.h"

#include <mitkITKImageImport.h>
#include <mitkImageWriteAccessor.h>
#include <mitkImageDataItem.h>
#include "mitkIOUtil.h"
#include <mitkLocaleSwitch.h>

#include <charconv>
#include <optional>

namespace
{
  using GradientDirectionType = mitk::DiffusionPropertyHelper::GradientDirectionType;
  using MeasurementFrameType = mitk::DiffusionPropertyHelper::MeasurementFrameType;
  using GradientDirectionContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

  /** Frames whose entries are all zero mean "no frame" in old MITK files. */
  bool IsZero(const MeasurementFrameType& frame)
  {
    return frame.absolute_value_max() <= 1e-9;
  }

  /**
   * Reads gradients, reference b-value and measurement frame from the NRRD key/value pairs.
   * Throws if they do not describe a DWI with one gradient per image component. Returns whether the
   * file was written by an earlier MITK version.
   */
  bool ParseDiffusionHeader(const itk::MetaDataDictionary& dictionary,
                            unsigned int numberOfComponents,
                            GradientDirectionContainerType* gradients,
                            MeasurementFrameType& measurementFrame,
                            float& referenceBValue)
  {
    const std::string gradientPrefix = "DWMRI_gradient_";
    std::map<unsigned int, GradientDirectionType> gradientsByIndex;
    std::optional<MeasurementFrameType> standardFrame;
    std::optional<MeasurementFrameType> mitkFrame;
    bool hasBValue = false;

    for (const auto& key : dictionary.GetKeys())
    {
      // ITK stores the standard NRRD field "measurement frame" (written by Slicer, DTIPrep, ITK) as
      // a matrix whose rows are the frame vectors, not as a string.
      if (key == "NRRD_measurement frame")
      {
        std::vector<std::vector<double>> rows;
        if (itk::ExposeMetaData(dictionary, key, rows) && rows.size() == 3 && rows[0].size() == 3 && rows[1].size() == 3 && rows[2].size() == 3)
        {
          MeasurementFrameType frame;
          for (unsigned int r = 0; r < 3; ++r)
            for (unsigned int c = 0; c < 3; ++c)
              frame(r, c) = rows[r][c];
          standardFrame = frame;
        }
        continue;
      }

      std::string value;
      if (!itk::ExposeMetaData<std::string>(dictionary, key, value))
        continue;

      if (key.starts_with(gradientPrefix))
      {
        const std::string indexString = key.substr(gradientPrefix.size());
        unsigned int index = 0;
        const auto [end, error] = std::from_chars(indexString.data(), indexString.data() + indexString.size(), index);
        if (error != std::errc() || end != indexString.data() + indexString.size())
          mitkThrow() << "Invalid gradient key " << key;

        double x, y, z;
        if (sscanf(value.c_str(), "%lf %lf %lf", &x, &y, &z) != 3)
          mitkThrow() << "Invalid gradient " << key << ":=" << value;
        GradientDirectionType g;
        g[0] = x;
        g[1] = y;
        g[2] = z;
        gradientsByIndex[index] = g;
      }
      else if (key == "DWMRI_b-value")
      {
        referenceBValue = atof(value.c_str());
        hasBValue = true;
      }
      else if (key == "DWMRI_NEX")
      {
        mitkThrow() << "DWMRI_NEX keys (repeated gradients) are not supported";
      }
      else if (key == "measurement frame")
      {
        // MITK's own key/value pair, written by the DWI NRRD writer.
        MeasurementFrameType frame;
        if (sscanf(value.c_str(), " ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ",
                   &frame(0, 0), &frame(0, 1), &frame(0, 2), &frame(1, 0), &frame(1, 1), &frame(1, 2), &frame(2, 0), &frame(2, 1), &frame(2, 2)) != 9)
          mitkThrow() << "Invalid measurement frame " << value;
        mitkFrame = frame;
      }
    }

    if (!hasBValue)
      mitkThrow() << "DWMRI_b-value is missing";
    if (gradientsByIndex.size() != numberOfComponents || (!gradientsByIndex.empty() && gradientsByIndex.rbegin()->first != numberOfComponents - 1))
      mitkThrow() << "Found " << gradientsByIndex.size() << " gradients DWMRI_gradient_0000.. for " << numberOfComponents << " image components";
    for (const auto& [index, g] : gradientsByIndex)
      gradients->InsertElement(index, g);

    // Earlier MITK versions wrote only their own key. Current MITK writes both, with identical values.
    const bool legacyMitkFile = mitkFrame && !standardFrame;
    const auto& frame = legacyMitkFile ? mitkFrame : standardFrame;
    if (frame && !IsZero(*frame))
      measurementFrame = *frame;
    return legacyMitkFile;
  }
}

namespace mitk
{

  DiffusionImageNrrdReader::
  DiffusionImageNrrdReader(const DiffusionImageNrrdReader & other)
    : AbstractFileReader(other)
  {
  }


  DiffusionImageNrrdReader* DiffusionImageNrrdReader::Clone() const
  {
    return new DiffusionImageNrrdReader(*this);
  }


  DiffusionImageNrrdReader::
  ~DiffusionImageNrrdReader()
  {}

  DiffusionImageNrrdReader::
  DiffusionImageNrrdReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_NAME() ), mitk::DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_DESCRIPTION() )
  {
    Options defaultOptions;
    defaultOptions["Apply image rotation to gradients"] = true;
    this->SetDefaultOptions(defaultOptions);
    m_ServiceReg = this->RegisterService();
  }

  std::vector<itk::SmartPointer<mitk::BaseData> >
  DiffusionImageNrrdReader::
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


  void DiffusionImageNrrdReader::InternalRead()
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

        MITK_INFO << "DiffusionImageNrrdReader: reading image information";
        VectorImageType::Pointer itkVectorImage;
        {
          typedef itk::ImageFileReader<VectorImageType> FileReaderType;
          FileReaderType::Pointer reader = FileReaderType::New();
          reader->SetFileName(this->GetInputLocation());
          itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
          reader->SetImageIO(io);
          reader->Update();
          itkVectorImage = reader->GetOutput();
        }

        GradientDirectionContainerType::Pointer DiffusionVectors = GradientDirectionContainerType::New();
        MeasurementFrameType MeasurementFrame;
        MeasurementFrame.set_identity();
        float BValue = -1;
        const bool legacyMitkFile = ParseDiffusionHeader(itkVectorImage->GetMetaDataDictionary(), itkVectorImage->GetVectorLength(), DiffusionVectors, MeasurementFrame, BValue);

        outputForCache = mitk::GrabItkImageMemory( itkVectorImage);

        // create BValueMap
        mitk::BValueMapProperty::BValueMap BValueMap = mitk::BValueMapProperty::CreateBValueMap(DiffusionVectors,BValue);
        mitk::DiffusionPropertyHelper::SetOriginalGradientContainer(outputForCache, DiffusionVectors);
        mitk::DiffusionPropertyHelper::SetMeasurementFrame(outputForCache, MeasurementFrame);
        mitk::DiffusionPropertyHelper::SetBValueMap(outputForCache, BValueMap);
        mitk::DiffusionPropertyHelper::SetReferenceBValue(outputForCache, BValue);
        // NRRD defines the gradients through the measurement frame alone, which makes the modified
        // gradients world vectors. Files written by earlier MITK versions are read as they always
        // were, with the image rotation applied before the frame.
        const bool applyRotation = legacyMitkFile && us::any_cast<bool>(this->GetOptions()["Apply image rotation to gradients"]);
        mitk::DiffusionPropertyHelper::SetApplyMatrixToGradients(outputForCache, applyRotation);
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
