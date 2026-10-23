/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkPeakImageWriter.h"

#include "mitkDiffusionIOMimeTypes.h"
#include <mitkImageToItk.h>
#include <mitkLocaleSwitch.h>
#include <mitkPeakImage.h>

#include <itkImageDuplicator.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMetaDataObject.h>
#include <itkNiftiImageIO.h>
#include <itkNrrdImageIO.h>
#include <itksys/SystemTools.hxx>

mitk::PeakImageWriter::PeakImageWriter()
  : AbstractFileWriter(PeakImage::GetStaticNameOfClass(),
                       CustomMimeType(DiffusionIOMimeTypes::PEAK_MIMETYPE_NAME()),
                       DiffusionIOMimeTypes::PEAK_MIMETYPE_DESCRIPTION())
{
  RegisterService();
}

mitk::PeakImageWriter::PeakImageWriter(const PeakImageWriter &other) : AbstractFileWriter(other)
{
}

mitk::PeakImageWriter::~PeakImageWriter() = default;

mitk::PeakImageWriter *mitk::PeakImageWriter::Clone() const
{
  return new PeakImageWriter(*this);
}

void mitk::PeakImageWriter::Write()
{
  const auto *input = dynamic_cast<const PeakImage *>(this->GetInput());
  if (input == nullptr)
    mitkThrow() << "Sorry, input to PeakImageWriter is not a peak image!";
  if (this->GetOutputLocation().empty())
    mitkThrow() << "Sorry, filename has not been set!";

  LocaleSwitch localeSwitch("C");
  auto caster = ImageToItk<PeakImage::ItkPeakImageType>::New();
  caster->SetInput(input);
  caster->Update();

  // A copy, so that the RAS conversion does not change the peaks in memory.
  auto duplicator = itk::ImageDuplicator<PeakImage::ItkPeakImageType>::New();
  duplicator->SetInputImage(caster->GetOutput());
  duplicator->Update();
  PeakImage::ItkPeakImageType::Pointer peaks = duplicator->GetOutput();

  auto writer = itk::ImageFileWriter<PeakImage::ItkPeakImageType>::New();
  const std::string ext = itksys::SystemTools::GetFilenameLastExtension(this->GetOutputLocation());
  if (ext == ".peak" || ext == ".nrrd")
  {
    auto io = itk::NrrdImageIO::New();
    io->UseCompressionOn();
    writer->SetImageIO(io);
  }
  else
  {
    // Peaks in NIfTI files are read as scanner RAS vectors by MRtrix and the tools that use its
    // peak files; MITK's world vectors are LPS. Every third value is z.
    itk::ImageRegionIteratorWithIndex<PeakImage::ItkPeakImageType> it(peaks, peaks->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      if (it.GetIndex()[3] % 3 != 2)
        it.Set(-it.Get());
    }
    itk::EncapsulateMetaData<std::string>(peaks->GetMetaDataDictionary(), "ITK_FileNotes", "MITK_peak_convention=RAS");
    writer->SetImageIO(itk::NiftiImageIO::New());
  }

  writer->SetInput(peaks);
  writer->SetFileName(this->GetOutputLocation());
  writer->SetUseCompression(true);
  writer->Update();
}
