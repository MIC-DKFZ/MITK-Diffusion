/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleDicomReader.h"
#include "mitkFiberBundleMimeTypes.h"
#include "mitkFiberBundleDicomReferenceTags.h"

#include <itksys/SystemTools.hxx>
#include <mitkFiberBundle.h>
#include <mitkLexicalCast.h>
#include <mitkLocaleSwitch.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <memory>

#include "dcmtk/ofstd/ofcond.h"
#include "dcmtk/dcmtract/trctractographyresults.h"
#include "dcmtk/dcmtract/trctrack.h"

mitk::FiberBundleDicomReader::FiberBundleDicomReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME(), "DICOM Fiber Bundle Reader" )
{
  this->RegisterService();
}

mitk::FiberBundleDicomReader * mitk::FiberBundleDicomReader::Clone() const
{
  return new FiberBundleDicomReader(*this);
}

namespace
{
  /** Series instance UID of a referenced instance. The Referenced Instance Sequence of the
   *  tractography results holds only SOP class and instance UIDs; the series is listed in the
   *  Common Instance Reference module, in this study or in another one. */
  std::string FindReferencedSeriesUid(TrcTractographyResults* trc, const OFString& sopInstanceUid)
  {
    auto search = [&sopInstanceUid](OFVector<IODSeriesAndInstanceReferenceMacro::ReferencedSeriesItem*>& seriesItems)
    {
      for (auto* series : seriesItems)
      {
        for (auto* instance : series->getReferencedInstanceItems())
        {
          OFString uid;
          if (instance->getReferencedSOPInstanceUID(uid).good() && uid == sopInstanceUid)
          {
            OFString seriesUid;
            series->getSeriesInstanceUID(seriesUid);
            return std::string(seriesUid.c_str());
          }
        }
      }
      return std::string();
    };

    auto& module = trc->getCommonInstanceReference();
    std::string seriesUid = search(module.getReferencedSeriesItems());
    for (auto* study : module.getStudiesContainingOtherReferences())
    {
      if (!seriesUid.empty())
        break;
      seriesUid = search(study->getReferencedSeriesAndInstanceReferences().getReferencedSeriesItems());
    }
    return seriesUid;
  }
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleDicomReader::DoRead()
{

  std::vector<itk::SmartPointer<mitk::BaseData> > output_fibs;
  // Restored on every exit, also when DCMTK or the code below throws.
  mitk::LocaleSwitch localeSwitch("C");

  std::string filename = this->GetInputLocation();


  OFCondition result;
  TrcTractographyResults *loaded = nullptr;
  result = TrcTractographyResults::loadFile(filename.c_str(), loaded);
  // loadFile allocates the object before reading it, so it is owned also when reading failed.
  std::unique_ptr<TrcTractographyResults> trc(loaded);
  if (result.bad())
    mitkThrow() << "Unable to load tractography dicom file: " << result.text();

  // Patient name and UIDs are not logged: the log file must not identify the patient.
  OFString val;
  size_t numTrackSets = trc->getNumberOfTrackSets();
  OFVector<TrcTrackSet*>& sets = trc->getTrackSets();
  for (size_t ts = 0; ts < numTrackSets; ts++)
  {
    size_t numTracks = sets[ts]->getNumberOfTracks();
    MITK_INFO << "  Track Set #" << ts << ": " << numTracks << " Tracks, "
         << sets[ts]->getNumberOfTrackSetStatistics() << " Track Set Statistics, "
         << sets[ts]->getNumberOfTrackStatistics() << " Track Statistics, "
         << sets[ts]->getNumberOfMeasurements() << " Measurements ";

    vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

    for (size_t t = 0; t < numTracks; t++)
    {
      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      TrcTrack* track = sets[ts]->getTracks()[t];
      const Float32* vals = nullptr;
      size_t numPoints = track->getTrackData(vals);

      for (size_t v = 0; v < numPoints; ++v)
      {
        vtkIdType id = vtkNewPoints->InsertNextPoint(vals[v*3],vals[v*3+1],vals[v*3+2]);
        container->GetPointIds()->InsertNextId(id);
      }

      vtkNewCells->InsertNextCell(container);
    }

    vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
    fiberPolyData->SetPoints(vtkNewPoints);
    fiberPolyData->SetLines(vtkNewCells);
    FiberBundle::Pointer fib = FiberBundle::New(fiberPolyData);

    // Type 1 in the standard, but DCMTK reads a file without it (with a warning); the
    // properties then keep the placeholder "0", as for the other missing codes.
    const auto& algorithms = sets[ts]->getTrackingAlgorithmIdentification();
    CodeSequenceMacro* algoCode = nullptr;
    if (!algorithms.empty() && algorithms.front() != nullptr)
      algoCode = &algorithms.front()->getAlgorithmFamilyCode();
    val = "0"; if (algoCode != nullptr) algoCode->getCodeValue(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::AlgorithmFamilyCodeValue,val.c_str());
    val = "0"; if (algoCode != nullptr) algoCode->getCodeMeaning(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::AlgorithmFamilyCodeMeaning,val.c_str());

    CodeSequenceMacro modelCode = sets[ts]->getDiffusionModelCode();
    val = "0"; modelCode.getCodeValue(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::ModelCodeValue,val.c_str());
    val = "0"; modelCode.getCodeMeaning(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::ModelCodeMeaning,val.c_str());

    CodeWithModifiers anatomy = sets[ts]->getTrackSetAnatomy();
    val = "0"; anatomy.getCodeValue(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::AnatomyValue,val.c_str());
    val = "0"; anatomy.getCodeMeaning(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::AnatomyMeaning,val.c_str());

    val = "0"; trc->getPatient().getPatientID(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::PatientId,val.c_str());
    val = "0"; trc->getPatient().getPatientName(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::PatientName,val.c_str());
    val = "0"; trc->getStudy().getStudyInstanceUID(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::StudyInstanceUid,val.c_str());
    // DICOM.series_instance_uid, DICOM.sop_instance_uid and DICOM.sop_class_uid name the
    // referenced image, as the writer and the tag editor use them, so that saving the bundle
    // again refers to the same image; the tractogram's own UIDs are kept for display.
    const auto& references = trc->getReferencedInstances().get();
    if (!references.empty())
    {
      const OFString& referencedInstance = references.front()->m_SOPInstanceUID;
      fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::SeriesInstanceUid, FindReferencedSeriesUid(trc.get(), referencedInstance).c_str());
      fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::SopInstanceUid, referencedInstance.c_str());
      fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::SopClassUid, references.front()->m_SOPClassUID.c_str());
    }
    val = "0"; trc->getSeries().getSeriesInstanceUID(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::TractogramSeriesInstanceUid,val.c_str());
    val = "0"; trc->getSOPCommon().getSOPInstanceUID(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::TractogramSopInstanceUid,val.c_str());
    val = "0"; trc->getFrameOfReference().getFrameOfReferenceUID(val);
    fib->GetPropertyList()->SetStringProperty(FiberBundleDicomTags::FrameOfReferenceUid,val.c_str());

    output_fibs.push_back(fib.GetPointer());
    MITK_INFO << "Fiber bundle read";
  }
  return output_fibs;
  return output_fibs;
}
