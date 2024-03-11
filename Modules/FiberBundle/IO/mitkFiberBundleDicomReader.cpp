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
#include <itkMetaDataObject.h>
#include <vtkPolyData.h>
#include <vtkDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <itksys/SystemTools.hxx>
#include <vtkCleanPolyData.h>
#include <mitkCustomMimeType.h>
#include "mitkFiberBundleMimeTypes.h"
#include <vtkTransformPolyDataFilter.h>
#include <mitkLexicalCast.h>

#include "dcmtk/ofstd/ofcond.h"
#include "dcmtk/dcmtract/trctractographyresults.h"
#include "dcmtk/dcmtract/trctrack.h"

mitk::FiberBundleDicomReader::FiberBundleDicomReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME(), "DICOM Fiber Bundle Reader" )
{
  m_ServiceReg = this->RegisterService();
}

mitk::FiberBundleDicomReader::FiberBundleDicomReader(const FiberBundleDicomReader &other)
  :mitk::AbstractFileReader(other)
{
}

mitk::FiberBundleDicomReader * mitk::FiberBundleDicomReader::Clone() const
{
  return new FiberBundleDicomReader(*this);
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleDicomReader::DoRead()
{

  std::vector<itk::SmartPointer<mitk::BaseData> > output_fibs;
  try
  {
    const std::string& locale = "C";
    const std::string& currLocale = setlocale( LC_ALL, nullptr );
    setlocale(LC_ALL, locale.c_str());

    std::string filename = this->GetInputLocation();


    OFCondition result;
    TrcTractographyResults *trc = nullptr;
    result = TrcTractographyResults::loadFile(filename.c_str(), trc);
    if (result.bad())
      mitkThrow() << "Unable to load tractography dicom file: " << result.text();

    OFString val = "-";
    trc->getPatient().getPatientName(val);
    MITK_INFO << "Patient Name: " <<  val; val = "-";
    trc->getStudy().getStudyInstanceUID(val);
    MITK_INFO << "Study       : " << val; val = "-";
    trc->getSeries().getSeriesInstanceUID(val);
    MITK_INFO << "Series      : " << val; val = "-";
    trc->getSOPCommon().getSOPInstanceUID(val);
    MITK_INFO << "Instance    : " << val; val = "-";
    MITK_INFO << "-------------------------------------------------------------------------";
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

      CodeSequenceMacro algoCode = sets[ts]->getTrackingAlgorithmIdentification().at(0)->getAlgorithmFamilyCode();
      val = "0"; algoCode.getCodeValue(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.algo_family_code.value",val.c_str());
      val = "0"; algoCode.getCodeMeaning(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.algo_family_code.meaning",val.c_str());

      CodeSequenceMacro modelCode = sets[ts]->getDiffusionModelCode();
      val = "0"; modelCode.getCodeValue(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.model_code.value",val.c_str());
      val = "0"; modelCode.getCodeMeaning(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.model_code.meaning",val.c_str());

      CodeWithModifiers anatomy = sets[ts]->getTrackSetAnatomy();
      val = "0"; anatomy.getCodeValue(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.anatomy.value",val.c_str());
      val = "0"; anatomy.getCodeMeaning(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.anatomy.meaning",val.c_str());

      val = "0"; trc->getPatient().getPatientID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.patient_id",val.c_str());
      val = "0"; trc->getPatient().getPatientName(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.patient_name",val.c_str());
      val = "0"; trc->getStudy().getStudyInstanceUID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.study_instance_uid",val.c_str());
      val = "0"; trc->getSeries().getSeriesInstanceUID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.series_instance_uid",val.c_str());
      val = "0"; trc->getSOPCommon().getSOPInstanceUID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.sop_instance_uid",val.c_str());
      val = "0"; trc->getSOPCommon().getSOPClassUID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.sop_class_uid",val.c_str());
      val = "0"; trc->getFrameOfReference().getFrameOfReferenceUID(val);
      fib->GetPropertyList()->SetStringProperty("DICOM.frame_of_reference_uid",val.c_str());

      output_fibs.push_back(fib.GetPointer());
      MITK_INFO << "Fiber bundle read";
    }
    delete trc;

    setlocale(LC_ALL, currLocale.c_str());
    return output_fibs;
  }
  catch(...)
  {
    throw;
  }
  return output_fibs;
}
