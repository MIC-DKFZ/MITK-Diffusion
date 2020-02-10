/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleDicomWriter.h"
#include <vtkSmartPointer.h>
#include <itksys/SystemTools.hxx>
#include <itkSize.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <mitkAbstractFileWriter.h>
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"
#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */
#include "dcmtk/dcmtract/trctractographyresults.h"

mitk::FiberBundleDicomWriter::FiberBundleDicomWriter()
  : mitk::AbstractFileWriter(mitk::FiberBundle::GetStaticNameOfClass(), mitk::DiffusionIOMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME(), "DICOM Fiber Bundle Writer")
{
  RegisterService();
}

mitk::FiberBundleDicomWriter::FiberBundleDicomWriter(const mitk::FiberBundleDicomWriter & other)
  :mitk::AbstractFileWriter(other)
{}

mitk::FiberBundleDicomWriter::~FiberBundleDicomWriter()
{}

mitk::FiberBundleDicomWriter * mitk::FiberBundleDicomWriter::Clone() const
{
  return new mitk::FiberBundleDicomWriter(*this);
}

void mitk::FiberBundleDicomWriter::Write()
{
  try
  {
    mitk::FiberBundle::ConstPointer fib = dynamic_cast<const mitk::FiberBundle*>(this->GetInput());
    vtkPolyData* poly = fib->GetFiberPolyData();
    mitk::PropertyList* p_list = fib->GetPropertyList();

    std::string patient_id = "";
    if (!p_list->GetStringProperty("DICOM.patient_id", patient_id))
      patient_id = "0";
    std::string patient_name = "";
    if (!p_list->GetStringProperty("DICOM.patient_name", patient_name))
      patient_name = "0";
    std::string study_instance_uid = "";
    if (!p_list->GetStringProperty("DICOM.study_instance_uid", study_instance_uid))
      study_instance_uid = "0";
    std::string series_instance_uid = "";
    if (!p_list->GetStringProperty("DICOM.series_instance_uid", series_instance_uid))
      series_instance_uid = "0";
    std::string sop_instance_uid = "";
    if (!p_list->GetStringProperty("DICOM.sop_instance_uid", sop_instance_uid))
      sop_instance_uid = "0";
    std::string sop_class_uid = "";
    if (!p_list->GetStringProperty("DICOM.sop_class_uid", sop_instance_uid))
      sop_class_uid = "0";
    std::string frame_of_reference_uid = "";
    if (!p_list->GetStringProperty("DICOM.frame_of_reference_uid", frame_of_reference_uid))
      frame_of_reference_uid = "0";
    std::string algo_code_value = "";
    if (!p_list->GetStringProperty("DICOM.algo_family_code.value", algo_code_value))
      algo_code_value = "0";
    std::string algo_code_meaning = "";
    if (!p_list->GetStringProperty("DICOM.algo_family_code.meaning", algo_code_meaning))
      algo_code_meaning = "0";
    std::string model_code_value = "";
    if (!p_list->GetStringProperty("DICOM.model_code.value", model_code_value))
      model_code_value = "0";
    std::string model_code_meaning = "";
    if (!p_list->GetStringProperty("DICOM.model_code.meaning", model_code_meaning))
      model_code_meaning = "0";
    std::string anatomy_value = "";
    if (!p_list->GetStringProperty("DICOM.anatomy.value", anatomy_value))
      anatomy_value = "0";
    std::string anatomy_meaning = "";
    if (!p_list->GetStringProperty("DICOM.anatomy.meaning", anatomy_meaning))
      anatomy_meaning = "0";

    const std::string& locale = "C";
    const std::string& currLocale = setlocale( LC_ALL, nullptr );
    setlocale(LC_ALL, locale.c_str());

    // Instance Number, Label, Description, Creator's Name
    ContentIdentificationMacro id("1", "TRACTOGRAM", "Tractogram processed with MITK Diffusion", "MIC@DKFZ");

    // Manufacturer, model name, serial number, software version(s)
    IODEnhGeneralEquipmentModule::EquipmentInfo equipment("MIC@DKFZ", "dcmtract library", "0815", OFFIS_DCMTK_VERSION_STRING);
    IODReferences refs;

    // We need at least one image reference this Tractography Results object is based on.
    // We provide: Patient ID, Study Instance UID, Series Instance UID, SOP Instance UID, SOP Class UID
    IODImageReference* ref = new IODImageReference(patient_id.c_str(), study_instance_uid.c_str(), series_instance_uid.c_str(), sop_instance_uid.c_str(), UID_MRImageStorage);
    refs.add(ref);

    std::time_t t = std::time(nullptr);
    char date_buffer[20];
    std::strftime(date_buffer, sizeof(date_buffer), "%Y%m%d", std::gmtime(&t));
    char time_buffer[20];
    std::strftime(time_buffer, sizeof(time_buffer), "%H%M%S", std::gmtime(&t));
    OFString contentDate(date_buffer);
    OFString contentTime(time_buffer);

    OFString val = "-";
    TrcTractographyResults *trc = nullptr;
    TrcTractographyResults::create(id, contentDate, contentTime, equipment, refs, trc);
    trc->getStudy().setStudyInstanceUID(study_instance_uid.c_str());
    trc->getSeries().setSeriesInstanceUID(series_instance_uid.c_str());
    trc->getSOPCommon().setSOPInstanceUID(sop_instance_uid.c_str());
    trc->getSeries().getSeriesInstanceUID(val);

    trc->getPatient().setPatientID(patient_id.c_str());
    trc->getPatient().setPatientName(patient_name.c_str());
    trc->getSeries().setSeriesDescription("Tractogram processed with MITK Diffusion");

    // Frame of Reference is required; could be the same as from related MR series
    trc->getFrameOfReference().setFrameOfReferenceUID(frame_of_reference_uid.c_str());

    DcmItem item;
    item.putAndInsertOFStringArray(DCM_PatientID, patient_id);
    item.putAndInsertOFStringArray(DCM_StudyInstanceUID, study_instance_uid);
    item.putAndInsertOFStringArray(DCM_SeriesInstanceUID, series_instance_uid);
    item.putAndInsertOFStringArray(DCM_SOPClassUID, sop_class_uid);
    item.putAndInsertOFStringArray(DCM_SOPInstanceUID, sop_instance_uid);
    IODReference* ref2 = new IODReference();
    ref2->readFromItem(item);
    trc->getReferencedInstances().add(ref2);

    // GET STUFF FROM DICOM FILE (disbled since file probably not available)
//    trc->importHierarchy("/media/neher/Rumpelkammer/Demo/Tumor/KOPF_UNTERSUCHUNG_20170419_164344_842000/T1_FL2D_TRA_5MM_0013/WISSING_MARISKA.MR.KOPF_UNTERSUCHUNG.0013.0022.2017.07.31.14.09.57.918625.254331277.IMA",
//                         true,
//                         true,
//                         true,
//                         true);

//    trc->getReferencedInstances().addFromFiles({"/media/neher/Rumpelkammer/Demo/Tumor/KOPF_UNTERSUCHUNG_20170419_164344_842000/T1_FL2D_TRA_5MM_0013/WISSING_MARISKA.MR.KOPF_UNTERSUCHUNG.0013.0022.2017.07.31.14.09.57.918625.254331277.IMA"});

    // Create track set
    CodeWithModifiers anatomy("");
    anatomy.set(anatomy_value.c_str(), "SRT", anatomy_meaning.c_str());

    // Every CodeSequenceMacro has: Code Value, Coding Scheme Designator, Code Meaning
    CodeSequenceMacro diffusionModel(model_code_value.c_str(), "DCM", model_code_meaning.c_str());

    AlgorithmIdentificationMacro algorithmId;
    algorithmId.getAlgorithmFamilyCode().set(algo_code_value, "DCM", algo_code_meaning);
    algorithmId.setAlgorithmName("0");
    algorithmId.setAlgorithmVersion("0");
    algorithmId.setAlgorithmParameters("0");
    algorithmId.setAlgorithmSource("0");

    TrcTrackSet *set = nullptr;
    trc->addTrackSet("TRACTOGRAM", "Tractogram processed with MITK Diffusion", anatomy, diffusionModel, algorithmId, set);

    // Create trackset
    Uint16 cieLabColor[3]; // color the whole track with this color; we use some blue
    cieLabColor[0] = 30000; // L
    cieLabColor[1] = 0 ; // a
    cieLabColor[2] = 0 ; // b
    std::vector< Float32* > tracts;

    for (unsigned int i=0; i<fib->GetNumFibers(); i++)
    {
        vtkCell* cell = poly->GetCell(i);
        int numPoints = cell->GetNumberOfPoints();
        vtkPoints* points = cell->GetPoints();

        Float32* pointData = new Float32[numPoints*3];
        for(int i=0; i<numPoints ;i++)
        {
            double* p = points->GetPoint(i);

            pointData[i*3] = p[0];
            pointData[i*3+1] = p[1];
            pointData[i*3+2] = p[2];
        }
        tracts.push_back(pointData);

        TrcTrack* track = nullptr;
        set->addTrack(pointData, numPoints, cieLabColor, 1 /* numColors */, track);
    }

    // Save file
    OFCondition result = trc->saveFile(this->GetOutputLocation().c_str());
    delete trc;
    if (result.bad())
      mitkThrow() << "Unable to save tractography as DICOM file: " << result.text();

    for (Float32* tract : tracts)
      delete [] tract;

    setlocale(LC_ALL, currLocale.c_str());
    MITK_INFO << "DICOM Fiber bundle written to " << this->GetOutputLocation();
  }
  catch(...)
  {
    throw;
  }
}
