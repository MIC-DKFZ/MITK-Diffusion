/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleVtkWriter.h"
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <itksys/SystemTools.hxx>
#include <itkSize.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <mitkAbstractFileWriter.h>
#include <mitkCustomMimeType.h>
#include <vtkXMLPolyDataWriter.h>
#include "mitkFiberBundleMimeTypes.h"
#include <vtkUnsignedCharArray.h>
#include <vtkTransformFilter.h>
#include <mitkGeometry3D.h>
#include <mitkLocaleSwitch.h>

mitk::FiberBundleVtkWriter::FiberBundleVtkWriter()
  : mitk::AbstractFileWriter(mitk::FiberBundle::GetStaticNameOfClass(), mitk::FiberBundleMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME(), "VTK Fiber Bundle Writer")
{
  Options defaultOptions;
  defaultOptions["Save as binary file"] = true;
  defaultOptions["Save as xml file (vtp style)"] = false;
  defaultOptions["Save color information"] = true;
  defaultOptions["Save fiber weights"] = true;
  this->SetDefaultOptions(defaultOptions);
  RegisterService();
}

mitk::FiberBundleVtkWriter::FiberBundleVtkWriter(const mitk::FiberBundleVtkWriter & other)
  :mitk::AbstractFileWriter(other)
{}

mitk::FiberBundleVtkWriter::~FiberBundleVtkWriter()
{}

mitk::FiberBundleVtkWriter * mitk::FiberBundleVtkWriter::Clone() const
{
  return new mitk::FiberBundleVtkWriter(*this);
}

void mitk::FiberBundleVtkWriter::Write()
{
  try
  {
    // ASCII files write numbers with the decimal point of the current locale
    mitk::LocaleSwitch localeSwitch("C");

    // A local file also for output streams; it is copied to the stream when it goes out of scope.
    LocalFile localFile(this);
    const std::string filename = localFile.GetFileName();

    mitk::FiberBundle::ConstPointer input = dynamic_cast<const mitk::FiberBundle*>(this->GetInput());
    if (input.IsNull())
      mitkThrow() << "The VTK fiber bundle writer needs a fiber bundle as input.";
    std::string ext = itksys::SystemTools::GetFilenameLastExtension(this->GetOutputLocation().c_str());
    Options options = this->GetOptions();

    vtkSmartPointer<vtkPolyData> fibPoly = input->GetFiberPolyData();

    bool ras = input->IsRAS();

    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    vtkSmartPointer< vtkMatrix4x4 > matrix = vtkSmartPointer< vtkMatrix4x4 >::New();
    matrix->Identity();
    if (ras)
    {
      matrix->SetElement(0,0,-matrix->GetElement(0,0));
      matrix->SetElement(0,1,-matrix->GetElement(0,1));
      matrix->SetElement(0,2,-matrix->GetElement(0,2));
      matrix->SetElement(0,3,-matrix->GetElement(0,3));

      matrix->SetElement(1,0,-matrix->GetElement(1,0));
      matrix->SetElement(1,1,-matrix->GetElement(1,1));
      matrix->SetElement(1,2,-matrix->GetElement(1,2));
      matrix->SetElement(1,3,-matrix->GetElement(1,3));

      geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

      vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
      transformFilter->SetInputData(fibPoly);
      transformFilter->SetTransform(geometry->GetVtkTransform());
      transformFilter->Update();
      fibPoly = transformFilter->GetPolyDataOutput();
    }

    // The bundle's polydata is shared with the mappers and must not change by saving, so the
    // weight and colour arrays are added to or removed from a shallow copy only.
    vtkSmartPointer<vtkPolyData> outputPoly = vtkSmartPointer<vtkPolyData>::New();
    outputPoly->ShallowCopy(fibPoly);
    fibPoly = outputPoly;


    if (us::any_cast<bool>(options["Save fiber weights"]))
    {
      MITK_INFO << "Adding fiber weight information";
      fibPoly->GetCellData()->AddArray(input->GetFiberWeights());
    }
    else if (fibPoly->GetCellData()->HasArray("FIBER_WEIGHTS"))
    {
      fibPoly->GetCellData()->RemoveArray("FIBER_WEIGHTS");
    }

    if (us::any_cast<bool>(options["Save color information"]))
    {
      MITK_INFO << "Adding color information";
      fibPoly->GetPointData()->AddArray(input->GetFiberColors());
    }
    else if (fibPoly->GetPointData()->HasArray("FIBER_COLORS"))
    {
      fibPoly->GetPointData()->RemoveArray("FIBER_COLORS");
    }

    // default extension is .fib
    if(ext == "")
    {
      ext = ".fib";
      this->SetOutputLocation(this->GetOutputLocation() + ext);
    }

    if (us::any_cast<bool>(options["Save as xml file (vtp style)"]))
    {
      vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetInputData(fibPoly);
      writer->SetFileName(filename.c_str());
      if (us::any_cast<bool>(options["Save as binary file"]))
      {
        MITK_INFO << "Writing fiber bundle as vtk binary file";
        writer->SetDataModeToBinary();
      }
      else
      {
        MITK_INFO << "Writing fiber bundle as vtk ascii file";
        writer->SetDataModeToAscii();
      }
      if (writer->Write() != 1 || writer->GetErrorCode() != 0)
        mitkThrow() << "Could not write fiber bundle to " << filename << " (VTK error code " << writer->GetErrorCode() << ").";
    }
    else
    {
      vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
      // VTK 9 writes legacy version 5.1 by default, which VTK 8 based readers (earlier MITK
      // Diffusion releases, other tools) cannot parse; 4.2 is read by both.
      writer->SetFileVersion(vtkDataWriter::VTK_LEGACY_READER_VERSION_4_2);
      writer->SetInputData(fibPoly);
      writer->SetFileName(filename.c_str());
      if (us::any_cast<bool>(options["Save as binary file"]))
      {
        MITK_INFO << "Writing fiber bundle as vtk binary file";
        writer->SetFileTypeToBinary();
      }
      else
      {
        MITK_INFO << "Writing fiber bundle as vtk ascii file";
        writer->SetFileTypeToASCII();
      }
      if (writer->Write() != 1 || writer->GetErrorCode() != 0)
        mitkThrow() << "Could not write fiber bundle to " << filename << " (VTK error code " << writer->GetErrorCode() << ").";
    }

    MITK_INFO << "VTK Fiber bundle written to " << filename;
  }
  catch(...)
  {
    throw;
  }
}
