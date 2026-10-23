/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleVtkReader.h"
#include "mitkFiberBundleMimeTypes.h"

#include <itksys/SystemTools.hxx>
#include <mitkFiberBundle.h>
#include <mitkGeometry3D.h>
#include <mitkLocaleSwitch.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkLogger.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkTransformFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLPolyDataReader.h>

#include <cstdint>
#include <string>
#include <thread>

namespace
{
  /** Collects the warnings and errors that VTK logs in the constructing thread while the object
   *  exists. VTK's legacy reader reports a file that ends inside its data only this way
   *  ("Error reading binary data!") and returns the points and arrays read so far, or allocated
   *  but partly filled. */
  class VtkProblemLog
  {
  public:
    VtkProblemLog()
      : m_Id("mitk::FiberBundleVtkReader " + std::to_string(reinterpret_cast<std::uintptr_t>(this)))
      , m_Thread(std::this_thread::get_id())
    {
      vtkLogger::AddCallback(m_Id.c_str(), &VtkProblemLog::OnMessage, this, vtkLogger::VERBOSITY_WARNING);
    }

    ~VtkProblemLog()
    {
      vtkLogger::RemoveCallback(m_Id.c_str());
    }

    VtkProblemLog(const VtkProblemLog&) = delete;
    VtkProblemLog& operator=(const VtkProblemLog&) = delete;

    const std::string& GetMessages() const { return m_Messages; }

  private:
    static void OnMessage(void* userData, const vtkLogger::Message& message)
    {
      auto self = static_cast<VtkProblemLog*>(userData);
      if (std::this_thread::get_id() == self->m_Thread && message.verbosity <= vtkLogger::VERBOSITY_WARNING)
        self->m_Messages += std::string(message.message != nullptr ? message.message : "") + " ";
    }

    std::string m_Id;
    std::thread::id m_Thread;
    std::string m_Messages;
  };
}


mitk::FiberBundleVtkReader::FiberBundleVtkReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME(), "VTK Fiber Bundle Reader" )
{
  Options defaultOptions;
  defaultOptions["Use RAS space"] = false;
  this->SetDefaultOptions(defaultOptions);

  this->RegisterService();
}

mitk::FiberBundleVtkReader * mitk::FiberBundleVtkReader::Clone() const
{
  return new FiberBundleVtkReader(*this);
}


std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleVtkReader::DoRead()
{
  std::vector<itk::SmartPointer<mitk::BaseData> > result;

  // ASCII files are parsed with the decimal point of the current locale
  mitk::LocaleSwitch localeSwitch("C");

  const std::string filename = this->GetInputLocation();
  Options options = this->GetOptions();
  const bool ras = us::any_cast<bool>(options["Use RAS space"]);

  MITK_INFO << "Loading tractogram (VTK format): " << itksys::SystemTools::GetFilenameName(filename);
  vtkSmartPointer<vtkPolyData> fiberPolyData;
  std::string problems;
  vtkSmartPointer<vtkPolyDataReader> legacyReader = vtkSmartPointer<vtkPolyDataReader>::New();
  legacyReader->SetFileName(filename.c_str());
  if (legacyReader->IsFilePolyData())
  {
    {
      VtkProblemLog log;
      legacyReader->Update();
      problems = log.GetMessages();
    }
    if (legacyReader->GetErrorCode() != 0)
      mitkThrow() << "Could not read " << filename << " (VTK error code " << legacyReader->GetErrorCode() << ").";
    fiberPolyData = legacyReader->GetOutput();
  }
  else
  {
    vtkSmartPointer<vtkXMLPolyDataReader> xmlReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    if (!xmlReader->CanReadFile(filename.c_str()))
      mitkThrow() << filename << " is no VTK polydata file (legacy .fib/.vtk or XML .vtp).";
    xmlReader->SetFileName(filename.c_str());
    {
      VtkProblemLog log;
      xmlReader->Update();
      problems = log.GetMessages();
    }
    if (xmlReader->GetErrorCode() != 0)
      mitkThrow() << "Could not read " << filename << " (VTK error code " << xmlReader->GetErrorCode() << ").";
    fiberPolyData = xmlReader->GetOutput();
  }
  if (fiberPolyData == nullptr)
    mitkThrow() << "Could not read " << filename;
  if (!problems.empty())
    mitkThrow() << "Could not read " << filename << " (truncated or corrupt): " << problems;
  // A legacy ASCII file that ends inside its points gives them without any fiber and no message.
  if (fiberPolyData->GetNumberOfPoints() > 0 && fiberPolyData->GetNumberOfCells() == 0)
    mitkThrow() << "Could not read " << filename << " (truncated or corrupt): points without fibers.";
  if (fiberPolyData->GetNumberOfPolys() > 0 || fiberPolyData->GetNumberOfStrips() > 0)
    mitkThrow() << filename << " contains polygons or triangle strips; a fiber bundle consists of lines.";

  if (ras)
  {
    vtkSmartPointer< vtkMatrix4x4 > matrix = vtkSmartPointer< vtkMatrix4x4 >::New();
    matrix->Identity();
    matrix->SetElement(0,0,-1);
    matrix->SetElement(1,1,-1);
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(fiberPolyData);
    transformFilter->SetTransform(geometry->GetVtkTransform());
    transformFilter->Update();
    fiberPolyData = transformFilter->GetPolyDataOutput();
  }

  FiberBundle::Pointer fiberBundle = FiberBundle::New(fiberPolyData);
  fiberBundle->setIsRAS(ras);

  // The arrays are taken from the bundle's own polydata: vtkCleanPolyData in the constructor drops
  // unused points and moves single-point fibers first, and carries the arrays along.
  vtkPolyData* cleaned = fiberBundle->GetFiberPolyData();

  vtkFloatArray* fileWeights = vtkFloatArray::SafeDownCast(cleaned->GetCellData()->GetArray("FIBER_WEIGHTS"));
  if (fileWeights != nullptr && fileWeights->GetNumberOfValues() != static_cast<vtkIdType>(fiberBundle->GetNumFibers()))
  {
    MITK_WARN << "Ignoring FIBER_WEIGHTS of " << filename << ": " << fileWeights->GetNumberOfValues() << " values for "
              << fiberBundle->GetNumFibers() << " fibers.";
  }
  else if (fileWeights != nullptr)
  {
    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
    weights->DeepCopy(fileWeights);
    for (vtkIdType i=0; i<weights->GetNumberOfValues(); i++)
    {
      if (!(weights->GetValue(i) >= 0.0f)) // negative or NaN
      {
        MITK_ERROR << "Fiber weight<0 or NaN detected! Setting value to 0.";
        weights->SetValue(i, 0);
      }
    }
    fiberBundle->SetFiberWeights(weights);
  }

  vtkUnsignedCharArray* fileColors = vtkUnsignedCharArray::SafeDownCast(cleaned->GetPointData()->GetArray("FIBER_COLORS"));
  if (fileColors != nullptr)
  {
    const int components = fileColors->GetNumberOfComponents();
    if (fileColors->GetNumberOfTuples() != cleaned->GetNumberOfPoints() || (components != 3 && components != 4))
    {
      MITK_WARN << "Ignoring FIBER_COLORS of " << filename << ": " << fileColors->GetNumberOfTuples() << " tuples of "
                << components << " components for " << cleaned->GetNumberOfPoints() << " points.";
    }
    else
    {
      // RGB colours get an opaque alpha
      vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
      colors->SetNumberOfComponents(4);
      colors->SetNumberOfTuples(fileColors->GetNumberOfTuples());
      for (vtkIdType i=0; i<fileColors->GetNumberOfTuples(); ++i)
        for (int c=0; c<4; ++c)
          colors->SetValue(4*i+c, c<components ? fileColors->GetValue(components*i+c) : 255);
      fiberBundle->SetFiberColors(colors);
    }
  }

  result.push_back(fiberBundle.GetPointer());
  return result;
}
