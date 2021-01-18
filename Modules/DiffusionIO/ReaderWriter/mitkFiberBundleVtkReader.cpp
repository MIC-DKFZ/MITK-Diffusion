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
#include <tinyxml2.h>
#include <vtkCleanPolyData.h>
#include <mitkTrackvis.h>
#include <mitkCustomMimeType.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLDataReader.h>
#include "mitkDiffusionIOMimeTypes.h"
#include <vtkUnsignedCharArray.h>


mitk::FiberBundleVtkReader::FiberBundleVtkReader()
  : mitk::AbstractFileReader( mitk::DiffusionIOMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME(), "VTK Fiber Bundle Reader" )
{
  m_ServiceReg = this->RegisterService();
}

mitk::FiberBundleVtkReader::FiberBundleVtkReader(const FiberBundleVtkReader &other)
  :mitk::AbstractFileReader(other)
{
}

mitk::FiberBundleVtkReader * mitk::FiberBundleVtkReader::Clone() const
{
  return new FiberBundleVtkReader(*this);
}


std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleVtkReader::DoRead()
{

  std::vector<itk::SmartPointer<mitk::BaseData> > result;

  const std::string& locale = "C";
  const std::string& currLocale = setlocale( LC_ALL, nullptr );
  setlocale(LC_ALL, locale.c_str());

  std::string filename = this->GetInputLocation();

  std::string ext = itksys::SystemTools::GetFilenameLastExtension(filename);
  ext = itksys::SystemTools::LowerCase(ext);

  try
  {
    MITK_INFO << "Loading tractogram (VTK format): " << itksys::SystemTools::GetFilenameName(filename);
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName( this->GetInputLocation().c_str() );

    if (reader->IsFilePolyData())
    {
      reader->Update();

      if ( reader->GetOutput() != nullptr )
      {
        vtkSmartPointer<vtkPolyData> fiberPolyData = reader->GetOutput();
        FiberBundle::Pointer fiberBundle = FiberBundle::New(fiberPolyData);

        vtkSmartPointer<vtkFloatArray> weights = vtkFloatArray::SafeDownCast(fiberPolyData->GetCellData()->GetArray("FIBER_WEIGHTS"));
        if (weights!=nullptr)
        {
          //                    float weight=0;
          for (int i=0; i<weights->GetNumberOfValues(); i++)
          {
            if (weights->GetValue(i)<0.0)
            {
              MITK_ERROR << "Fiber weight<0 detected! Setting value to 0.";
              weights->SetValue(i,0);
            }
          }
          //                        if (!mitk::Equal(weights->GetValue(i),weight,0.00001))
          //                        {
          //                            MITK_INFO << "Weight: " << weights->GetValue(i);
          //                            weight = weights->GetValue(i);
          //                        }
          fiberBundle->SetFiberWeights(weights);
        }

        vtkSmartPointer<vtkUnsignedCharArray> fiberColors = vtkUnsignedCharArray::SafeDownCast(fiberPolyData->GetPointData()->GetArray("FIBER_COLORS"));
        if (fiberColors!=nullptr)
          fiberBundle->SetFiberColors(fiberColors);

        result.push_back(fiberBundle.GetPointer());
        return result;
      }
    }
    else
      MITK_INFO << "File is not VTK format.";
  }
  catch(...)
  {
    throw;
  }

  try
  {
    MITK_INFO << "Trying to load fiber file as VTP format.";
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName( this->GetInputLocation().c_str() );

    if ( reader->CanReadFile(this->GetInputLocation().c_str()) )
    {
      reader->Update();

      if ( reader->GetOutput() != nullptr )
      {
        vtkSmartPointer<vtkPolyData> fiberPolyData = reader->GetOutput();
        FiberBundle::Pointer fiberBundle = FiberBundle::New(fiberPolyData);

        vtkSmartPointer<vtkFloatArray> weights = vtkFloatArray::SafeDownCast(fiberPolyData->GetCellData()->GetArray("FIBER_WEIGHTS"));

        if (weights!=nullptr)
        {
          //                                float weight=0;
          //                                for (int i=0; i<weights->GetSize(); i++)
          //                                    if (!mitk::Equal(weights->GetValue(i),weight,0.00001))
          //                                    {
          //                                        MITK_INFO << "Weight: " << weights->GetValue(i);
          //                                        weight = weights->GetValue(i);
          //                                    }
          fiberBundle->SetFiberWeights(weights);
        }

        vtkSmartPointer<vtkUnsignedCharArray> fiberColors = vtkUnsignedCharArray::SafeDownCast(fiberPolyData->GetPointData()->GetArray("FIBER_COLORS"));
        if (fiberColors!=nullptr)
          fiberBundle->SetFiberColors(fiberColors);

        result.push_back(fiberBundle.GetPointer());
        return result;
      }
    }
    else
      MITK_INFO << "File is not VTP format.";
  }
  catch(...)
  {
    throw;
  }

  throw "Selected file is no vtk readable fiber format (binary or ascii vtk or vtp file).";

  return result;
}
