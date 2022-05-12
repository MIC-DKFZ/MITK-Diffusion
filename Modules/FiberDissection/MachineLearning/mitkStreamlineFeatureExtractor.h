/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef StreamlineFeatureExtractor_h
#define StreamlineFeatureExtractor_h

// MITK
#include <MitkFiberTrackingExports.h>
#include <mitkPlanarEllipse.h>
#include <mitkFiberBundle.h>
#include <mitkClusteringMetric.h>

// ITK
#include <itkProcessObject.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

namespace mitk{

/**
* \brief    */

class MITKFIBERDISSECTION_EXPORT StreamlineFeatureExtractor
{
public:


  StreamlineFeatureExtractor();
  ~StreamlineFeatureExtractor();

  typedef itk::Image< float, 3 >                              FloatImageType;
  typedef itk::Image< unsigned char, 3 >                      UcharImageType;

  void Update(){
    this->GenerateData();
  }



protected:

  void GenerateData();
  std::vector< vnl_matrix<float> > ResampleFibers(FiberBundle::Pointer tractogram);



  unsigned int                                m_NumPoints;
  mitk::FiberBundle::Pointer                  m_Tractogram;
  std::vector<vnl_matrix<float> >             T_prototypes;
  std::vector<vnl_matrix<float> >             T_Tractogram;

};
}

#endif
