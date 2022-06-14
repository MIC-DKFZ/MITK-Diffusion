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


#include "MitkFiberDissectionExports.h"
// MITK
#include <mitkPlanarEllipse.h>
#include <mitkFiberBundle.h>

// ITK
#include <itkProcessObject.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

// OpenCV
#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

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

  void SetTractogramPlus(const mitk::FiberBundle::Pointer &Tractogram);
  void SetTractogramMinus(const mitk::FiberBundle::Pointer &Tractogram);

  void SetTractogramTest(const mitk::FiberBundle::Pointer &Tractogram, std::string TractogramTestName);
  void CreateClassifier();
  std::vector<std::vector<unsigned int>>  GetData();

//  void CreatePrediction(std::vector<unsigned int> &index);
  mitk::FiberBundle::Pointer CreatePrediction(std::vector<unsigned int> &index);

  mitk::FiberBundle::Pointer                  m_Prediction;
  mitk::FiberBundle::Pointer                  m_ToLabel;

  std::vector<std::vector<unsigned int>>      m_index;


protected:

  void GenerateData();
  std::vector< vnl_matrix<float> > ResampleFibers(FiberBundle::Pointer tractogram);
  std::vector<vnl_matrix<float> > CalculateDmdf(std::vector<vnl_matrix<float> > tractogram, std::vector<vnl_matrix<float> > prototypes);



  unsigned int                                m_NumPoints;
  mitk::FiberBundle::Pointer                  m_TractogramPlus;
  mitk::FiberBundle::Pointer                  m_TractogramMinus;
  mitk::FiberBundle::Pointer                  m_TractogramTest;
  std::string                                 m_DistancesTestName;
  std::vector<vnl_matrix<float> >             T_Prototypes;
  std::vector<vnl_matrix<float> >             T_TractogramPlus;
  std::vector<vnl_matrix<float> >             T_TractogramMinus;
  std::vector<vnl_matrix<float> >             T_TractogramTest;
  std::vector<vnl_matrix<float> >             m_DistancesPlus;
  std::vector<vnl_matrix<float> >             m_DistancesMinus;
  std::vector<vnl_matrix<float> >             m_DistancesTest;
  cv::Ptr<cv::ml::TrainData>                  m_traindata;
};
}

#endif
