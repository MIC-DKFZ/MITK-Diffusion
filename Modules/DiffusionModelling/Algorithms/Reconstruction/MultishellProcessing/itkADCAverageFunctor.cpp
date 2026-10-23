/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "itkADCAverageFunctor.h"
#include <cmath>
#include <iostream>
#include <iomanip>

void itk::ADCAverageFunctor::operator()(vnl_matrix<double> & newSignal, const vnl_matrix<double> & SignalMatrix, const double & S0)
{

  vnl_matrix<double> ADCMatrix(SignalMatrix.rows(),SignalMatrix.cols());

  // Calculate ADC for each measurement
  for(unsigned int i = 0 ; i < SignalMatrix.rows(); ++i)
    for(unsigned int j = 0; j < SignalMatrix.cols(); ++j)
      ADCMatrix(i,j) = std::log(SignalMatrix(i,j) / S0) / (-m_BValueList[j]);// D = ln(S/S0)/-b


  // Calculate new Signal using Average ADC
  for(unsigned int i = 0 ; i < SignalMatrix.rows(); ++i){
    double averageADC = ADCMatrix.get_row(i).mean();
    newSignal.put(i,0, S0 * std::exp(-m_TargetBvalue * averageADC) ); // S = S0*exp(-b*D)

    // Root Mean Squares Error
    double error = 0;
    for(unsigned int j = 0 ; j < SignalMatrix.cols(); ++j)
      error += std::pow(SignalMatrix(i,j) - S0 * std::exp(-m_BValueList[j] * averageADC),2); // sum of squres
    error /= (double)SignalMatrix.cols(); // mean
    error = std::sqrt(error);

    newSignal.put(i, 1, error ); // RMS Error
  }
}
