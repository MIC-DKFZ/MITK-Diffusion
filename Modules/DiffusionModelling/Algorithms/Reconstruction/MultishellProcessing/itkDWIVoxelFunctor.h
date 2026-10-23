/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkDWIVoxelFunctor_h
#define itkDWIVoxelFunctor_h

#include <MitkDiffusionModellingExports.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include <itkObject.h>
#include <itkObjectFactory.h>

namespace itk
{

/**
 * \brief The DWIVoxelFunctor class
 * Abstract basisclass for voxelprocessing of Diffusion Weighted Images
 */
class MITKDIFFUSIONMODELLING_EXPORT DWIVoxelFunctor: public Object
{
protected:
  DWIVoxelFunctor(){}
  ~DWIVoxelFunctor() override{}

public:
  typedef DWIVoxelFunctor                         Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef itk::Object                                  Superclass;

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(DWIVoxelFunctor)
  /** Method for creation through the object factory. */
  //itkFactorylessNewMacro(Self)
  //itkCloneMacro(Self)

  //The first column of the matrix is reserved for the new calculated signal (other columns can hold e.g. the RMS-error)
  virtual void operator()(vnl_matrix<double> & /*newSignal*/,const vnl_matrix<double> & /*SignalMatrix*/, const double & /*S0*/)=0;

  /** b-value of the signal the functor computes for every direction */
  void setTargetBValue(const double & targetBValue){m_TargetBvalue = targetBValue;}
  /** b-value of each shell, in the order of the columns of the signal matrix */
  void setListOfBValues(const vnl_vector<double> & BValueList){m_BValueList = BValueList;}

  /** Number of model parameters fitted per direction */
  virtual unsigned int GetNumberOfUnknowns() const = 0;

  /**
   * \brief Checks the settings against the shells of the input, before any voxel is processed.
   * \throw itk::ExceptionObject without a positive target b-value, when the list of b-values does
   * not have one entry per shell, or with fewer shells than unknowns (Levenberg-Marquardt does not
   * start then, and every voxel would get the start values).
   */
  void CheckShells(unsigned int numberOfShells) const
  {
    if (!(m_TargetBvalue > 0))
      itkExceptionMacro(<< "No positive target b-value set (" << m_TargetBvalue << ").");
    if (m_BValueList.size() != numberOfShells)
      itkExceptionMacro(<< m_BValueList.size() << " b-values set for " << numberOfShells << " shells.");
    if (numberOfShells < GetNumberOfUnknowns())
      itkExceptionMacro(<< "The model has " << GetNumberOfUnknowns() << " unknowns per direction and needs as many shells; the input has "
                        << numberOfShells << ".");
  }

protected:
  double m_TargetBvalue = 0;
  vnl_vector<double> m_BValueList;
};

}

#endif
