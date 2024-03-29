/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _itk_DWIVoxelFunctor_h_
#define _itk_DWIVoxelFunctor_h_

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
  itkTypeMacro(DWIVoxelFunctor, Object)
  /** Method for creation through the object factory. */
  //itkFactorylessNewMacro(Self)
  //itkCloneMacro(Self)

  //The first column of the matrix is reserved for the new calculated signal (other columns can hold e.g. the RMS-error)
  virtual void operator()(vnl_matrix<double> & /*newSignal*/,const vnl_matrix<double> & /*SignalMatrix*/, const double & /*S0*/)=0;

};

}

#endif
