/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkOrientationDistributionFunction_h
#define itkOrientationDistributionFunction_h

#include "itkFixedArray.h"
#include "itkMatrix.h"
#include <map>
#include <mutex>
#include "itkDiffusionTensor3D.h"

#include "vtkPolyData.h"
#include <vtkSmartPointer.h>
#include "mitkPointShell.h"

// generate by n-fold subdivisions of an icosahedron
// 1 - 12
// 2 - 42
// 3 - 92
// 4 - 162
// 5 - 252
// 6 - 362
// 7 - 492
// 8 - 642
// 9 - 812
// 10 - 1002
#define ODF_SAMPLING_SIZE 252

namespace itk
{

/** \class OrientationDistributionFunction
 * \brief Represents an ODF
 *
 * Reference: David S. Tuch, Q-ball imaging,
 * Magnetic Resonance in Medicine Volume 52 Issue 6, Pages 1358 - 1372
 *
 * \author Klaus Fritzsche, MBI
 *
 */

template < typename TComponent, unsigned int NOdfDirections >
class OrientationDistributionFunction: public
        FixedArray<TComponent, NOdfDirections>
{
public:

  /** Standard class typedefs. */
  typedef OrientationDistributionFunction  Self;
  typedef FixedArray<TComponent,NOdfDirections> Superclass;

  /** Dimension of the vector space. */
  itkStaticConstMacro(InternalDimension, unsigned int, NOdfDirections);

  /** Convenience typedefs. */
  typedef FixedArray<TComponent,
    itkGetStaticConstMacro(InternalDimension)> BaseArray;

  /**  Define the component type. */
  typedef TComponent ComponentType;
  typedef typename Superclass::ValueType ValueType;
  typedef typename NumericTraits<ValueType>::RealType AccumulateValueType;
  typedef typename NumericTraits<ValueType>::RealType RealValueType;


  typedef vnl_matrix_fixed<double, 3, NOdfDirections> DirectionsType;

  /** Default constructor has nothing to do. */
  OrientationDistributionFunction() {this->Fill(0); m_EigenAnalysisCalculated = false; }

  OrientationDistributionFunction (const ComponentType& r) { this->Fill(r); m_EigenAnalysisCalculated = false; }

  typedef ComponentType ComponentArrayType[ itkGetStaticConstMacro(InternalDimension) ];

  /** Pass-through constructor for the Array base class. The eigen analysis of InitFromEllipsoid
   * is copied too, so a copy reports the same principal direction as the original. */
  OrientationDistributionFunction(const Self& r): BaseArray(r) { CopyEigenAnalysis(r); }
  OrientationDistributionFunction(const ComponentArrayType r): BaseArray(r) { m_EigenAnalysisCalculated = false; }

  /** Pass-through assignment operator for the Array base class. */
  Self& operator= (const Self& r);
  Self& operator= (const ComponentType& r);
  Self& operator= (const ComponentArrayType r);

  /** Aritmetic operations between pixels. Return a new OrientationDistributionFunction. */
  Self operator+(const Self &vec) const;
  Self operator-(const Self &vec) const;
  const Self & operator+=(const Self &vec);
  const Self & operator-=(const Self &vec);

  /** Arithmetic operations between ODFs and scalars */
  Self operator*(const RealValueType & scalar ) const;
  Self operator/(const RealValueType & scalar ) const;
  const Self & operator*=(const RealValueType & scalar );
  const Self & operator/=(const RealValueType & scalar );

  /** The N sampling directions, as the columns of a 3 x N matrix (mitk::pointshell). */
  static const DirectionsType &GetDirections()
  {
    return mitk::pointshell::Directions<NOdfDirections>();
  }

  /** Return the number of components. */
  static unsigned int GetNumberOfComponents()
  {
    return itkGetStaticConstMacro(InternalDimension);
  }

  /** Return the value for the Nth component. */
  ComponentType GetNthComponent(int c) const
  { return this->operator[](c); }

  /** Set the Nth component to v. */
  void SetNthComponent(int c, const ComponentType& v)
  {  this->operator[](c) = v; }

  void InitFromTensor(itk::DiffusionTensor3D<TComponent> tensor);

  /** Evaluate diffusion tensor as ellipsoid. */
  void InitFromEllipsoid(itk::DiffusionTensor3D<TComponent> tensor);

  void Normalize();

  Self MinMaxNormalize() const;

  Self MaxNormalize() const;


  int GetPrincipalDiffusionDirectionIndex() const;

  vnl_vector_fixed<double,3> GetPrincipalDiffusionDirection() const;


  TComponent GetGeneralizedFractionalAnisotropy() const;

  TComponent GetGeneralizedGFA(int k, int p) const;

  TComponent GetNormalizedEntropy() const;

  TComponent GetNematicOrderParameter() const;

  TComponent GetStdDevByMaxValue() const;

  ComponentType GetMaxValue() const;

  ComponentType GetMeanValue() const;

  TComponent GetPrincipleCurvature(double alphaMinDegree, double alphaMaxDegree, int invert) const;

  /** The directions that share a triangle of the base mesh with direction idx, sorted. */
  static const std::vector<int> &GetNeighbors(int idx);

  /** Triangulation of the sampling directions, shared and not to be changed. */
  static vtkPolyData *GetBaseMesh();

  static vnl_vector_fixed<double,3> GetDirection(int i);

private:

  static vtkSmartPointer<vtkPolyData> BuildBaseMesh();

  /** Per pair of angles (GetPrincipleCurvature) and direction: the directions in that angular range */
  static std::map< std::pair<double, double>, std::vector< std::vector<int> > > m_AngularRangeIdxs;

  static std::mutex m_MutexAngularRange;
  typename itk::DiffusionTensor3D<TComponent>::EigenValuesArrayType   m_EigenValues;
  typename itk::DiffusionTensor3D<TComponent>::EigenVectorsMatrixType m_EigenVectors;
  bool  m_EigenAnalysisCalculated;

  void CopyEigenAnalysis(const Self& r)
  {
    m_EigenAnalysisCalculated = r.m_EigenAnalysisCalculated;
    if (m_EigenAnalysisCalculated)
    {
      m_EigenValues = r.m_EigenValues;
      m_EigenVectors = r.m_EigenVectors;
    }
  }

};

///** This extra typedef is necessary for preventing an Internal Compiler Error in
// * Microsoft Visual C++ 6.0. This typedef is not needed for any other compiler. */
//typedef std::ostream               OutputStreamType;
//typedef std::istream               InputStreamType;

//template< typename TComponent, unsigned int NOdfDirections  >
//MITKDIFFUSIONCORE_EXPORT OutputStreamType& operator<<(OutputStreamType& os,
//              const OrientationDistributionFunction<TComponent,NOdfDirections> & c);
//template< typename TComponent, unsigned int NOdfDirections  >
//MITKDIFFUSIONCORE_EXPORT InputStreamType& operator>>(InputStreamType& is,
//                    OrientationDistributionFunction<TComponent,NOdfDirections> & c);



} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOrientationDistributionFunction.txx"
#endif


#endif //itkOrientationDistributionFunction_h
