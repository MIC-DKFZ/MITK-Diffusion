/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkOrientationDistributionFunction_h
#define __itkOrientationDistributionFunction_h

#include "itkIndent.h"
#include "itkFixedArray.h"
#include "itkMatrix.h"
#include "itkSymmetricEigenAnalysis.h"
#include <mutex>
#include "itkDiffusionTensor3D.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDelaunay2D.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkPlane.h"

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

  enum InterpolationMethods {
    ODF_NEAREST_NEIGHBOR_INTERP,
    ODF_TRILINEAR_BARYCENTRIC_INTERP,
    ODF_SPHERICAL_GAUSSIAN_BASIS_FUNCTIONS
  };

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

  typedef Matrix<TComponent, NOdfDirections, NOdfDirections> MatrixType;

  typedef vnl_matrix_fixed<double, 3, NOdfDirections> DirectionsType;

  /** Default constructor has nothing to do. */
  OrientationDistributionFunction() {this->Fill(0); m_EigenAnalysisCalculated = false; }

  OrientationDistributionFunction (const ComponentType& r) { this->Fill(r); m_EigenAnalysisCalculated = false; }

  typedef ComponentType ComponentArrayType[ itkGetStaticConstMacro(InternalDimension) ];

  /** Pass-through constructor for the Array base class. */
  OrientationDistributionFunction(const Self& r): BaseArray(r) { m_EigenAnalysisCalculated = false; }
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

  /** Return the number of components. */
  static DirectionsType* GetDirections()
  {
    return itkGetStaticConstMacro(m_Directions);
  }

  /** Return the number of components. */
  static unsigned int GetNumberOfComponents()
  {
    return itkGetStaticConstMacro(InternalDimension);
  }

  /** Return the value for the Nth component. */
  ComponentType GetNthComponent(int c) const
  { return this->operator[](c); }

  /** Return the value for the Nth component. */
  ComponentType GetInterpolatedComponent( vnl_vector_fixed<double,3> dir, InterpolationMethods method ) const;

  /** Set the Nth component to v. */
  void SetNthComponent(int c, const ComponentType& v)
  {  this->operator[](c) = v; }

  /** Matrix notation, in const and non-const forms. */
  ValueType & operator()( unsigned int row, unsigned int col );
  const ValueType & operator()( unsigned int row, unsigned int col ) const;

  /** Set the distribution to isotropic.*/
  void SetIsotropic();

  void InitFromTensor(itk::DiffusionTensor3D<TComponent> tensor);

  /** Evaluate diffusion tensor as ellipsoid. */
  void InitFromEllipsoid(itk::DiffusionTensor3D<TComponent> tensor);

  /** Pre-Multiply by a Matrix as ResultingTensor = Matrix * ThisTensor. */
  Self PreMultiply( const MatrixType & m ) const;

  /** Post-Multiply by a Matrix as ResultingTensor = ThisTensor * Matrix. */
  Self PostMultiply( const MatrixType & m ) const;

  void Normalize();

  Self MinMaxNormalize() const;

  Self MaxNormalize() const;

  void L2Normalize();

  int GetPrincipalDiffusionDirectionIndex() const;

  vnl_vector_fixed<double,3> GetPrincipalDiffusionDirection() const;

  int GetNthDiffusionDirection(int n, vnl_vector_fixed<double,3> rndVec) const;

  TComponent GetGeneralizedFractionalAnisotropy() const;

  TComponent GetGeneralizedGFA(int k, int p) const;

  TComponent GetNormalizedEntropy() const;

  TComponent GetNematicOrderParameter() const;

  TComponent GetStdDevByMaxValue() const;

  ComponentType GetMaxValue() const;

  ComponentType GetMinValue() const;

  ComponentType GetMeanValue() const;

  TComponent GetPrincipleCurvature(double alphaMinDegree, double alphaMaxDegree, int invert) const;

  static std::vector<int> GetNeighbors(int idx);

  static vtkPolyData* GetBaseMesh(){ComputeBaseMesh(); return m_BaseMesh;}

  static void ComputeBaseMesh();

  static double GetMaxChordLength();

  static vnl_vector_fixed<double,3> GetDirection(int i);

private:

  static vtkPolyData* m_BaseMesh;

  static double m_MaxChordLength;

  static DirectionsType* m_Directions;

  static std::vector< std::vector<int>* >* m_NeighborIdxs;

  static std::vector< std::vector<int>* >* m_AngularRangeIdxs;

  static std::vector<int>* m_HalfSphereIdxs;

  static std::mutex m_MutexBaseMesh;
  static std::mutex m_MutexHalfSphereIdxs;
  static std::mutex m_MutexNeighbors;
  static std::mutex m_MutexAngularRange;
  typename itk::DiffusionTensor3D<TComponent>::EigenValuesArrayType   m_EigenValues;
  typename itk::DiffusionTensor3D<TComponent>::EigenVectorsMatrixType m_EigenVectors;
  bool  m_EigenAnalysisCalculated;

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


#endif //__itkOrientationDistributionFunction_h
