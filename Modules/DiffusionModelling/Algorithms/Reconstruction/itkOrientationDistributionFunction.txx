/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _itkOrientationDistributionFunction_txx
#define _itkOrientationDistributionFunction_txx

#include <vnl/algo/vnl_matrix_inverse.h>
#include <itkPointShell.h>
#include <itkMatrix.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_inverse.h>
#include <ciso646>
#include <vnl/vnl_vector_fixed.hxx>
#include <vnl/vnl_cross.h>
#include <mitkLogMacros.h>
#include <vtkDelaunay2D.h>
#include <vtkPlane.h>

namespace itk
{

template<class T, unsigned int N>
vtkPolyData* itk::OrientationDistributionFunction<T,N>::m_BaseMesh = nullptr;

template<class T, unsigned int N>
double itk::OrientationDistributionFunction<T,N>::m_MaxChordLength = -1.0;

template<class T, unsigned int N>
vnl_matrix_fixed<double, 3, N>* itk::OrientationDistributionFunction<T,N>::m_Directions
= itk::PointShell<N, vnl_matrix_fixed<double, 3, N> >::DistributePointShell();

template<class T, unsigned int N>
std::vector< std::vector<int>* >* itk::OrientationDistributionFunction<T,N>::m_NeighborIdxs = nullptr;

template<class T, unsigned int N>
std::vector< std::vector<int>* >* itk::OrientationDistributionFunction<T,N>::m_AngularRangeIdxs = nullptr;

template<class T, unsigned int N>
std::vector<int>* itk::OrientationDistributionFunction<T,N>::m_HalfSphereIdxs = nullptr;

template<class T, unsigned int N>
std::mutex itk::OrientationDistributionFunction<T,N>::m_MutexBaseMesh;
template<class T, unsigned int N>
std::mutex itk::OrientationDistributionFunction<T,N>::m_MutexHalfSphereIdxs;
template<class T, unsigned int N>
std::mutex itk::OrientationDistributionFunction<T,N>::m_MutexNeighbors;
template<class T, unsigned int N>
std::mutex itk::OrientationDistributionFunction<T,N>::m_MutexAngularRange;

/**
  * Assignment Operator
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>&
OrientationDistributionFunction<T, NOdfDirections>
::operator= (const Self& r)
{
  BaseArray::operator=(r);
  return *this;
}

/**
  * Assignment Operator from a scalar constant
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>&
OrientationDistributionFunction<T, NOdfDirections>
::operator= (const ComponentType & r)
{
  BaseArray::operator=(&r);
  return *this;
}

/**
  * Assigment from a plain array
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>&
OrientationDistributionFunction<T, NOdfDirections>
::operator= (const ComponentArrayType r )
{
  BaseArray::operator=(r);
  return *this;
}

/**
  * Returns a temporary copy of a vector
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::operator+(const Self & r) const
{
  Self result;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    result[i] = (*this)[i] + r[i];
  }
  return result;
}

/**
  * Returns a temporary copy of a vector
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::operator-(const Self & r) const
{
  Self result;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    result[i] = (*this)[i] - r[i];
  }
  return result;
}

/**
  * Performs addition in place
  */
template<class T, unsigned int NOdfDirections>
const OrientationDistributionFunction<T, NOdfDirections> &
OrientationDistributionFunction<T, NOdfDirections>
::operator+=(const Self & r)
{
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    (*this)[i] += r[i];
  }
  return *this;
}

/**
  * Performs subtraction in place
  */
template<class T, unsigned int NOdfDirections>
const OrientationDistributionFunction<T, NOdfDirections> &
OrientationDistributionFunction<T, NOdfDirections>
::operator-=(const Self & r)
{
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    (*this)[i] -= r[i];
  }
  return *this;
}

/**
  * Performs multiplication by a scalar, in place
  */
template<class T, unsigned int NOdfDirections>
const OrientationDistributionFunction<T, NOdfDirections> &
OrientationDistributionFunction<T, NOdfDirections>
::operator*=(const RealValueType & r)
{
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    (*this)[i] *= r;
  }
  return *this;
}

/**
  * Performs division by a scalar, in place
  */
template<class T, unsigned int NOdfDirections>
const OrientationDistributionFunction<T, NOdfDirections> &
OrientationDistributionFunction<T, NOdfDirections>
::operator/=(const RealValueType & r)
{
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    (*this)[i] /= r;
  }
  return *this;
}

/**
  * Performs multiplication with a scalar
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::operator*(const RealValueType & r) const
{
  Self result;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    result[i] = (*this)[i] * r;
  }
  return result;
}

/**
  * Performs division by a scalar
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::operator/(const RealValueType & r) const
{
  Self result;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    result[i] = (*this)[i] / r;
  }
  return result;
}

/**
  * Matrix notation access to elements
  */
template<class T, unsigned int NOdfDirections>
const typename OrientationDistributionFunction<T, NOdfDirections>::ValueType &
OrientationDistributionFunction<T, NOdfDirections>
::operator()(unsigned int row, unsigned int col) const
{
  unsigned int k;

  if( row < col )
  {
    k = row * InternalDimension + col - row * ( row + 1 ) / 2;
  }
  else
  {
    k = col * InternalDimension + row - col * ( col + 1 ) / 2;
  }


  if( k >= InternalDimension )
  {
    k = 0;
  }

  return (*this)[k];
}

/**
  * Matrix notation access to elements
  */
template<class T, unsigned int NOdfDirections>
typename OrientationDistributionFunction<T, NOdfDirections>::ValueType &
OrientationDistributionFunction<T, NOdfDirections>
::operator()(unsigned int row, unsigned int col)
{
  unsigned int k;

  if( row < col )
  {
    k = row * InternalDimension + col - row * ( row + 1 ) / 2;
  }
  else
  {
    k = col * InternalDimension + row - col * ( col + 1 ) / 2;
  }


  if( k >= InternalDimension )
  {
    k = 0;
  }

  return (*this)[k];
}

/**
  * Set the Tensor to an Identity.
  * Set ones in the diagonal and zeroes every where else.
  */
template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::SetIsotropic()
{
  this->Fill(NumericTraits< T >::One / NOdfDirections);
}

/**
  * InitFromTensor()
  */
template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::InitFromTensor(itk::DiffusionTensor3D<T> tensor)
{
  m_EigenAnalysisCalculated = false;
  for(unsigned int i=0; i<NOdfDirections; i++)
  {
    /*
      *               | t0  t1  t2 |    g0
      *  g0 g1 g2  *  | t1  t3  t4 | *  g1
      *               | t2  t4  t5 |    g2
      *
      * =   g0 * (t0g0*t1g1*t2g2)
      *   + g1 * (t1g0+t3g1+t4g2)
      *   + g2 * (t2g0+t4g1+t5g2)
      */
    T g0 = (*m_Directions)(0,i);
    T g1 = (*m_Directions)(1,i);
    T g2 = (*m_Directions)(2,i);
    T t0 = tensor[0];
    T t1 = tensor[1];
    T t2 = tensor[2];
    T t3 = tensor[3];
    T t4 = tensor[4];
    T t5 = tensor[5];
    (*this)[i] = g0 * (t0*g0+t1*g1+t2*g2)
                 + g1 * (t1*g0+t3*g1+t4*g2)
                 + g2 * (t2*g0+t4*g1+t5*g2);

    if ((*this)[i]<0 || (*this)[i]!=(*this)[i])
      (*this)[i] = 0;
  }
}

/**
  * InitFromEllipsoid()
  */
template<class T, unsigned int NOdfDirections>
void OrientationDistributionFunction<T, NOdfDirections>::
InitFromEllipsoid( itk::DiffusionTensor3D<T> tensor )
{
  m_EigenAnalysisCalculated = false;
  FixedArray<T, 6> nulltensor;
  nulltensor.Fill(0.0);
  if( tensor == nulltensor )
  {
    for ( unsigned int it=0; it < NOdfDirections; ++it ){ (*this)[it] = (T)0; }
    MITK_DEBUG << "OrientationDistributionFunction<" << typeid(T).name() << ", " << NOdfDirections
               << ">::InitFromEllipsoid(" << typeid(tensor).name()
               << ") encountered a nulltensor as dti input point and ignorend it.";
    return;
  }

  tensor.ComputeEigenAnalysis( m_EigenValues, this->m_EigenVectors ); // gives normalized eigenvectors as lines i.e. rows.
  m_EigenAnalysisCalculated = true;

  double a = m_EigenValues[0]; // those eigenvalues are the 3 |axes of the ellipsoid|,
  double b = m_EigenValues[1]; // ComputeEigenAnalysis gives eigenValues in ascending < order per default,
  double c = m_EigenValues[2]; // therefor the third eigenVector is the main direction of diffusion.

  if( a <= 0.0 || b <= 0.0 || c <= 0.0 )
  {
    for ( unsigned int it=0; it < NOdfDirections; ++it ){ (*this)[it] = (T)0; }
    MITK_DEBUG << "OrientationDistributionFunction<" << typeid(T).name() << ", " << NOdfDirections
               << ">::InitFromEllipsoid(" << typeid(tensor).name()
               << ") encountered an eigenvalue <= 0 and ignored this input point.";
    return;
  }

  // check magnitude and scale towards 1 to minimize numerical condition kappa:
#ifdef _MSC_VER
#if _MSC_VER <= 1700
  int exponent_a = floor(std::log(a)/std::log(2));
  int exponent_b = floor(std::log(b)/std::log(2));
  int exponent_c = floor(std::log(c)/std::log(2));
#else
  int exponent_a = std::ilogb(a);
  int exponent_b = std::ilogb(b);
  int exponent_c = std::ilogb(c);
#endif
#else
  int exponent_a = std::ilogb(a);
  int exponent_b = std::ilogb(b);
  int exponent_c = std::ilogb(c);
#endif

  T min_exponent= fmin(exponent_a, fmin(exponent_b, exponent_c) );
  T max_exponent= fmax(exponent_c, fmax(exponent_b, exponent_a) );
  int scale_exponent = floor(0.5 * (min_exponent + max_exponent));
  double scaling = pow(2, scale_exponent);
  a= a/scaling;
  b= b/scaling;
  c= c/scaling;

  vnl_matrix_fixed<double, 3, 3> eigenBase; // for change of base system.
  for (int row = 0 ; row < 3; ++row) // Transposing because ComputeEigenAnalysis(,) gave us _row_ vectors.
  {
    for (int col = 0; col < 3; ++col)
    {
      eigenBase(row, col) = this->m_EigenVectors(col, row);
    }
  }
  eigenBase= vnl_inverse(eigenBase); // Assuming canonical orthonormal system x=[1;0;0];y=[0;1;0];z=[0;0;1] for original DT.
  eigenBase.assert_finite();

#ifndef NDEBUG
  double kappa=1.0;
  { // calculate numerical condition kappa= ||f(x)-f(x~)|| approximately:
    double gxaa = pow( a, -2.0);  double gybb = pow( b, -2.0);  double gzcc = pow( c, -2.0);
    kappa = sqrt( pow( a, 2.0)+ pow( b, 2.0)+ pow( c, 2.0) + 1 ) / (gxaa + gybb + gzcc)
            * sqrt( pow( a, -6.0)+ pow( b, -6.0)+ pow( c, -6.0) + pow( a, -4.0)+ pow( b, -4.0)+ pow( c, -4.0) );
    MITK_DEBUG <<"kappa= "<< kappa << ", eigenvalues= [" << a <<", "<< b <<", "<< c <<"], eigenbase= ["
              <<eigenBase(0,0)<<", "<<eigenBase(1,0)<<", "<<eigenBase(2,0)<<"; "
             <<eigenBase(0,1)<<", "<<eigenBase(1,1)<<", "<<eigenBase(2,1)<<"; "
            <<eigenBase(0,2)<<", "<<eigenBase(1,2)<<", "<<eigenBase(2,2)<<"] ";
    if( std::isnan(kappa) )
    {
      MITK_DEBUG << "oh noes: kappa was NaN, setting kappa to 1e5"; // typical value kappa=1e5 for 1e-6<=a,b,c<=1e-4.
      kappa=1e5;
    };
  }
#endif

  for( unsigned int i=0; i < NOdfDirections; ++i )
  {
    /// calculate probability p(g)=r as ellipsoids magnitude in direction of g' = B^-1 * g:
    /// (g0*r/evx)² + (g1*r/evy)² + (g2*r/evz)² = 1, |g'|=1.

    vnl_vector_fixed<double, 3> g( (*m_Directions)(0,i), (*m_Directions)(1,i), (*m_Directions)(2,i) );
    g = eigenBase*g;   // passive change of base system of g.
    g = g.normalize(); // unit vectors necessary.
    (*this)[i] = scaling / sqrt( (g[0]/a)*(g[0]/a) + (g[1]/b)*(g[1]/b) + (g[2]/c)*(g[2]/c) );

#ifndef NDEBUG
    { // boundary check for numerical stability, assuming sigma=6, ||f(x~)-f~(x~)|| <= eps*kappa*sigma.
      T min_ev= fmin(a, fmin(b, c) );   T max_ev= fmax(c, fmax(b, a) );
      double eps= std::numeric_limits<T>::epsilon();
      assert( scaling*min_ev <= ((*this)[i] + eps*kappa*6.0) ); // we should be between smallest and
      assert( (*this)[i] <= (scaling*max_ev + eps*kappa*6.0) ); // biggest eigenvalue.
    }
#endif
    if ( (*this)[i] < T(0) || (*this)[i] > T(1) || std::isnan((*this)[i]) ) // P∈[0;1] sanity check.
    { // C: NaN != NaN, C++11: isnan((*this)[i]).
      MITK_DEBUG << "OrientationDistributionFunction<" << typeid(T).name() << ", " << NOdfDirections
                 << ">::InitFromEllipsoid(" << typeid(tensor).name()
                 << ") encountered a probability value out of range [0;1] and set it to zero: (*this)["
                 << i <<"]= " << (*this)[i];
      (*this)[i] = T(0);
    }
  }
}

/**
  * L2-Normalization
  */
template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::L2Normalize()
{
  T sum = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    sum += (*this)[i]*(*this)[i];
  }
  sum = std::sqrt(sum);
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    (*this)[i] = (*this)[i] / sum;
  }
}

/**
  * Normalization to PDF
  */
template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::Normalize()
{
  T sum = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    sum += (*this)[i];
  }
  if (sum>0)
  {
    for( unsigned int i=0; i<InternalDimension; i++)
    {
      (*this)[i] = (*this)[i] / sum;
    }
  }
}

/**
  * Min/Max-Normalization
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::MinMaxNormalize() const
{
  T max = NumericTraits<T>::NonpositiveMin();
  T min = NumericTraits<T>::max();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    max = (*this)[i] > max ? (*this)[i] : max;
    min = (*this)[i] < min ? (*this)[i] : min;
  }
  Self retval;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    retval[i] = ((*this)[i] - min) / (max - min);
  }
  return retval;
}

/**
  * Max-Normalization
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::MaxNormalize() const
{
  T max = NumericTraits<T>::NonpositiveMin();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    max = (*this)[i] > max ? (*this)[i] : max;
  }
  Self retval;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    retval[i] = (*this)[i] / max;
  }
  return retval;
}

template<class T, unsigned int NOdfDirections>
T
OrientationDistributionFunction<T, NOdfDirections>
::GetMaxValue() const
{
  T max = NumericTraits<T>::NonpositiveMin();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    if((*this)[i] >= max )
    {
      max = (*this)[i];
    }
  }
  return max;
}

template<class T, unsigned int NOdfDirections>
T
OrientationDistributionFunction<T, NOdfDirections>
::GetMinValue() const
{
  T min = NumericTraits<T>::max();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    if((*this)[i] >= min )
    {
      min = (*this)[i];
    }
  }
  return min;
}

template<class T, unsigned int NOdfDirections>
T
OrientationDistributionFunction<T, NOdfDirections>
::GetMeanValue() const
{
  T sum = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
    sum += (*this)[i];
  return sum/InternalDimension;
}

template<class T, unsigned int NOdfDirections>
double
OrientationDistributionFunction<T, NOdfDirections>
::GetMaxChordLength()
{
  if(m_MaxChordLength<0.0)
  {
    ComputeBaseMesh();
    double max_dist = -1;
    vtkPoints* points = m_BaseMesh->GetPoints();
    for(int i=0; i<NOdfDirections; i++)
    {
      double p[3];
      points->GetPoint(i,p);
      std::vector<int> neighbors = GetNeighbors(i);
      for(std::size_t j=0; j<neighbors.size(); j++)
      {
        double n[3];
        points->GetPoint(neighbors[j],n);
        double d = sqrt(
                     (p[0]-n[0])*(p[0]-n[0]) +
            (p[1]-n[1])*(p[1]-n[1]) +
            (p[2]-n[2])*(p[2]-n[2]));
        max_dist = d>max_dist ? d : max_dist;
      }
    }
    m_MaxChordLength = max_dist;
  }
  return m_MaxChordLength;
}

template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::ComputeBaseMesh()
{

  m_MutexBaseMesh.lock();
  if(m_BaseMesh == nullptr)
  {

    vtkPoints* points = vtkPoints::New();
    for(unsigned int j=0; j<NOdfDirections; j++){
      double x = (*m_Directions)(0,j);
      double y = (*m_Directions)(1,j);
      double z = (*m_Directions)(2,j);
      double az = atan2(y,x);
      double elev = atan2(z,sqrt(x*x+y*y));
      double r = sqrt(x*x+y*y+z*z);
      points->InsertNextPoint(az,elev,r);
    }

    vtkPolyData* polydata = vtkPolyData::New();
    polydata->SetPoints( points );
    vtkDelaunay2D *delaunay = vtkDelaunay2D::New();
    delaunay->SetInputData( polydata );
    delaunay->Update();

    vtkCellArray* vtkpolys = delaunay->GetOutput()->GetPolys();
    vtkCellArray* vtknewpolys = vtkCellArray::New();
    vtkIdType npts; vtkIdType const *pts;
    while(vtkpolys->GetNextCell(npts,pts))
    {
      bool insert = true;
      for(int i=0; i<npts; i++)
      {
        double *tmpPoint = points->GetPoint(pts[i]);
        double az = tmpPoint[0];
        double elev = tmpPoint[1];
        if((std::abs(az)>itk::Math::pi-0.5) || (std::abs(elev)>itk::Math::pi/2-0.5))
          insert = false;
      }
      if(insert)
        vtknewpolys->InsertNextCell(npts, pts);
    }

    vtkPoints* points2 = vtkPoints::New();
    for(unsigned int j=0; j<NOdfDirections; j++){
      double x = -(*m_Directions)(0,j);
      double y = -(*m_Directions)(2,j);
      double z = -(*m_Directions)(1,j);
      double az = atan2(y,x);
      double elev = atan2(z,sqrt(x*x+y*y));
      double r = sqrt(x*x+y*y+z*z);
      points2->InsertNextPoint(az,elev,r);
    }

    vtkPolyData* polydata2 = vtkPolyData::New();
    polydata2->SetPoints( points2 );
    vtkDelaunay2D *delaunay2 = vtkDelaunay2D::New();
    delaunay2->SetInputData( polydata2 );
    delaunay2->Update();

    vtkpolys = delaunay2->GetOutput()->GetPolys();
    while(vtkpolys->GetNextCell(npts,pts))
    {
      bool insert = true;
      for(int i=0; i<npts; i++)
      {
        double *tmpPoint = points2->GetPoint(pts[i]);
        double az = tmpPoint[0];
        double elev = tmpPoint[1];
        if((std::abs(az)>itk::Math::pi-0.5) || (std::abs(elev)>itk::Math::pi/2-0.5))
          insert = false;
      }
      if(insert)
        vtknewpolys->InsertNextCell(npts, pts);
    }

    polydata->SetPolys(vtknewpolys);

    for (unsigned int p = 0; p < NOdfDirections; p++)
    {
      points->SetPoint(p,m_Directions->get_column(p).data_block());
    }
    polydata->SetPoints( points );

    m_BaseMesh = polydata;
  }
  m_MutexBaseMesh.unlock();
}

/**
  * Extract the index of the principal diffusion direction
  */
template<class T, unsigned int NOdfDirections>
int OrientationDistributionFunction<T, NOdfDirections>::GetPrincipalDiffusionDirectionIndex() const
{
  T max = NumericTraits<T>::NonpositiveMin();
  int maxidx = -1;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    if((*this)[i] >= max )
    {
      max = (*this)[i];
      maxidx = i;
    }
  }

  return maxidx;
}

/**
  * Extract the principal diffusion direction
  */
template<class T, unsigned int NOdfDirections>
vnl_vector_fixed<double,3> OrientationDistributionFunction<T, NOdfDirections>::GetPrincipalDiffusionDirection() const
{
  if (m_EigenAnalysisCalculated)
  {
    vnl_vector_fixed<double,3> vec;
    vec[0] = this->m_EigenVectors(2,0);
    vec[1] = this->m_EigenVectors(2,1);
    vec[2] = this->m_EigenVectors(2,2);
    vec.normalize();
    return vec;
  }
  else
  {
    int idx = GetPrincipalDiffusionDirectionIndex();
    if (idx>0 && idx<(int)NOdfDirections)
      return OrientationDistributionFunction<T, NOdfDirections>::GetDirection(idx);
    vnl_vector_fixed<double,3> vec; vec.fill(0);
    return vec;
  }
}

template<class T, unsigned int NOdfDirections>
std::vector<int>
OrientationDistributionFunction<T, NOdfDirections>
::GetNeighbors(int idx)
{
  ComputeBaseMesh();

  m_MutexNeighbors.lock();
  if(m_NeighborIdxs == nullptr)
  {
    m_NeighborIdxs = new std::vector< std::vector<int>* >();
    vtkCellArray* polys = m_BaseMesh->GetPolys();

    for(unsigned int i=0; i<NOdfDirections; i++)
    {
      auto      idxs = new std::vector<int>();
      polys->InitTraversal();
      vtkIdType npts; vtkIdType const *pts;
      while(polys->GetNextCell(npts,pts))
      {
        if( pts[0] == i )
        {
          idxs->push_back(pts[1]);
          idxs->push_back(pts[2]);
        }
        else if( pts[1] == i )
        {
          idxs->push_back(pts[0]);
          idxs->push_back(pts[2]);
        }
        else if( pts[2] == i )
        {
          idxs->push_back(pts[0]);
          idxs->push_back(pts[1]);
        }
      }
      std::sort(idxs->begin(), idxs->end());
      std::vector< int >::iterator endLocation;
      endLocation = std::unique( idxs->begin(), idxs->end() );
      idxs->erase(endLocation, idxs->end());
      m_NeighborIdxs->push_back(idxs);
    }
  }
  m_MutexNeighbors.unlock();

  return *m_NeighborIdxs->at(idx);
}

/**
  * Extract the n-th diffusion direction
  */
template<class T, unsigned int NOdfDirections>
int
OrientationDistributionFunction<T, NOdfDirections>
::GetNthDiffusionDirection(int n, vnl_vector_fixed<double,3> rndVec) const
{

  if( n == 0 )
    return GetPrincipalDiffusionDirectionIndex();

  m_MutexHalfSphereIdxs.lock();
  if( !m_HalfSphereIdxs )
  {
    m_HalfSphereIdxs = new std::vector<int>();
    for( unsigned int i=0; i<InternalDimension; i++)
    {
      if(dot_product(m_Directions->get_column(i),rndVec) > 0.0)
      {
        m_HalfSphereIdxs->push_back(i);
      }
    }
  }
  m_MutexHalfSphereIdxs.unlock();

  // collect indices of directions
  // that are local maxima
  std::vector<int> localMaxima;
  std::vector<int>::iterator it;
  for( it=m_HalfSphereIdxs->begin();
       it!=m_HalfSphereIdxs->end();
       it++)
  {
    std::vector<int> nbs = GetNeighbors(*it);
    std::vector<int>::iterator it2;
    bool max = true;
    for(it2 = nbs.begin();
        it2 != nbs.end();
        it2++)
    {
      if((*this)[*it2] > (*this)[*it])
      {
        max = false;
        break;
      }
    }
    if(max)
      localMaxima.push_back(*it);
  }

  // delete n highest local maxima from list
  // and return remaining highest
  int maxidx = -1;
  for( int i=0; i<=n; i++ )
  {
    maxidx = -1;
    T max = NumericTraits<T>::NonpositiveMin();
    for(it = localMaxima.begin();
        it != localMaxima.end();
        it++)
    {
      if((*this)[*it]>max)
      {
        max = (*this)[*it];
        maxidx = *it;
      }
    }
    it = find(localMaxima.begin(), localMaxima.end(), maxidx);
    if(it!=localMaxima.end())
      localMaxima.erase(it);
  }

  return maxidx;
}

template < typename TComponent, unsigned int NOdfDirections >
vnl_vector_fixed<double,3> itk::OrientationDistributionFunction<TComponent, NOdfDirections>
::GetDirection( int i )
{
  return m_Directions->get_column(i);
}

/**
  * Interpolate a position between sampled directions
  */
template<class T, unsigned int NOdfDirections>
T
OrientationDistributionFunction<T, NOdfDirections>
::GetInterpolatedComponent(vnl_vector_fixed<double,3> dir, InterpolationMethods method) const
{

  ComputeBaseMesh();
  double retval = -1.0;

  switch(method)
  {
  case ODF_NEAREST_NEIGHBOR_INTERP:
  {

    vtkPoints* points = m_BaseMesh->GetPoints();
    double current_min = NumericTraits<double>::max();
    int current_min_idx = -1;
    for(int i=0; i<NOdfDirections; i++)
    {
      vnl_vector_fixed<double,3> P(points->GetPoint(i));
      double dist = (dir-P).two_norm();
      current_min_idx = dist<current_min ? i : current_min_idx;
      current_min = dist<current_min ? dist : current_min;
    }
    retval = this->GetNthComponent(current_min_idx);
    break;

  }
  case ODF_TRILINEAR_BARYCENTRIC_INTERP:
  {

    double maxChordLength = GetMaxChordLength();
    vtkCellArray* polys = m_BaseMesh->GetPolys();
    vtkPoints* points = m_BaseMesh->GetPoints();
    vtkIdType npts; vtkIdType const *pts;
    double current_min = NumericTraits<double>::max();
    polys->InitTraversal();
    while(polys->GetNextCell(npts,pts))
    {
      vnl_vector_fixed<double,3> A(points->GetPoint(pts[0]));
      vnl_vector_fixed<double,3> B(points->GetPoint(pts[1]));
      vnl_vector_fixed<double,3> C(points->GetPoint(pts[2]));

      vnl_vector_fixed<double,3> d1;
      d1.put(0,(dir-A).two_norm());
      d1.put(1,(dir-B).two_norm());
      d1.put(2,(dir-C).two_norm());
      double maxval = d1.max_value();
      if(maxval>maxChordLength)
      {
        continue;
      }

      // Compute vectors
      vnl_vector_fixed<double,3> v0 = C - A;
      vnl_vector_fixed<double,3> v1 = B - A;

      // Project direction to plane ABC
      vnl_vector_fixed<double,3> v6 = dir;
      vnl_vector_fixed<double,3> cross = vnl_cross_3d(v0, v1);
      cross = cross.normalize();
      vtkPlane::ProjectPoint(v6.data_block(),A.data_block(),cross.data_block(),v6.data_block());
      v6 = v6-A;

      // Calculate barycentric coords
      vnl_matrix_fixed<double,3,2> mat;
      mat.set_column(0, v0);
      mat.set_column(1, v1);
      vnl_matrix_inverse<double> inv(mat);
      vnl_matrix_fixed<double,2,3> inver = inv.pinverse();
      vnl_vector<double> uv = inv.pinverse()*v6;

      // Check if point is in triangle
      double eps = 0.01;
      if( (uv(0) >= 0-eps) && (uv(1) >= 0-eps) && (uv(0) + uv(1) <= 1+eps) )
      {
        // check if minimum angle is the max so far
        if(d1.two_norm() < current_min)
        {
          current_min = d1.two_norm();
          vnl_vector<double> barycentricCoords(3);
          barycentricCoords[2] = uv[0]<0 ? 0 : (uv[0]>1?1:uv[0]);
          barycentricCoords[1] = uv[1]<0 ? 0 : (uv[1]>1?1:uv[1]);
          barycentricCoords[0] = 1-(barycentricCoords[1]+barycentricCoords[2]);
          retval =  barycentricCoords[0]*this->GetNthComponent(pts[0]) +
              barycentricCoords[1]*this->GetNthComponent(pts[1]) +
              barycentricCoords[2]*this->GetNthComponent(pts[2]);
        }
      }
    }

    break;
  }
  case ODF_SPHERICAL_GAUSSIAN_BASIS_FUNCTIONS:
  {
    double maxChordLength = GetMaxChordLength();
    double sigma = asin(maxChordLength/2);

    // this is the contribution of each kernel to each sampling point on the
    // equator
    vnl_vector<double> contrib;
    contrib.set_size(NOdfDirections);

    vtkPoints* points = m_BaseMesh->GetPoints();
    double sum = 0;
    for(int i=0; i<NOdfDirections; i++)
    {
      vnl_vector_fixed<double,3> P(points->GetPoint(i));
      double stv =  dir[0]*P[0]
                    + dir[1]*P[1]
                    + dir[2]*P[2];
      stv = (stv<-1.0) ? -1.0 : ( (stv>1.0) ? 1.0 : stv);
      double x = acos(stv);
      contrib[i] = (1.0/(sigma*sqrt(2.0*itk::Math::pi)))
                   *exp((-x*x)/(2*sigma*sigma));
      sum += contrib[i];
    }

    retval = 0;
    for(int i=0; i<NOdfDirections; i++)
    {
      retval += (contrib[i] / sum)*this->GetNthComponent(i);
    }
    break;
  }

  }

  if(retval==-1)
  {
    std::cout << "Interpolation failed" << std::endl;
    return 0;
  }

  return retval;

}

/**
  * Calculate Generalized Fractional Anisotropy
  */
template<class T, unsigned int NOdfDirections>
T
OrientationDistributionFunction<T, NOdfDirections>
::GetGeneralizedFractionalAnisotropy() const
{
  double mean = 0;
  double std = 0;
  double rms = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    T val = (*this)[i];
    mean += val;
  }
  mean /= NOdfDirections;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    T val = (*this)[i];
    std += (val - mean) * (val - mean);
    rms += val*val;
  }
  std *= NOdfDirections;
  rms *= NOdfDirections - 1;

  if(rms == 0)
  {
    return 0;
  }
  else
  {
    return sqrt(std/rms);
  }
}


template < typename T, unsigned int N>
T itk::OrientationDistributionFunction<T, N>
::GetGeneralizedGFA( int k, int p ) const
{
  double mean = 0;
  double std = 0;
  double rms = 0;
  double max = NumericTraits<double>::NonpositiveMin();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    double val = (double)(*this)[i];
    mean += pow(val,(double)p);
    max = val > max ? val : max;
  }
  max = pow(max,(double)p);
  mean /= N;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    double val = (double)(*this)[i];
    std += (pow(val,(double)p) - mean) * (pow(val,(double)p) - mean);
    if(k>0)
    {
      rms += pow(val,(double)(p*k));
    }
  }
  std /= N - 1;
  std = sqrt(std);

  if(k>0)
  {
    rms /= N;
    rms = pow(rms,(double)(1.0/k));
  }
  else if(k<0) // lim k->inf gives us the maximum
  {
    rms = max;
  }
  else // k==0 undefined, we define zeros root from 1 as 1
  {
    rms = 1;
  }

  if(rms == 0)
  {
    return 0;
  }
  else
  {
    return (T)(std/rms);
  }
}

/**
  * Calculate Nematic Order Parameter
  */
template < typename T, unsigned int N >
T itk::OrientationDistributionFunction<T, N>
::GetNematicOrderParameter() const
{
  // not yet implemented
  return 0;
}

/**
  * Calculate StdDev by MaxValue
  */
template < typename T, unsigned int N >
T itk::OrientationDistributionFunction<T, N>
::GetStdDevByMaxValue() const
{
  double mean = 0;
  double std = 0;
  T max = NumericTraits<T>::NonpositiveMin();
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    T val = (*this)[i];
    mean += val;
    max = (*this)[i] > max ? (*this)[i] : max;
  }
  mean /= InternalDimension;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    T val = (*this)[i];
    std += (val - mean) * (val - mean);
  }
  std /= InternalDimension-1;

  if(max == 0)
  {
    return 0;
  }
  else
  {
    return (sqrt(std)/max);
  }
}

template < typename T, unsigned int N >
T itk::OrientationDistributionFunction<T, N>
::GetPrincipleCurvature(double alphaMinDegree, double alphaMaxDegree, int invert) const
{
  // following loop only performed once
  // (computing indices of each angular range)
  m_MutexAngularRange.lock();
  if(m_AngularRangeIdxs == nullptr)
  {
    m_AngularRangeIdxs = new std::vector< std::vector<int>* >();
    for(unsigned int i=0; i<N; i++)
    {
      vnl_vector_fixed<double,3> pDir = GetDirection(i);
      auto      idxs = new std::vector<int>();
      for(unsigned int j=0; j<N; j++)
      {
        vnl_vector_fixed<double,3> cDir = GetDirection(j);
        double angle = ( 180 / itk::Math::pi ) * acos( dot_product(pDir, cDir) );
        if( (angle < alphaMaxDegree) && (angle > alphaMinDegree) )
        {
          idxs->push_back(j);
        }
      }
      m_AngularRangeIdxs->push_back(idxs);
    }
  }
  m_MutexAngularRange.unlock();

  // find the maximum (or minimum) direction (remember index and value)
  T mode;
  int pIdx = -1;
  if(invert == 0)
  {
    pIdx = GetPrincipalDiffusionDirectionIndex();
    mode = (*this)[pIdx];
  }
  else
  {
    mode = NumericTraits<T>::max();
    for( unsigned int i=0; i<N; i++)
    {
      if((*this)[i] < mode )
      {
        mode = (*this)[i];
        pIdx = i;
      }
    }
  }
  //////////////////////
  //////// compute median of mode and its neighbors to become more stable to noise
  //////// compared to simply using the mode
  //////////////////////

  //////// values of mode and its neighbors
  //////std::vector<int> nbs = GetNeighbors(pIdx);
  //////std::vector<T> modeAndNeighborVals;
  //////modeAndNeighborVals.push_back(mode);
  //////int numNeighbors = nbs.size();
  //////for(int i=0; i<numNeighbors; i++)
  //////{
  //////  modeAndNeighborVals.push_back((*this)[nbs[i]]);
  //////}

  //////// sort by value
  //////std::sort( modeAndNeighborVals.begin(), modeAndNeighborVals.end() );

  //////// median of mode and neighbors
  //////mode = modeAndNeighborVals[floor(0.5*(double)(numNeighbors+1)+0.5)];

  ////////////////
  // computing a quantile of the angular range
  ////////////////

  // define quantile
  double quantile = 0.00;

  // collect all values in angular range of mode
  std::vector<T> odfValuesInAngularRange;
  int numInRange = m_AngularRangeIdxs->at(pIdx)->size();
  for(int i=0; i<numInRange; i++)
  {
    odfValuesInAngularRange.push_back((*this)[(*m_AngularRangeIdxs->at(pIdx))[i] ]);
  }

  // sort them by value
  std::sort( odfValuesInAngularRange.begin(), odfValuesInAngularRange.end() );

  // median of angular range
  T median = odfValuesInAngularRange[floor(quantile*(double)numInRange+0.5)];

  // compute and return final value
  if(mode > median)
  {
    return mode/median - 1.0;
  }
  else
  {
    return median/mode - 1.0;
  }
}

/**
  * Calculate Normalized Entropy
  */
template < typename T, unsigned int N >
T itk::OrientationDistributionFunction<T, N>
::GetNormalizedEntropy() const
{
  double mean = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
  {
    T val = (*this)[i];
    if( val != 0 )
    {
      val = log(val);
    }
    else
    {
      val = log(0.0000001);
    }
    mean += val;
  }
  double _n = (double) InternalDimension;
  mean /= _n;
  return (T) (-_n / log(_n) * mean);
}

/**
  * Pre-multiply the Tensor by a Matrix
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::PreMultiply( const MatrixType & m ) const
{
  Self result;
  typedef typename NumericTraits<T>::AccumulateType  AccumulateType;
  for(unsigned int r=0; r<NOdfDirections; r++)
  {
    for(unsigned int c=0; c<NOdfDirections; c++)
    {
      AccumulateType sum = NumericTraits<AccumulateType>::ZeroValue();
      for(unsigned int t=0; t<NOdfDirections; t++)
      {
        sum += m(r,t) * (*this)(t,c);
      }
      result(r,c) = static_cast<T>( sum );
    }
  }
  return result;
}

/**
  * Post-multiply the Tensor by a Matrix
  */
template<class T, unsigned int NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
OrientationDistributionFunction<T, NOdfDirections>
::PostMultiply( const MatrixType & m ) const
{
  Self result;
  typedef typename NumericTraits<T>::AccumulateType  AccumulateType;
  for(unsigned int r=0; r<NOdfDirections; r++)
  {
    for(unsigned int c=0; c<NOdfDirections; c++)
    {
      AccumulateType sum = NumericTraits<AccumulateType>::ZeroValue();
      for(unsigned int t=0; t<NOdfDirections; t++)
      {
        sum += (*this)(r,t) * m(t,c);
      }
      result(r,c) = static_cast<T>( sum );
    }
  }
  return result;
}

/**
  * Print content to an ostream
  */
template<class T, unsigned int NOdfDirections>
std::ostream &
operator<<(std::ostream& os,const OrientationDistributionFunction<T, NOdfDirections> & c )
{
  for(unsigned int i=0; i<c.GetNumberOfComponents(); i++)
  {
    os <<  static_cast<typename NumericTraits<T>::PrintType>(c[i]) << "  ";
  }
  return os;
}

/**
  * Read content from an istream
  */
template<class T, unsigned int NOdfDirections>
std::istream &
operator>>(std::istream& is, OrientationDistributionFunction<T, NOdfDirections> & dt )
{
  for(unsigned int i=0; i < dt.GetNumberOfComponents(); i++)
  {
    is >> dt[i];
  }
  return is;
}
} // end namespace itk

#endif
