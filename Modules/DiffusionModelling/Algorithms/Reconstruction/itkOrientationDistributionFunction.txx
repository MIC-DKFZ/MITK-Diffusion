/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkOrientationDistributionFunction_txx
#define itkOrientationDistributionFunction_txx

#include <vnl/algo/vnl_matrix_inverse.h>
#include <itkMatrix.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_vector_fixed.hxx>
#include <vnl/vnl_cross.h>
#include <mitkLog.h>
#include <vtkCellArrayIterator.h>
#include <vtkDelaunay2D.h>
#include <vtkNew.h>

#include <algorithm>
#include <array>
#include <functional>
#include <set>

namespace itk
{

template<class T, unsigned int N>
std::map< std::pair<double, double>, std::vector< std::vector<int> > > itk::OrientationDistributionFunction<T,N>::m_AngularRangeIdxs;

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
  CopyEigenAnalysis(r);
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
  // Fill, not BaseArray::operator=(&r), which would read N values from the address of the one
  // scalar.
  this->Fill(r);
  m_EigenAnalysisCalculated = false;
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
  m_EigenAnalysisCalculated = false;
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
  * InitFromTensor()
  */
template<class T, unsigned int NOdfDirections>
void
OrientationDistributionFunction<T, NOdfDirections>
::InitFromTensor(itk::DiffusionTensor3D<T> tensor)
{
  m_EigenAnalysisCalculated = false;
  const DirectionsType &directions = GetDirections();
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
    T g0 = directions(0,i);
    T g1 = directions(1,i);
    T g2 = directions(2,i);
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
  int exponent_a = std::ilogb(a);
  int exponent_b = std::ilogb(b);
  int exponent_c = std::ilogb(c);

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

  const DirectionsType &directions = GetDirections();
  for( unsigned int i=0; i < NOdfDirections; ++i )
  {
    /// calculate probability p(g)=r as ellipsoids magnitude in direction of g' = B^-1 * g:
    /// (g0*r/evx)² + (g1*r/evy)² + (g2*r/evz)² = 1, |g'|=1.

    vnl_vector_fixed<double, 3> g( directions(0,i), directions(1,i), directions(2,i) );
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
    // The radius lies between the smallest and the largest eigenvalue, which exceed 1 for tensors
    // not given in mm^2/s; only negative and NaN radii are invalid.
    if ( (*this)[i] < T(0) || std::isnan((*this)[i]) )
    {
      MITK_DEBUG << "OrientationDistributionFunction<" << typeid(T).name() << ", " << NOdfDirections
                 << ">::InitFromEllipsoid(" << typeid(tensor).name()
                 << ") encountered a negative or NaN radius and set it to zero: (*this)["
                 << i <<"]= " << (*this)[i];
      (*this)[i] = T(0);
    }
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
  // A constant ODF (e.g. all zero) has no range; it stays zero instead of becoming NaN.
  if (!(max > min))
    return retval;
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
  // A zero ODF stays zero instead of becoming NaN.
  if (max == 0)
    return retval;
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
::GetMeanValue() const
{
  T sum = 0;
  for( unsigned int i=0; i<InternalDimension; i++)
    sum += (*this)[i];
  return sum/InternalDimension;
}

template<class T, unsigned int NOdfDirections>
vtkSmartPointer<vtkPolyData>
OrientationDistributionFunction<T, NOdfDirections>
::BuildBaseMesh()
{
  // Two Delaunay triangulations in azimuth/elevation, of the directions and of the directions
  // rotated so that the first one's poles and seam lie in the middle; each keeps the triangles
  // away from its poles and seam. Both keep the triangles of the band between, which the mesh holds
  // once: each glyph draws every triangle, and repeated ones would weigh twice in its normals.
  vtkNew<vtkCellArray> polys;
  std::set<std::array<vtkIdType, 3>> inserted;
  const auto triangulate = [&polys, &inserted](const std::function<vnl_vector_fixed<double, 3>(unsigned int)> &direction) {
    vtkNew<vtkPoints> points;
    for (unsigned int j = 0; j < NOdfDirections; j++)
    {
      const auto d = direction(j);
      points->InsertNextPoint(atan2(d[1], d[0]), atan2(d[2], sqrt(d[0] * d[0] + d[1] * d[1])), d.two_norm());
    }
    vtkNew<vtkPolyData> projected;
    projected->SetPoints(points);
    vtkNew<vtkDelaunay2D> delaunay;
    delaunay->SetInputData(projected);
    delaunay->Update();

    vtkCellArray *triangles = delaunay->GetOutput()->GetPolys();
    vtkIdType npts;
    vtkIdType const *pts;
    for (triangles->InitTraversal(); triangles->GetNextCell(npts, pts);)
    {
      bool insert = true;
      for (int i = 0; i < npts; i++)
      {
        const double *point = points->GetPoint(pts[i]);
        if ((std::abs(point[0]) > itk::Math::pi - 0.5) || (std::abs(point[1]) > itk::Math::pi / 2 - 0.5))
          insert = false;
      }
      std::array<vtkIdType, 3> triangle = {pts[0], pts[1], pts[2]};
      std::sort(triangle.begin(), triangle.end());
      if (insert && inserted.insert(triangle).second)
        polys->InsertNextCell(npts, pts);
    }
  };
  const auto &directions = GetDirections();
  triangulate([&directions](unsigned int j) { return vnl_vector_fixed<double, 3>(directions(0, j), directions(1, j), directions(2, j)); });
  triangulate([&directions](unsigned int j) { return vnl_vector_fixed<double, 3>(-directions(0, j), -directions(2, j), -directions(1, j)); });

  vtkNew<vtkPoints> points;
  for (unsigned int p = 0; p < NOdfDirections; p++)
    points->InsertNextPoint(directions.get_column(p).data_block());
  auto mesh = vtkSmartPointer<vtkPolyData>::New();
  mesh->SetPoints(points);
  mesh->SetPolys(polys);
  return mesh;
}

template<class T, unsigned int NOdfDirections>
vtkPolyData *
OrientationDistributionFunction<T, NOdfDirections>
::GetBaseMesh()
{
  // Built once, on first use from any thread, and never changed.
  static const vtkSmartPointer<vtkPolyData> mesh = BuildBaseMesh();
  return mesh;
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
    if (idx>=0 && idx<(int)NOdfDirections)
      return OrientationDistributionFunction<T, NOdfDirections>::GetDirection(idx);
    vnl_vector_fixed<double,3> vec; vec.fill(0);
    return vec;
  }
}

template<class T, unsigned int NOdfDirections>
const std::vector<int> &
OrientationDistributionFunction<T, NOdfDirections>
::GetNeighbors(int idx)
{
  // Built once, on first use from any thread, and never changed: for each direction the sorted
  // indices of the directions it shares a base mesh triangle with.
  static const std::vector<std::vector<int>> neighbors = []() {
    std::vector<std::vector<int>> result(NOdfDirections);
    // An iterator of its own: the traversal position of the shared cell array belongs to its users.
    auto triangles = vtk::TakeSmartPointer(GetBaseMesh()->GetPolys()->NewIterator());
    for (triangles->GoToFirstCell(); !triangles->IsDoneWithTraversal(); triangles->GoToNextCell())
    {
      vtkIdType npts;
      vtkIdType const *pts;
      triangles->GetCurrentCell(npts, pts);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          if (a != b)
            result[pts[a]].push_back(pts[b]);
    }
    for (auto &idxs : result)
    {
      std::sort(idxs.begin(), idxs.end());
      idxs.erase(std::unique(idxs.begin(), idxs.end()), idxs.end());
    }
    return result;
  }();
  return neighbors.at(idx);
}

template < typename TComponent, unsigned int NOdfDirections >
vnl_vector_fixed<double,3> itk::OrientationDistributionFunction<TComponent, NOdfDirections>
::GetDirection( int i )
{
  return GetDirections().get_column(i);
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
  // The angular ranges depend on the two angles (node properties of the glyphs, parameters of the
  // GFA filter), so they are computed once per pair of angles. std::map keeps the references valid
  // when other pairs are added.
  const std::vector< std::vector<int> > *angularRanges = nullptr;
  {
    std::lock_guard<std::mutex> lock(m_MutexAngularRange);
    auto &ranges = m_AngularRangeIdxs[std::make_pair(alphaMinDegree, alphaMaxDegree)];
    if (ranges.empty())
    {
      ranges.resize(N);
      for(unsigned int i=0; i<N; i++)
      {
        vnl_vector_fixed<double,3> pDir = GetDirection(i);
        for(unsigned int j=0; j<N; j++)
        {
          vnl_vector_fixed<double,3> cDir = GetDirection(j);
          double angle = ( 180 / itk::Math::pi ) * acos( dot_product(pDir, cDir) );
          if( (angle < alphaMaxDegree) && (angle > alphaMinDegree) )
          {
            ranges[i].push_back(j);
          }
        }
      }
    }
    angularRanges = &ranges;
  }

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
  int numInRange = angularRanges->at(pIdx).size();
  // An empty angular range or an ODF without a positive mode or median has no curvature.
  if (numInRange == 0)
    return 0;
  for(int i=0; i<numInRange; i++)
  {
    odfValuesInAngularRange.push_back((*this)[angularRanges->at(pIdx)[i] ]);
  }

  // sort them by value
  std::sort( odfValuesInAngularRange.begin(), odfValuesInAngularRange.end() );

  // median of angular range
  T median = odfValuesInAngularRange[floor(quantile*(double)numInRange+0.5)];

  // compute and return final value
  if (!(mode > 0) || !(median > 0))
    return 0;
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
