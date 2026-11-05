/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkShEvaluator.h>

#include <itkMath.h>
#include <mitkExceptionMacro.h>

#include <vnl/vnl_cross.h>

#include <algorithm>
#include <cmath>

mitk::ShEvaluator::ShEvaluator(unsigned int order, bool mrtrix)
  : m_Order(order)
{
  if (order % 2 != 0 || order > MaxOrder)
    mitkThrow() << "SH order " << order << " is not an even order up to " << MaxOrder << ".";

  // The basis of mitk::sh::Yj, written with Pbar_k^a, the associated Legendre function without
  // the Condon-Shortley phase (Boost's legendre_p and spherical_harmonic include it), and
  // N = sqrt((2k+1)/(4 pi) (k-a)!/(k+a)!):
  //   MRtrix (legacy): legendre_p(k, a, -cos(theta)) = Pbar_k^a(cos(theta)) for even k, so
  //     +a: N Pbar cos(a phi), -a: N Pbar sin(a phi);
  //   FSL: spherical_harmonic(k, a) = (-1)^a N Pbar e^{i a phi}, so
  //     +a: sqrt(2) N Pbar sin(a phi), -a: (-1)^a sqrt(2) N Pbar cos(a phi).
  // Pbar_k^a(cos(theta)) cos(a phi) = R(z) Re((x + iy)^a), and the same with sin and Im.
  const double sqrt2 = std::sqrt(2.0);
  for (unsigned int k = 0; k <= order; k += 2)
  {
    const unsigned int centre = (k * k + k + 2) / 2 - 1;
    for (unsigned int a = 0; a <= k; ++a)
    {
      double factorialRatio = 1.0; // (k-a)! / (k+a)!
      for (unsigned int j = k - a + 1; j <= k + a; ++j)
        factorialRatio /= j;
      const double n = std::sqrt((2.0 * k + 1.0) / (4.0 * itk::Math::pi) * factorialRatio);

      Term term;
      term.k = k;
      term.a = a;
      if (a == 0)
      {
        term.realIndex = centre;
        term.realFactor = n;
        term.imagIndex = centre;
        term.imagFactor = 0.0;
      }
      else if (mrtrix)
      {
        term.realIndex = centre + a;
        term.realFactor = n;
        term.imagIndex = centre - a;
        term.imagFactor = n;
      }
      else
      {
        term.realIndex = centre - a;
        term.realFactor = (a % 2 == 0 ? 1.0 : -1.0) * sqrt2 * n;
        term.imagIndex = centre + a;
        term.imagFactor = sqrt2 * n;
      }
      m_Terms.push_back(term);
    }
  }
}

template <bool Derivatives>
double mitk::ShEvaluator::Evaluate(const double *coefficients, const Vector3 &direction, double *derivatives) const
{
  const unsigned int order = m_Order;
  const double x = direction[0];
  const double y = direction[1];
  const double z = direction[2];

  // r[l][m] = d^m P_l / dz^m, so that d r[l][m] / dz = r[l][m+1]; the entries m = l+1 and l+2
  // are 0. Upward recurrence in l from r[m][m] = (2m-1)!!.
  double r[MaxOrder + 1][MaxOrder + 3];
  double doubleFactorial = 1.0;
  for (unsigned int m = 0; m <= order; ++m)
  {
    if (m > 0)
      doubleFactorial *= 2.0 * m - 1.0;
    r[m][m] = doubleFactorial;
    if (m + 1 <= order)
      r[m + 1][m] = z * (2.0 * m + 1.0) * doubleFactorial;
    for (unsigned int l = m + 2; l <= order; ++l)
      r[l][m] = (z * (2.0 * l - 1.0) * r[l - 1][m] - (l + m - 1.0) * r[l - 2][m]) / (l - m);
  }
  for (unsigned int l = 0; l <= order; ++l)
  {
    r[l][l + 1] = 0.0;
    r[l][l + 2] = 0.0;
  }

  // Real and imaginary part of (x + iy)^a.
  double re[MaxOrder + 1];
  double im[MaxOrder + 1];
  re[0] = 1.0;
  im[0] = 0.0;
  for (unsigned int a = 1; a <= order; ++a)
  {
    re[a] = re[a - 1] * x - im[a - 1] * y;
    im[a] = re[a - 1] * y + im[a - 1] * x;
  }

  double f = 0.0;
  double fx = 0.0, fy = 0.0, fz = 0.0;
  double fxx = 0.0, fxy = 0.0, fxz = 0.0, fyy = 0.0, fyz = 0.0, fzz = 0.0;
  for (const Term &term : m_Terms)
  {
    const unsigned int a = term.a;
    const double cr = coefficients[term.realIndex] * term.realFactor;
    const double ci = coefficients[term.imagIndex] * term.imagFactor;
    const double t = cr * re[a] + ci * im[a];
    const double rk = r[term.k][a];
    f += rk * t;
    if constexpr (Derivatives)
    {
      // d/dx (x + iy)^a = a (x + iy)^(a-1), d/dy = i a (x + iy)^(a-1).
      double tx = 0.0, ty = 0.0, txx = 0.0, txy = 0.0;
      if (a >= 1)
      {
        tx = a * (cr * re[a - 1] + ci * im[a - 1]);
        ty = a * (ci * re[a - 1] - cr * im[a - 1]);
      }
      if (a >= 2)
      {
        txx = a * (a - 1.0) * (cr * re[a - 2] + ci * im[a - 2]);
        txy = a * (a - 1.0) * (ci * re[a - 2] - cr * im[a - 2]);
      }
      const double rz = r[term.k][a + 1];
      fx += rk * tx;
      fy += rk * ty;
      fz += rz * t;
      fxx += rk * txx;
      fxy += rk * txy;
      fyy -= rk * txx;  // the terms are harmonic in x and y
      fxz += rz * tx;
      fyz += rz * ty;
      fzz += r[term.k][a + 2] * t;
    }
  }
  if constexpr (Derivatives)
  {
    derivatives[0] = fx;
    derivatives[1] = fy;
    derivatives[2] = fz;
    derivatives[3] = fxx;
    derivatives[4] = fxy;
    derivatives[5] = fxz;
    derivatives[6] = fyy;
    derivatives[7] = fyz;
    derivatives[8] = fzz;
  }
  return f;
}

double mitk::ShEvaluator::GetValue(const double *coefficients, const Vector3 &direction) const
{
  return Evaluate<false>(coefficients, direction, nullptr);
}

double mitk::ShEvaluator::GetValue(const double *coefficients, const Vector3 &direction, Vector3 &gradient, Matrix3 &hessian) const
{
  double d[9];
  const double value = Evaluate<true>(coefficients, direction, d);

  // The polynomial is defined in all of R^3; on the sphere its gradient is the tangential part,
  // and its Riemannian Hessian P H P - (u . g) P with the projection P = I - u u^T.
  const Vector3 g(d[0], d[1], d[2]);
  Matrix3 h;
  h(0, 0) = d[3]; h(0, 1) = d[4]; h(0, 2) = d[5];
  h(1, 0) = d[4]; h(1, 1) = d[6]; h(1, 2) = d[7];
  h(2, 0) = d[5]; h(2, 1) = d[7]; h(2, 2) = d[8];

  const double radial = dot_product(direction, g);
  Matrix3 projection;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      projection(i, j) = (i == j ? 1.0 : 0.0) - direction[i] * direction[j];

  gradient = g - radial * direction;
  hessian = projection * h * projection - radial * projection;
  return value;
}

double mitk::ShEvaluator::FindLocalMaximum(const double *coefficients, Vector3 &direction, unsigned int maxIterations, unsigned int *iterations) const
{
  const unsigned int limit = maxIterations > 0 ? maxIterations : 1000;
  Vector3 u = direction / direction.magnitude();
  Vector3 gradient;
  Matrix3 hessian;
  double value = GetValue(coefficients, u, gradient, hessian);

  unsigned int steps = 0;
  for (; steps < limit; ++steps)
  {
    const double gradientNorm = gradient.magnitude();
    if (gradientNorm <= 1e-9 * std::fabs(value))
      break;

    // Orthonormal basis of the tangent plane; gradient and Hessian in it.
    Vector3 e1 = vnl_cross_3d(u, std::fabs(u[0]) < 0.6 ? Vector3(1, 0, 0) : Vector3(0, 1, 0));
    e1.normalize();
    const Vector3 e2 = vnl_cross_3d(u, e1);
    const double g1 = dot_product(e1, gradient);
    const double g2 = dot_product(e2, gradient);
    const Vector3 he1 = hessian * e1;
    const Vector3 he2 = hessian * e2;
    const double h11 = dot_product(e1, he1);
    const double h12 = dot_product(e1, he2);
    const double h22 = dot_product(e2, he2);

    // Newton step with the eigenvalues of the 2x2 Hessian replaced by -max(|lambda|, floor):
    // the Newton step where the Hessian is negative definite; elsewhere (near a saddle or on a
    // flat ridge) it still moves uphill along each eigenvector, by the gradient over the
    // curvature's magnitude.
    const double mean = 0.5 * (h11 + h22);
    const double radius = std::hypot(0.5 * (h11 - h22), h12);
    const double lambdas[2] = {mean - radius, mean + radius};
    double eigen[2][2];
    if (radius <= 1e-14 * std::fabs(mean))
    {
      eigen[0][0] = 1; eigen[0][1] = 0;
      eigen[1][0] = 0; eigen[1][1] = 1;
    }
    else
    {
      // Eigenvector of lambdas[0]: (h12, lambda - h11) or (lambda - h22, h12), the longer one.
      double a1 = h12, a2 = lambdas[0] - h11;
      const double b1 = lambdas[0] - h22, b2 = h12;
      if (a1 * a1 + a2 * a2 < b1 * b1 + b2 * b2)
      {
        a1 = b1;
        a2 = b2;
      }
      const double n = std::hypot(a1, a2);
      eigen[0][0] = a1 / n; eigen[0][1] = a2 / n;
      eigen[1][0] = -a2 / n; eigen[1][1] = a1 / n;
    }
    const double curvatureFloor = 1e-3 * std::fabs(value) + 1e-300;
    double s1 = 0, s2 = 0;
    for (unsigned int i = 0; i < 2; ++i)
    {
      const double along = (eigen[i][0] * g1 + eigen[i][1] * g2) / std::max(std::fabs(lambdas[i]), curvatureFloor);
      s1 += along * eigen[i][0];
      s2 += along * eigen[i][1];
    }

    double length = std::min(std::hypot(s1, s2), MaxStep);
    const Vector3 tangent = (s1 * e1 + s2 * e2).normalize();
    bool accepted = false;
    Vector3 newU, newGradient;
    Matrix3 newHessian;
    double newValue = value;
    for (; length >= 1e-8; length *= 0.5)
    {
      // Along the great circle through u in the step's direction.
      newU = std::cos(length) * u + std::sin(length) * tangent;
      newU.normalize();
      newValue = GetValue(coefficients, newU, newGradient, newHessian);
      if (newValue > value)
      {
        accepted = true;
        break;
      }
    }
    if (!accepted)
      break;
    u = newU;
    value = newValue;
    gradient = newGradient;
    hessian = newHessian;
  }

  direction = u;
  if (iterations != nullptr)
    *iterations = steps;
  return value;
}
