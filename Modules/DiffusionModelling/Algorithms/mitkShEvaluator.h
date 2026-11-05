/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkShEvaluator_h
#define mitkShEvaluator_h

#include <MitkDiffusionModellingExports.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

#include <vector>

namespace mitk
{
  /**
   * Evaluates a real SH function of even order in double, with its gradient and Hessian on the
   * unit sphere, in the basis of mitk::sh::Yj (MITK's legacy "MRtrix" basis or MITK's FSL
   * definition).
   *
   * Each basis function of order k and |m| = a is N * R(z) * Re or Im of (x + iy)^a, where R is
   * the a-th derivative of the Legendre polynomial P_k and N the normalization; this is a
   * polynomial in the direction's x, y and z, so the derivatives have no singularity at the poles.
   * The normalization constants are computed once in the constructor, the Legendre derivatives
   * per evaluation by the three-term recurrence, without Boost calls and without allocation.
   * The methods are const and keep no state, so one evaluator serves several threads.
   */
  class MITKDIFFUSIONMODELLING_EXPORT ShEvaluator
  {
  public:
    using Vector3 = vnl_vector_fixed<double, 3>;
    using Matrix3 = vnl_matrix_fixed<double, 3, 3>;

    /** Highest order supported; the evaluation works on arrays of this size on the stack. */
    static constexpr unsigned int MaxOrder = 20;

    /** Throws for an odd order or one above MaxOrder. */
    ShEvaluator(unsigned int order, bool mrtrix);

    unsigned int GetOrder() const { return m_Order; }
    unsigned int GetNumberOfCoefficients() const { return (m_Order + 1) * (m_Order + 2) / 2; }

    /** Value at the unit vector direction, with GetNumberOfCoefficients() coefficients. */
    double GetValue(const double *coefficients, const Vector3 &direction) const;

    /**
     * Value at the unit vector direction, its gradient on the sphere (a tangent vector) and its
     * Riemannian Hessian on the sphere: for tangent vectors t, t^T hessian t is the second
     * derivative of the value along the great circle through direction in the direction of t.
     */
    double GetValue(const double *coefficients, const Vector3 &direction, Vector3 &gradient, Matrix3 &hessian) const;

    /**
     * Moves the unit vector direction uphill on the sphere to a local maximum of the SH function
     * and returns the value there. Each step is a Newton step in the tangent plane with the
     * eigenvalues of the Hessian replaced by -max(|eigenvalue|, 1e-3 |value|), so it is the Newton
     * step where the Hessian is negative definite and still goes uphill near saddles and on flat
     * ridges; it is at most MaxStep long and halved until the value increases, so the result is
     * never below the value at the start. Stops when the gradient is
     * at most 1e-9 of the value, when no step of at least 1e-8 rad increases the value, or after
     * maxIterations steps (0: at most 1000). iterations, if given, gets the number of steps taken.
     */
    double FindLocalMaximum(const double *coefficients, Vector3 &direction, unsigned int maxIterations, unsigned int *iterations = nullptr) const;

    /** Longest step of FindLocalMaximum in rad (about 3 degrees, a quarter of the spacing of the
      * 252 ODF directions). A step is accepted when the value at its end is higher, so a longer
      * one can cross a saddle into the basin of another maximum: with 0.2 rad, the candidate of a
      * flat secondary peak of an order-4 ODF ended on the primary peak 70 degrees away. */
    static constexpr double MaxStep = 0.05;

  private:
    /** Coefficient index and factor (normalization and convention sign) of the Re and Im part of
      * order k and |m| = a. */
    struct Term
    {
      unsigned int k;
      unsigned int a;
      unsigned int realIndex;
      double realFactor;
      unsigned int imagIndex;  ///< unused for a = 0
      double imagFactor;       ///< 0 for a = 0
    };

    template <bool Derivatives>
    double Evaluate(const double *coefficients, const Vector3 &direction, double *derivatives) const;

    unsigned int m_Order;
    std::vector<Term> m_Terms;
  };
}

#endif
