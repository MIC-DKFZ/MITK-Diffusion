/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkShRegularization_h
#define mitkShRegularization_h

namespace mitk
{
  /**
   * \brief Diagonal element of the Laplace-Beltrami regularization for the coefficient (l, m) of
   * MITK's SH basis (mitk::sh::Yj with mrtrix = true), as the Q-ball filters fit it.
   *
   * The regularized fit (Descoteaux et al. 2007) minimizes |B c - s|^2 + lambda * R(c), with R the
   * integral of the squared Laplace-Beltrami operator of the fitted function over the sphere:
   * R(c) = sum_j (l_j (l_j + 1))^2 |Y_j|^2 c_j^2. MITK's basis is MRtrix 0.2's, whose functions with
   * m != 0 lack the factor sqrt(2) of an orthonormal basis and have the squared norm |Y_j|^2 = 1/2.
   * Only with that weight is the fit rotation invariant; it then equals a fit in an orthonormal basis
   * (as dipy's QballModel and CsaOdfModel do it) converted to MITK's basis.
   */
  inline double LaplaceBeltramiPenalty(int l, int m)
  {
    const double eigenvalue = l * (l + 1);
    return eigenvalue * eigenvalue * (m == 0 ? 1.0 : 0.5);
  }
}

#endif
