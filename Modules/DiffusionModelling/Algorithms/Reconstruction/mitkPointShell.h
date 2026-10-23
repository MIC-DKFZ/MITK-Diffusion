/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPointShell_h
#define mitkPointShell_h

#include <MitkDiffusionModellingExports.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
// The member definitions of vnl_matrix_fixed, which ITK instantiates only for small sizes.
#include <vnl/vnl_matrix_fixed.hxx>

#include <vector>

namespace mitk::pointshell
{
  /** \brief The numbers of directions with a table: 12, 42, 92, 162, 252, 362, 492, 642, 812 and 1002. */
  MITKDIFFUSIONMODELLING_EXPORT const std::vector<unsigned int> &SupportedSizes();

  /**
   * \brief The n sampling directions of the sphere, as the columns of a 3 x n matrix.
   *
   * The tables are the icosahedron (12) and its subdivisions, unit vectors in antipodal pairs (to
   * about 1.5e-6). Their values and their order are part of the ODF image format (one value per
   * direction) and of Fiberfox's stick directions, so they are stored, not computed. Each table is
   * built once and shared.
   *
   * \throw mitk::Exception for n without a table.
   */
  MITKDIFFUSIONMODELLING_EXPORT const vnl_matrix<double> &Directions(unsigned int n);

  /** \brief Directions(N) as fixed-size matrix, built once. */
  template <unsigned int N>
  const vnl_matrix_fixed<double, 3, N> &Directions()
  {
    static const vnl_matrix_fixed<double, 3, N> directions(Directions(N).data_block());
    return directions;
  }
}

#endif
