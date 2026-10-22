/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkB0Indices_h
#define mitkB0Indices_h

#include <mitkBValueMapProperty.h>

#include <vector>

namespace mitk
{
  /**
   * \brief Whether a volume is a b0 volume: its b-value, reference b-value times the squared
   * gradient length, rounds to 0.
   *
   * This is the definition of the b-value map (BValueMapProperty::GetBValue), so the filters of
   * this module agree with the shells the views and the property helper show.
   */
  inline bool IsB0(const BValueMapProperty::GradientDirectionType& gradient, double referenceBValue)
  {
    return BValueMapProperty::GetBValue(gradient, referenceBValue) == 0;
  }

  /** \brief Indices of the b0 volumes (see IsB0), in ascending order. */
  inline std::vector<unsigned int> GetB0Indices(const BValueMapProperty::GradientDirectionsContainerType* gradients, double referenceBValue)
  {
    std::vector<unsigned int> indices;
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      if (IsB0(gradients->ElementAt(i), referenceBValue))
        indices.push_back(i);
    return indices;
  }
}

#endif
