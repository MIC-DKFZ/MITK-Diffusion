/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkSingleShellGradients_h
#define mitkSingleShellGradients_h

#include <mitkBValueMapProperty.h>
#include <mitkDiffusionPropertyHelper.h>

#include <itkMacro.h>

#include <set>
#include <sstream>
#include <vector>

namespace mitk
{
  /** \brief The volumes of a single-shell DWI, as the single-shell Q-ball filters use them. */
  struct SingleShellGradients
  {
    std::vector<unsigned int> baselineIndices;
    std::vector<unsigned int> diffusionWeightedIndices;
    /** Unit directions of the diffusion-weighted volumes, in the order of diffusionWeightedIndices. */
    std::vector<vnl_vector_fixed<double, 3>> directions;
    /** b-value of the diffusion-weighted volumes as the b-value map gives it (one shell, so all share
     * it); below the reference b-value when e.g. other shells were removed. */
    double bValue = 0;
  };

  /**
   * \brief Splits the gradients of a DWI into b0 and diffusion-weighted volumes and checks that the
   * diffusion-weighted volumes form one shell.
   *
   * A volume is a b0 volume when its b-value rounds to 0 (DiffusionPropertyHelper::IsB0), as in the
   * b-value map the views show. The single-shell models normalize the gradients and use one
   * b-value, so a second shell would be fitted with the wrong b-value.
   *
   * \throw itk::ExceptionObject without a b0 volume, without a diffusion-weighted volume, or when the
   * diffusion-weighted volumes have more than one b-value.
   */
  inline SingleShellGradients GetSingleShellGradients(const BValueMapProperty::GradientDirectionsContainerType *gradients,
                                                      double referenceBValue)
  {
    SingleShellGradients result;
    std::set<double> shells;
    for (unsigned int i = 0; i < gradients->Size(); ++i)
    {
      const auto &g = gradients->ElementAt(i);
      if (DiffusionPropertyHelper::IsB0(g, referenceBValue))
      {
        result.baselineIndices.push_back(i);
      }
      else
      {
        result.diffusionWeightedIndices.push_back(i);
        result.directions.push_back(g / g.two_norm());
        shells.insert(BValueMapProperty::GetBValue(g, referenceBValue));
      }
    }

    if (result.baselineIndices.empty())
      itkGenericExceptionMacro("No baseline image detected (no b-zero image): no volume has a b-value that rounds to 0.");
    if (result.diffusionWeightedIndices.empty())
      itkGenericExceptionMacro("No diffusion-weighted volume: every b-value rounds to 0.");
    if (shells.size() > 1)
    {
      std::ostringstream bValues;
      for (auto b : shells)
        bValues << (b == *shells.begin() ? "" : ", ") << b;
      itkGenericExceptionMacro("This reconstruction models one shell, but the diffusion-weighted volumes have "
                               << shells.size() << " b-values (" << bValues.str()
                               << "). Extract one shell, or round the b-values if they differ only by noise.");
    }
    result.bValue = *shells.begin();
    return result;
  }
}

#endif
