/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_OdfNormalizationMethodProperty__H_
#define _MITK_OdfNormalizationMethodProperty__H_

#include "mitkEnumerationProperty.h"
#include "mitkEnumerationPropertySerializer.h"
#include <MitkDiffusionModellingExports.h>

#include "mitkSerializerMacros.h"
#include <tinyxml2.h>

namespace mitk
{

  enum OdfNormalizationMethod
  {
    ODFN_MINMAX,
    ODFN_MAX,
    ODFN_NONE,
    ODFN_GLOBAL_MAX
  };

/**
 * Encapsulates the enumeration for ODF normalization. Valid values are
 * ODFN_MINMAX, ODFN_MAX, ODFN_NONE ODFN_GLOBAL_MAX
 * Default is ODFN_MINMAX
 */
class MITKDIFFUSIONMODELLING_EXPORT OdfNormalizationMethodProperty : public EnumerationProperty
{
public:

  mitkClassMacro( OdfNormalizationMethodProperty, EnumerationProperty );
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  mitkNewMacro1Param(OdfNormalizationMethodProperty, const IdType&);

  mitkNewMacro1Param(OdfNormalizationMethodProperty, const std::string&);

  /**
   * Returns the current interpolation value as defined by VTK constants.
   */
  virtual int GetNormalization();

  /**
   * Sets the interpolation type to ODFN_MINMAX.
   */
  virtual void SetNormalizationToMinMax();

  /**
   * Sets the interpolation type to ODFN_MAX.
   */
  virtual void SetNormalizationToMax();

  /**
   * Sets the interpolation type to ODFN_NONE.
   */
  virtual void SetNormalizationToNone();

  /**
   * Sets the interpolation type to ODFN_GLOBAL_MAX.
   */
  virtual void SetNormalizationToGlobalMax();

protected:

  /** Sets reslice interpolation mode to default (VTK_RESLICE_NEAREST).
   */
  OdfNormalizationMethodProperty( );

  /**
   * Constructor. Sets reslice interpolation to the given value.
   */
  OdfNormalizationMethodProperty( const IdType& value );

  /**
   * Constructor. Sets reslice interpolation to the given value.
   */
  OdfNormalizationMethodProperty( const std::string& value );

  /**
   * this function is overridden as protected, so that the user may not add
   * additional invalid interpolation types.
   */
  bool AddEnum( const std::string& name, const IdType& id ) override;

  /**
   * Adds the enumeration types as defined by vtk to the list of known
   * enumeration values.
   */
  virtual void AddInterpolationTypes();
};

class MITKDIFFUSIONMODELLING_EXPORT OdfNormalizationMethodPropertySerializer : public EnumerationPropertySerializer
{
  public:

    mitkClassMacro( OdfNormalizationMethodPropertySerializer, EnumerationPropertySerializer );
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    BaseProperty::Pointer Deserialize(const tinyxml2::XMLElement* element) override
    {
      if (!element) return nullptr;

      std::string s = "";
      s = std::string(element->Attribute("value"));

      OdfNormalizationMethodProperty::Pointer property = OdfNormalizationMethodProperty::New();
      property->SetValue( s );
      return property.GetPointer();
    }

  protected:

    OdfNormalizationMethodPropertySerializer () {}
    ~OdfNormalizationMethodPropertySerializer () override {}
};


} // end of namespace mitk

#endif


