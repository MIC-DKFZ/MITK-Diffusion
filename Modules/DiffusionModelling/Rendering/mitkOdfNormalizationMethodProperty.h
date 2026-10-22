/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkOdfNormalizationMethodProperty_h
#define mitkOdfNormalizationMethodProperty_h

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
 * How the ODF glyphs are normalized before they are scaled: "Min-Max" (ODFN_MINMAX) maps each
 * ODF to [0, 1], "Maximum" (ODFN_MAX) divides it by its maximum, "None" (ODFN_NONE) keeps the
 * values. "Global Maximum" (ODFN_GLOBAL_MAX) is kept for scene files that hold it; it normalizes
 * like "Maximum". Default is ODFN_MAX.
 */
class MITKDIFFUSIONMODELLING_EXPORT OdfNormalizationMethodProperty : public EnumerationProperty
{
public:

  mitkClassMacro( OdfNormalizationMethodProperty, EnumerationProperty );
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  mitkCloneMacro(Self)

  mitkNewMacro1Param(OdfNormalizationMethodProperty, const IdType&);

  mitkNewMacro1Param(OdfNormalizationMethodProperty, const std::string&);

  /** The current OdfNormalizationMethod. */
  virtual int GetNormalization() const;

  virtual void SetNormalizationToMinMax();
  virtual void SetNormalizationToMax();
  virtual void SetNormalizationToNone();
  virtual void SetNormalizationToGlobalMax();

protected:

  /** Normalization ODFN_MAX. */
  OdfNormalizationMethodProperty( );

  OdfNormalizationMethodProperty(const OdfNormalizationMethodProperty&) = default;

  /** The given normalization; ODFN_MINMAX if the value is not one of OdfNormalizationMethod. */
  OdfNormalizationMethodProperty( const IdType& value );

  /** The normalization of the given name; ODFN_MINMAX if the name is unknown. */
  OdfNormalizationMethodProperty( const std::string& value );

  /** Protected, so that no other values can be added. */
  bool AddEnum( const std::string& name, const IdType& id ) override;

  /** Adds the values of OdfNormalizationMethod with their names. */
  virtual void AddNormalizationMethods();
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

      // A malformed scene file may lack the attribute.
      const char* value = element->Attribute("value");
      if (value == nullptr) return nullptr;
      std::string s(value);

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


