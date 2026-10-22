/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkOdfScaleByProperty_h
#define mitkOdfScaleByProperty_h

#include <MitkDiffusionModellingExports.h>
#include "mitkEnumerationProperty.h"
#include "mitkEnumerationPropertySerializer.h"
#include <tinyxml2.h>

namespace mitk
{

  enum OdfScaleBy
  {
    ODFSB_NONE,
    ODFSB_GFA,
    ODFSB_PC
  };

/**
 * By what the ODF glyphs are scaled in addition to the "Scaling" property: "None" (ODFSB_NONE),
 * the generalized fractional anisotropy "GFA" (ODFSB_GFA) or the "Principal Curvature"
 * (ODFSB_PC) of each ODF. Default is ODFSB_NONE.
 */
class MITKDIFFUSIONMODELLING_EXPORT OdfScaleByProperty : public EnumerationProperty
{
public:

  mitkClassMacro( OdfScaleByProperty, EnumerationProperty );
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  mitkCloneMacro(Self)

  mitkNewMacro1Param(OdfScaleByProperty, const IdType&);

  mitkNewMacro1Param(OdfScaleByProperty, const std::string&);

  /** The current OdfScaleBy value. */
  virtual int GetScaleBy() const;

  virtual void SetScaleByNothing();
  virtual void SetScaleByGFA();
  virtual void SetScaleByPrincipalCurvature();

protected:

  /** Scaling ODFSB_NONE. */
  OdfScaleByProperty( );

  OdfScaleByProperty(const OdfScaleByProperty&) = default;

  /** The given scaling; ODFSB_NONE if the value is not one of OdfScaleBy. */
  OdfScaleByProperty( const IdType& value );

  /** The scaling of the given name; ODFSB_NONE if the name is unknown. */
  OdfScaleByProperty( const std::string& value );

  /** Protected, so that no other values can be added. */
  bool AddEnum( const std::string& name, const IdType& id ) override;

  /** Adds the values of OdfScaleBy with their names. */
  virtual void AddScaleByTypes();
};

class MITKDIFFUSIONMODELLING_EXPORT OdfScaleByPropertySerializer : public EnumerationPropertySerializer
{
  public:

    mitkClassMacro( OdfScaleByPropertySerializer, EnumerationPropertySerializer );
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    BaseProperty::Pointer Deserialize(const tinyxml2::XMLElement* element) override
    {
      if (!element) return nullptr;

      // A malformed scene file may lack the attribute.
      const char* value = element->Attribute("value");
      if (value == nullptr) return nullptr;
      std::string s(value);

      OdfScaleByProperty::Pointer property = OdfScaleByProperty::New();
      property->SetValue( s );
      return property.GetPointer();
    }

  protected:

    OdfScaleByPropertySerializer () {}
    ~OdfScaleByPropertySerializer () override {}
};

} // end of namespace mitk

#endif


