/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_OdfScaleByProperty__H_
#define _MITK_OdfScaleByProperty__H_

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
 * Encapsulates the enumeration for ODF normalization. Valid values are
 *     ODFSB_NONE, ODFSB_GFA, ODFSB_PC
 * Default is ODFSB_NONE
 */
class MITKDIFFUSIONMODELLING_EXPORT OdfScaleByProperty : public EnumerationProperty
{
public:

  mitkClassMacro( OdfScaleByProperty, EnumerationProperty );
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  mitkNewMacro1Param(OdfScaleByProperty, const IdType&);

  mitkNewMacro1Param(OdfScaleByProperty, const std::string&);

  /**
   * Returns the current interpolation value as defined by VTK constants.
   */
  virtual int GetScaleBy();

  virtual void SetScaleByNothing();
  virtual void SetScaleByGFA();
  virtual void SetScaleByPrincipalCurvature();

protected:

  /** Sets reslice interpolation mode to default (ODFSB_NONE).
   */
  OdfScaleByProperty( );

  /**
   * Constructor. Sets reslice interpolation to the given value.
   */
  OdfScaleByProperty( const IdType& value );

  /**
   * Constructor. Sets reslice interpolation to the given value.
   */
  OdfScaleByProperty( const std::string& value );

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

class MITKDIFFUSIONMODELLING_EXPORT OdfScaleByPropertySerializer : public EnumerationPropertySerializer
{
  public:

    mitkClassMacro( OdfScaleByPropertySerializer, EnumerationPropertySerializer );
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    BaseProperty::Pointer Deserialize(const tinyxml2::XMLElement* element) override
    {
      if (!element) return nullptr;

      std::string s = "";
      s = std::string(element->Attribute("value"));

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


