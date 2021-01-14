/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkMeasurementFramePropertySerializer_h_included
#define mitkMeasurementFramePropertySerializer_h_included

#include "mitkBasePropertySerializer.h"
#include "mitkMeasurementFrameProperty.h"
#include <MitkDiffusionCoreExports.h>
#include <tinyxml2.h>

namespace mitk
{

class MITKDIFFUSIONCORE_EXPORT MeasurementFramePropertySerializer : public BasePropertySerializer
{
  public:

    mitkClassMacro( MeasurementFramePropertySerializer, BasePropertySerializer )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    tinyxml2::XMLElement* Serialize(tinyxml2::XMLDocument& doc) override
    {
      if (const MeasurementFrameProperty* prop = dynamic_cast<const MeasurementFrameProperty*>(m_Property.GetPointer()))
      {

        typedef mitk::MeasurementFrameProperty::MeasurementFrameType MeasurementFrameType;
        const MeasurementFrameType & mft = prop->GetMeasurementFrame();

        if(mft.is_zero()) return nullptr;

        auto  element = doc.NewElement("measurementframe");

        auto  child = doc.NewElement("entry");
        std::stringstream ss;
        ss << mft;
        child->SetAttribute("value", ss.str().c_str());
        element->InsertEndChild(child);

        return element;
      }
      else return nullptr;
    }

    BaseProperty::Pointer Deserialize(const tinyxml2::XMLElement* element) override
    {
      if (!element) return nullptr;

      const tinyxml2::XMLElement* entry = element->FirstChildElement( "entry" )->ToElement();

      std::stringstream ss;
      std::string value;

	  value = std::string(entry->Attribute("value"));
      ss << value;

      MeasurementFrameProperty::MeasurementFrameType matrix;
      matrix.read_ascii(ss);

      return MeasurementFrameProperty::New(matrix).GetPointer();
    }

  protected:

    MeasurementFramePropertySerializer() {}
    ~MeasurementFramePropertySerializer() override {}
};

} // namespace

// important to put this into the GLOBAL namespace (because it starts with 'namespace mitk')
MITK_REGISTER_SERIALIZER(MeasurementFramePropertySerializer)

#endif
