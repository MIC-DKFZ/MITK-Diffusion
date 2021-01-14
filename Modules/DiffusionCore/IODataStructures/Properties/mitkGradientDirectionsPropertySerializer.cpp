/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkGradientDirectionPropertySerializer_h_included
#define mitkGradientDirectionPropertySerializer_h_included

#include "mitkBasePropertySerializer.h"
#include "mitkGradientDirectionsProperty.h"
#include "MitkDiffusionCoreExports.h"
#include <tinyxml2.h>

namespace mitk
{

class MITKDIFFUSIONCORE_EXPORT GradientDirectionsPropertySerializer : public BasePropertySerializer
{
  public:

    mitkClassMacro( GradientDirectionsPropertySerializer, BasePropertySerializer )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    tinyxml2::XMLElement* Serialize(tinyxml2::XMLDocument& doc) override
    {
      if (const GradientDirectionsProperty* prop = dynamic_cast<const GradientDirectionsProperty*>(m_Property.GetPointer()))
      {

        typedef mitk::GradientDirectionsProperty::GradientDirectionsContainerType GradientDirectionsContainerType;
        GradientDirectionsContainerType::ConstPointer gdc = prop->GetGradientDirectionsContainer().GetPointer();

        if(gdc.IsNull() || gdc->Size() == 0) return nullptr;


        GradientDirectionsContainerType::ConstIterator it = gdc->Begin();
        GradientDirectionsContainerType::ConstIterator end = gdc->End();

        auto  element = doc.NewElement("gradientdirections");

        while (it != end) {
          auto  child = doc.NewElement("entry");
          std::stringstream ss;
          ss << it.Value();
          child->SetAttribute("value", ss.str().c_str());
          element->InsertEndChild(child);

          ++it;
        }

        return element;
      }
      else return nullptr;
    }

    BaseProperty::Pointer Deserialize(const tinyxml2::XMLElement* element) override
    {
      if (!element) return nullptr;

      mitk::GradientDirectionsProperty::GradientDirectionsContainerType::Pointer gdc;
      gdc = mitk::GradientDirectionsProperty::GradientDirectionsContainerType::New();

      const tinyxml2::XMLElement* entry = element->FirstChildElement( "entry" )->ToElement();
      while(entry != nullptr){

        std::stringstream ss;
        std::string value;

	    value = std::string(entry->Attribute("value"));
        ss << value;

        vnl_vector_fixed<double, 3> vector;
        vector.read_ascii(ss);

        gdc->push_back(vector);

        entry = entry->NextSiblingElement( "entry" );
      }

      return GradientDirectionsProperty::New(gdc).GetPointer();
    }

  protected:

    GradientDirectionsPropertySerializer() {}
    ~GradientDirectionsPropertySerializer() override {}
};

} // namespace

// important to put this into the GLOBAL namespace (because it starts with 'namespace mitk')
MITK_REGISTER_SERIALIZER(GradientDirectionsPropertySerializer)

#endif
