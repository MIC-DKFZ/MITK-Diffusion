/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkIOUtil.h"
#include "mitkTestingMacros.h"
#include "mitkTestFixture.h"

#include <mitkBasePropertySerializer.h>

#include <mitkBValueMapProperty.h>

#include <mitkGradientDirectionsProperty.h>
#include <mitkMeasurementFrameProperty.h>

#include <mitkDiffusionPropertyHelper.h>
#include <tinyxml2.h>

class mitkDiffusionPropertySerializerTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkDiffusionPropertySerializerTestSuite);
  MITK_TEST(Equal_SerializeandDeserialize_ReturnsTrue);

  //MITK_TEST(Equal_DifferentChannels_ReturnFalse);


  CPPUNIT_TEST_SUITE_END();

private:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/
  mitk::PropertyList::Pointer propList; //represet image propertylist
  mitk::BValueMapProperty::Pointer bvaluemap_prop;
  mitk::GradientDirectionsProperty::Pointer gradientdirection_prop;
  mitk::MeasurementFrameProperty::Pointer measurementframe_prop;

public:

  /**
* @brief Setup Always call this method before each Test-case to ensure correct and new intialization of the used members for a new test case. (If the members are not used in a test, the method does not need to be called).
*/
  void setUp() override
  {

    propList = mitk::PropertyList::New();

    mitk::BValueMapProperty::BValueMap map;
    std::vector<unsigned int> indices1;
    indices1.push_back(1);
    indices1.push_back(2);
    indices1.push_back(3);
    indices1.push_back(4);

    map[0] = indices1;
    std::vector<unsigned int> indices2;
    indices2.push_back(4);
    indices2.push_back(3);
    indices2.push_back(2);
    indices2.push_back(1);

    map[1000] = indices2;
    bvaluemap_prop = mitk::BValueMapProperty::New(map).GetPointer();
    propList->SetProperty(mitk::DiffusionPropertyHelper::GetBvaluePropertyName().c_str(), bvaluemap_prop);

    mitk::GradientDirectionsProperty::GradientDirectionsContainerType::Pointer gdc;
    gdc = mitk::GradientDirectionsProperty::GradientDirectionsContainerType::New();

    double a[3] = {3.0,4.0,1.4};
    vnl_vector_fixed<double,3> vec1;
    vec1.set(a);
    gdc->push_back(vec1);

    double b[3] = {1.0,5.0,123.4};
    vnl_vector_fixed<double,3> vec2;
    vec2.set(b);
    gdc->push_back(vec2);

    double c[3] = {13.0,84.02,13.4};
    vnl_vector_fixed<double,3> vec3;
    vec3.set(c);
    gdc->push_back(vec3);


    gradientdirection_prop = mitk::GradientDirectionsProperty::New(gdc).GetPointer();
    propList->ReplaceProperty(mitk::DiffusionPropertyHelper::GetGradientContainerPropertyName().c_str(), gradientdirection_prop);

    mitk::MeasurementFrameProperty::MeasurementFrameType mft;

    double row0[3] = {1,0,0};
    double row1[3] = {0,1,0};
    double row2[3] = {0,0,1};

    mft.set_row(0,row0);
    mft.set_row(1,row1);
    mft.set_row(2,row2);

    measurementframe_prop = mitk::MeasurementFrameProperty::New(mft).GetPointer();
    propList->ReplaceProperty(mitk::DiffusionPropertyHelper::GetGradientContainerPropertyName().c_str(), measurementframe_prop);
  }

  void tearDown() override
  {

  }

  void Equal_SerializeandDeserialize_ReturnsTrue()
  {

    assert(propList);

    /* try to serialize each property in the list, then deserialize again and check for equality */
    for (mitk::PropertyList::PropertyMap::const_iterator it = propList->GetMap()->begin(); it != propList->GetMap()->end(); ++it)
    {
      const mitk::BaseProperty* property = it->second;

      std::string serializername = std::string(property->GetNameOfClass()) + "Serializer";

      MITK_INFO << "Testing " << serializername;

      std::list<itk::LightObject::Pointer> allSerializers =
          itk::ObjectFactoryBase::CreateAllInstance(serializername.c_str());
      CPPUNIT_ASSERT_EQUAL(size_t(1), allSerializers.size());

      auto *serializer =
          dynamic_cast<mitk::BasePropertySerializer *>(allSerializers.begin()->GetPointer());
      CPPUNIT_ASSERT(serializer != nullptr);

      serializer->SetProperty(property);
      tinyxml2::XMLDocument doc;
      tinyxml2::XMLElement *serialization(nullptr);
      try
      {
        serialization = serializer->Serialize(doc);
      }
      catch (...)
      {
      }
      CPPUNIT_ASSERT(serialization != nullptr);

      mitk::BaseProperty::Pointer restoredProperty = serializer->Deserialize(serialization);
      CPPUNIT_ASSERT(restoredProperty.IsNotNull());

      CPPUNIT_ASSERT(*(restoredProperty.GetPointer()) == *property);

      MITK_INFO << "Testing " << serializername << " done";
    }
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionPropertySerializer)
