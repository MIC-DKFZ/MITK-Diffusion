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
  MITK_TEST(Equal_DifferentGradientDirections_ReturnsFalse);
  MITK_TEST(Deserialize_MissingEntriesOrAttributes_DoesNotCrash);
  MITK_TEST(Json_RoundTripThroughText);


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

    // Values with more than six significant digits catch serialization at default stream precision.
    double a[3] = {0.70710678118654752,-0.57735026918962576,1.4};
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

    // Not symmetric, so a transposed round trip fails.
    double row0[3] = {0.8660254037844386,-0.5,0};
    double row1[3] = {0.5,0.8660254037844386,0};
    double row2[3] = {0,0,1};

    mft.set_row(0,row0);
    mft.set_row(1,row1);
    mft.set_row(2,row2);

    measurementframe_prop = mitk::MeasurementFrameProperty::New(mft).GetPointer();
    propList->ReplaceProperty(mitk::DiffusionPropertyHelper::GetMeasurementFramePropertyName().c_str(), measurementframe_prop);
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

  void Equal_DifferentGradientDirections_ReturnsFalse()
  {
    auto other = mitk::GradientDirectionsProperty::GradientDirectionsContainerType::New();
    for (const auto& g : gradientdirection_prop->GetGradientDirectionsContainer()->CastToSTLConstContainer())
      other->push_back(-g);
    auto otherProp = mitk::GradientDirectionsProperty::New(other);
    CPPUNIT_ASSERT(!(*otherProp == *gradientdirection_prop));
  }

  void Deserialize_MissingEntriesOrAttributes_DoesNotCrash()
  {
    tinyxml2::XMLDocument doc;
    auto noEntries = doc.NewElement("element");
    auto entryWithoutAttributes = doc.NewElement("element");
    entryWithoutAttributes->InsertEndChild(doc.NewElement("entry"));

    for (std::string name : {"GradientDirectionsPropertySerializer", "MeasurementFramePropertySerializer", "BValueMapPropertySerializer"})
    {
      auto serializers = itk::ObjectFactoryBase::CreateAllInstance(name.c_str());
      CPPUNIT_ASSERT_EQUAL(size_t(1), serializers.size());
      auto serializer = dynamic_cast<mitk::BasePropertySerializer *>(serializers.begin()->GetPointer());

      // Without entries, gradients and b-value map are empty; a frame cannot be built.
      auto fromNoEntries = serializer->Deserialize(noEntries);
      CPPUNIT_ASSERT_EQUAL(name != "MeasurementFramePropertySerializer", fromNoEntries.IsNotNull());
      CPPUNIT_ASSERT(serializer->Deserialize(entryWithoutAttributes).IsNull());
    }
  }

  /**
   * The REST API serializes node properties to JSON (mitk::ConvertPropertyToSelfContainedJson):
   * each diffusion property gives a JSON value that reads back to an equal property, also after a
   * round trip through text, and the deserialization service creates the types by name.
   */
  void Json_RoundTripThroughText()
  {
    const auto roundTrip = [](const mitk::BaseProperty* property, mitk::BaseProperty* target) {
      nlohmann::json j;
      CPPUNIT_ASSERT(property->ToJSON(j));
      CPPUNIT_ASSERT(target->FromJSON(nlohmann::json::parse(j.dump())));
      CPPUNIT_ASSERT(*target == *property);
      return j;
    };

    const auto bvalues = roundTrip(bvaluemap_prop, mitk::BValueMapProperty::New());
    CPPUNIT_ASSERT_EQUAL(std::string("[1,2,3,4]"), bvalues.at("0").dump());
    CPPUNIT_ASSERT_EQUAL(std::string("[4,3,2,1]"), bvalues.at("1000").dump());

    const auto gradients = roundTrip(gradientdirection_prop, mitk::GradientDirectionsProperty::New());
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), gradients.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.70710678118654752, gradients.at(0).at(0).get<double>(), 1e-16);

    const auto frame = roundTrip(measurementframe_prop, mitk::MeasurementFrameProperty::New());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, frame.at(0).at(1).get<double>(), 1e-16);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, frame.at(1).at(0).get<double>(), 1e-16);

    // The self-contained form names the type, and reads back through the deserialization service.
    const auto self = mitk::ConvertPropertyToSelfContainedJson(bvaluemap_prop);
    CPPUNIT_ASSERT_EQUAL(std::string("BValueMapProperty"), self.at("type").get<std::string>());
    CPPUNIT_ASSERT(self.at("value") == bvalues);
    const auto back = mitk::ConvertPropertyFromSelfContainedJson(nlohmann::json::parse(self.dump()));
    CPPUNIT_ASSERT(back.IsNotNull() && *back == *bvaluemap_prop);
    const auto frameBack = mitk::ConvertPropertyFromSelfContainedJson(mitk::ConvertPropertyToSelfContainedJson(measurementframe_prop));
    CPPUNIT_ASSERT(frameBack.IsNotNull() && *frameBack == *measurementframe_prop);

    // Malformed values are refused with an exception, not taken in parts.
    auto broken = mitk::MeasurementFrameProperty::New();
    CPPUNIT_ASSERT_THROW(broken->FromJSON(nlohmann::json::parse("[[1,2,3],[4,5]]")), nlohmann::json::exception);
    CPPUNIT_ASSERT_THROW(mitk::GradientDirectionsProperty::New()->FromJSON(nlohmann::json::parse("[[1,2]]")), nlohmann::json::exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionPropertySerializer)
