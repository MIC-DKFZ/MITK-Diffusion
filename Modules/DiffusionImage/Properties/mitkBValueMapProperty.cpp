/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkBValueMapProperty.h"

#include <nlohmann/json.hpp>
#include <string>

#include <cmath>

mitk::BValueMapProperty::BValueMapProperty() = default;

mitk::BValueMapProperty::BValueMapProperty(const BValueMapProperty& other)
  : mitk::BaseProperty(other)
{
  m_BValueMap = other.GetBValueMap();
}

mitk::BValueMapProperty::BValueMapProperty(const BValueMap& bValueMap)
{
  m_BValueMap = bValueMap;
}

mitk::BValueMapProperty::~BValueMapProperty() = default;

const mitk::BValueMapProperty::BValueMap & mitk::BValueMapProperty::GetBValueMap() const
{
  return m_BValueMap;
}

void mitk::BValueMapProperty::SetBValueMap(const BValueMap & map)
{
  this->m_BValueMap = map;
}

bool mitk::BValueMapProperty::IsEqual(const BaseProperty& property) const
{
  return this->m_BValueMap == static_cast<const Self&>(property).m_BValueMap;
}

bool mitk::BValueMapProperty::Assign(const BaseProperty& property)
{
  this->m_BValueMap = static_cast<const Self&>(property).m_BValueMap;
  return true;
}

// {"<b-value>": [volume indices], ...}: the form clients of the REST API see.
bool mitk::BValueMapProperty::ToJSON(nlohmann::json& j) const
{
  j = nlohmann::json::object();

  for (const auto& [bValue, indices] : m_BValueMap)
    j[std::to_string(bValue)] = indices;

  return true;
}

bool mitk::BValueMapProperty::FromJSON(const nlohmann::json& j)
{
  BValueMap map;

  for (const auto& [key, indices] : j.items())
    map[static_cast<unsigned int>(std::stoul(key))] = indices.get<std::vector<unsigned int>>();

  m_BValueMap = map;
  return true;
}

mitk::BValueMapProperty::BValueMap mitk::BValueMapProperty::CreateBValueMap(const mitk::BValueMapProperty::GradientDirectionsContainerType * gdc, float referenceBValue)
{
  mitk::BValueMapProperty::BValueMap map;

  mitk::BValueMapProperty::GradientDirectionsContainerType::ConstIterator gdcit;
  for( gdcit = gdc->Begin(); gdcit != gdc->End(); ++gdcit)
    map[GetBValue(gdcit.Value(), referenceBValue)].push_back(gdcit.Index());

  return map;
}

double mitk::BValueMapProperty::GetBValue(const GradientDirectionType& gradient, double referenceBValue)
{
  if (gradient.one_norm() <= 0.0)
    return 0;

  const double twonorm = gradient.two_norm();
  const double bValue = referenceBValue*twonorm*twonorm;
  return bValue < 0 ? std::ceil(bValue - 0.5) : std::floor(bValue + 0.5);
}

itk::LightObject::Pointer  mitk::BValueMapProperty::InternalClone() const
{
  itk::LightObject::Pointer result(new Self(*this));
  result->UnRegister();
  return result;
}
