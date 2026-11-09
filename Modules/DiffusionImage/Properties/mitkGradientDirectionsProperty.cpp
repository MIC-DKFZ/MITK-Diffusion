/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkGradientDirectionsProperty.h"

#include <nlohmann/json.hpp>

mitk::GradientDirectionsProperty::GradientDirectionsProperty()
{
  m_GradientDirectionsContainer = mitk::GradientDirectionsProperty::GradientDirectionsContainerType::New();
}

mitk::GradientDirectionsProperty::GradientDirectionsProperty(const GradientDirectionsProperty& other)
  : mitk::BaseProperty(other)
{
  m_GradientDirectionsContainer = other.GetGradientDirectionsContainer();
}

mitk::GradientDirectionsProperty::GradientDirectionsProperty(const GradientDirectionsContainerType::ConstPointer gradientDirectionsContainer)
{
  m_GradientDirectionsContainer = gradientDirectionsContainer;
}

mitk::GradientDirectionsProperty::GradientDirectionsProperty(const GradientDirectionsContainerType::Pointer gradientDirectionsContainer)
{
  m_GradientDirectionsContainer = gradientDirectionsContainer;
}

mitk::GradientDirectionsProperty::~GradientDirectionsProperty() = default;

const mitk::GradientDirectionsProperty::GradientDirectionsContainerType::Pointer mitk::GradientDirectionsProperty::GetGradientDirectionsContainerCopy() const
{
  auto temp = mitk::GradientDirectionsProperty::GradientDirectionsContainerType::New();

  for(unsigned int index = 0; index < m_GradientDirectionsContainer->Size(); ++index)
  {
    GradientDirectionType newDirection;
    newDirection[0] = m_GradientDirectionsContainer->at(index)[0];
    newDirection[1] = m_GradientDirectionsContainer->at(index)[1];
    newDirection[2] = m_GradientDirectionsContainer->at(index)[2];
    temp->InsertElement( index, newDirection);
  }
  return temp;
}

const mitk::GradientDirectionsProperty::GradientDirectionsContainerType::ConstPointer mitk::GradientDirectionsProperty::GetGradientDirectionsContainer() const
{
  return m_GradientDirectionsContainer;
}

bool mitk::GradientDirectionsProperty::IsEqual(const BaseProperty& property) const
{

  const auto& lhs = this->m_GradientDirectionsContainer->CastToSTLConstContainer();
  const auto& rhs = static_cast<const Self&>(property).m_GradientDirectionsContainer->CastToSTLConstContainer();
  return lhs == rhs;
}

bool mitk::GradientDirectionsProperty::Assign(const BaseProperty& property)
{
  this->m_GradientDirectionsContainer = static_cast<const Self&>(property).m_GradientDirectionsContainer;
  return true;
}

// [[x, y, z], ...], one entry per volume.
bool mitk::GradientDirectionsProperty::ToJSON(nlohmann::json& j) const
{
  j = nlohmann::json::array();

  for (unsigned int index = 0; index < m_GradientDirectionsContainer->Size(); ++index)
  {
    const auto& direction = m_GradientDirectionsContainer->ElementAt(index);
    j.push_back({ direction[0], direction[1], direction[2] });
  }

  return true;
}

bool mitk::GradientDirectionsProperty::FromJSON(const nlohmann::json& j)
{
  auto container = GradientDirectionsContainerType::New();

  for (const auto& item : j)
  {
    GradientDirectionType direction;

    for (unsigned int i = 0; i < 3; ++i)
      direction[i] = item.at(i).get<double>();

    container->push_back(direction);
  }

  m_GradientDirectionsContainer = container.GetPointer();
  return true;
}

itk::LightObject::Pointer  mitk::GradientDirectionsProperty::InternalClone() const
{
  itk::LightObject::Pointer result(new Self(*this));
  result->UnRegister();
  return result;
}
