/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkMeasurementFrameProperty.h"

#include <nlohmann/json.hpp>

mitk::MeasurementFrameProperty::MeasurementFrameProperty()
  : m_MeasurementFrame()
{
}

mitk::MeasurementFrameProperty::MeasurementFrameProperty(const MeasurementFrameProperty& other)
  : mitk::BaseProperty(other)
{
  m_MeasurementFrame = other.GetMeasurementFrame();
}

mitk::MeasurementFrameProperty::MeasurementFrameProperty(const MeasurementFrameType& measurementFrame)
{
  m_MeasurementFrame = measurementFrame;
}

mitk::MeasurementFrameProperty::~MeasurementFrameProperty()
{
}

const mitk::MeasurementFrameProperty::MeasurementFrameType & mitk::MeasurementFrameProperty::GetMeasurementFrame() const
{
  return m_MeasurementFrame;
}

void mitk::MeasurementFrameProperty::SetMeasurementFrame(const MeasurementFrameType & frame)
{
  m_MeasurementFrame = frame;
}

bool mitk::MeasurementFrameProperty::IsEqual(const BaseProperty& property) const
{
  return this->m_MeasurementFrame == static_cast<const Self&>(property).m_MeasurementFrame;
}

bool mitk::MeasurementFrameProperty::Assign(const BaseProperty& property)
{
  this->m_MeasurementFrame = static_cast<const Self&>(property).m_MeasurementFrame;
  return true;
}

// Three rows of three values.
bool mitk::MeasurementFrameProperty::ToJSON(nlohmann::json& j) const
{
  j = nlohmann::json::array();

  for (unsigned int row = 0; row < 3; ++row)
    j.push_back({ m_MeasurementFrame(row, 0), m_MeasurementFrame(row, 1), m_MeasurementFrame(row, 2) });

  return true;
}

bool mitk::MeasurementFrameProperty::FromJSON(const nlohmann::json& j)
{
  MeasurementFrameType frame;

  for (unsigned int row = 0; row < 3; ++row)
    for (unsigned int column = 0; column < 3; ++column)
      frame(row, column) = j.at(row).at(column).get<double>();

  m_MeasurementFrame = frame;
  return true;
}

itk::LightObject::Pointer  mitk::MeasurementFrameProperty::InternalClone() const
{
  itk::LightObject::Pointer result(new Self(*this));
  result->UnRegister();
  return result;
}
