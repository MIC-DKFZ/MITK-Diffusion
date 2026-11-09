/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberWeightsFile.h>

#include <mitkExceptionMacro.h>
#include <mitkFiberBundle.h>

#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>

std::vector<float> mitk::FiberWeightsFile::Read(std::istream& stream)
{
  std::vector<float> weights;
  std::string line;
  unsigned int lineNumber = 0;
  while (std::getline(stream, line))
  {
    ++lineNumber;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#')
      continue;

    std::istringstream tokens(line);
    std::string token;
    while (tokens >> token)
    {
      // std::from_chars is locale independent and, unlike operator>>, rejects trailing characters.
      float value = 0;
      const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
      if (error != std::errc() || end != token.data() + token.size() || !std::isfinite(value))
        mitkThrow() << "Invalid fiber weight '" << token << "' in line " << lineNumber << ".";
      weights.push_back(value);
    }
  }
  return weights;
}

std::vector<float> mitk::FiberWeightsFile::Read(const std::string& filename)
{
  std::ifstream stream(filename);
  if (!stream.is_open())
    mitkThrow() << "Could not open fiber weights file " << filename << ".";
  return Read(stream);
}

void mitk::FiberWeightsFile::Apply(const std::vector<float>& weights, FiberBundle* fib)
{
  if (fib == nullptr)
    mitkThrow() << "No fiber bundle given.";
  if (weights.size() != fib->GetNumFibers())
    mitkThrow() << "The weights file holds " << weights.size() << " weights, but the tractogram has "
                << fib->GetNumFibers() << " fibers.";

  for (std::size_t i = 0; i < weights.size(); ++i)
    fib->SetFiberWeight(static_cast<unsigned int>(i), weights[i]);
}
