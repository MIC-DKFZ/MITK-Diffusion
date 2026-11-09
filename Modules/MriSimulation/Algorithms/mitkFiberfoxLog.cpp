/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberfoxLog.h"

#include <itkEventObject.h>
#include <algorithm>
#include <cmath>
#include <iostream>

mitk::FiberfoxLog::FiberfoxLog()
  : m_StartTime(std::chrono::steady_clock::now())
{
}

mitk::FiberfoxLog::~FiberfoxLog() = default;

bool mitk::FiberfoxLog::OpenFile(const std::string& fileName)
{
  std::lock_guard<std::mutex> lock(m_TextMutex);
  if (m_File.is_open())
    m_File.close();
  m_File.open(fileName);
  return m_File.is_open();
}

void mitk::FiberfoxLog::CloseFile()
{
  std::lock_guard<std::mutex> lock(m_TextMutex);
  if (m_File.is_open())
    m_File.close();
}

void mitk::FiberfoxLog::StartTimer()
{
  std::lock_guard<std::mutex> lock(m_TextMutex);
  m_StartTime = std::chrono::steady_clock::now();
}

std::string mitk::FiberfoxLog::GetTime() const
{
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - m_StartTime;
  const auto total = static_cast<unsigned long>(std::round(std::max(elapsed.count(), 0.0)));
  return std::to_string(total/3600) + ":" + std::to_string((total%3600)/60) + ":" + std::to_string(total%60);
}

void mitk::FiberfoxLog::Print(const std::string& m, bool addTime, bool lineBreak, bool stdOut)
{
  std::lock_guard<std::mutex> lock(m_TextMutex);

  std::string text;
  if (addTime)
    text = this->GetTime() + " > ";
  text += m;
  if (lineBreak)
    text += "\n";

  m_StatusText += text;
  if (m_File.is_open())
  {
    m_File << text;
    m_File.flush();
  }
  // flushed: the stars of a progress bar have no line break, and a terminal shows stdout line by line
  if (stdOut)
    std::cout << text << std::flush;
}

std::string mitk::FiberfoxLog::GetStatusText() const
{
  std::lock_guard<std::mutex> lock(m_TextMutex);
  return m_StatusText;
}

void mitk::FiberfoxLog::StartProgress(std::uint64_t total)
{
  {
    std::lock_guard<std::mutex> lock(m_ProgressMutex);
    m_ProgressTotal = std::max<std::uint64_t>(total, 1);
    m_ProgressCount = 0;
    m_ProgressStars = 0;
    m_Progress = 0.0;
  }
  this->Print("0%   10   20   30   40   50   60   70   80   90   100%", false, true, true);
  this->Print("|----|----|----|----|----|----|----|----|----|----|\n*", false, false, true);
}

void mitk::FiberfoxLog::AdvanceProgress(std::uint64_t steps)
{
  std::lock_guard<std::mutex> lock(m_ProgressMutex);
  m_ProgressCount += steps;
  const std::uint64_t done = std::min(m_ProgressCount, m_ProgressTotal);
  const std::uint64_t stars = 50*done/m_ProgressTotal;
  m_Progress = static_cast<double>(done)/m_ProgressTotal;
  if (stars == m_ProgressStars)
    return;

  // Printed and reported under the progress lock, so the stars appear in order and the observers
  // run one at a time.
  this->Print(std::string(stars - m_ProgressStars, '*'), false, false, true);
  m_ProgressStars = stars;
  this->InvokeEvent(itk::ProgressEvent());
}

double mitk::FiberfoxLog::GetProgress() const
{
  return m_Progress;
}
