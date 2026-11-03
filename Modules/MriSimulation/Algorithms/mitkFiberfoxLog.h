/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberfoxLog_h
#define mitkFiberfoxLog_h

#include <MitkMriSimulationExports.h>
#include <itkObject.h>
#include <itkObjectFactory.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <string>

namespace mitk
{
  /**
  * \brief Log of a Fiberfox simulation.
  *
  * Every message goes to the status text (shown by the Fiberfox view while the simulation runs),
  * to the log file next to the result if one is open, and optionally to stdout, with the time
  * since StartTimer(). Each simulation step shows a progress bar of 50 stars; whenever it grows,
  * the log invokes itk::ProgressEvent with GetProgress() as the fraction of the step done.
  *
  * The parallel loops of the simulation write to the log and the view reads the status text from
  * the GUI thread, so all methods are thread-safe. Progress observers run in the simulating
  * thread that advanced the bar, one at a time; they may call GetProgress() and GetStatusText(),
  * but must not write to the log.
  */
  class MITKMRISIMULATION_EXPORT FiberfoxLog : public itk::Object
  {
  public:
    using Self = FiberfoxLog;
    using Superclass = itk::Object;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    itkFactorylessNewMacro(Self);
    itkOverrideGetNameOfClassMacro(FiberfoxLog);

    /** Opens the log file, replacing its content; false if it cannot be opened. */
    bool OpenFile(const std::string& fileName);
    void CloseFile();

    /** The time printed before messages counts from here. */
    void StartTimer();

    /** Appends m, after the elapsed time if addTime and followed by a line break if lineBreak. */
    void Print(const std::string& m, bool addTime = true, bool lineBreak = true, bool stdOut = true);

    /** Everything printed so far. */
    std::string GetStatusText() const;

    /** Starts the progress bar of a simulation step with total steps. */
    void StartProgress(std::uint64_t total);

    /** Steps of the progress bar, one by default. */
    void AdvanceProgress(std::uint64_t steps = 1);

    /** Fraction of the current progress bar done, 0 to 1. */
    double GetProgress() const;

  protected:
    FiberfoxLog();
    ~FiberfoxLog() override;

  private:
    /** Elapsed time as h:m:s; called with m_TextMutex held. */
    std::string GetTime() const;

    mutable std::mutex m_TextMutex;  ///< guards m_StatusText, m_File, m_StartTime and stdout
    std::string m_StatusText;
    std::ofstream m_File;
    std::chrono::steady_clock::time_point m_StartTime;

    std::mutex m_ProgressMutex;  ///< guards the counters and serializes the progress events
    // 64 bits: a simulation counts up to fibers x slices x volumes steps
    std::uint64_t m_ProgressTotal = 1;
    std::uint64_t m_ProgressCount = 0;
    std::uint64_t m_ProgressStars = 0;  ///< stars printed so far
    std::atomic<double> m_Progress{0.0};
  };
}

#endif
