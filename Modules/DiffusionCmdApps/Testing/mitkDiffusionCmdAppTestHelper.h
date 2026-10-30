/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionCmdAppTestHelper_h
#define mitkDiffusionCmdAppTestHelper_h

#include <mitkDiffusionTestData.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestingMacros.h>

#include <itksys/Process.h>
#include <tinyxml2.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

/**
 * Runs the command line apps of mitk-diffusion as separate processes for the smoke tests of
 * Modules/DiffusionCmdApps. The apps are in the folder of the test driver (MITK_DIFFUSION_CMDAPP_DIR,
 * set by Testing/CMakeLists.txt); each test case runs in a fresh folder below the temporary folder.
 */
namespace mitk::CmdAppTest
{
  struct Result
  {
    /** Exit code, or -1 if the process did not exit normally (crash, exception, timeout). */
    int exitCode = -1;
    /** stdout and stderr together. */
    std::string output;
  };

  inline std::string AppPath(const std::string &app)
  {
    return (std::filesystem::path(MITK_DIFFUSION_CMDAPP_DIR) / ("Mitk" + app)).string();
  }

  /** A fresh, empty working folder for one test case. */
  inline std::string WorkDir(const std::string &name)
  {
    const auto dir = std::filesystem::temp_directory_path() / "mitkDiffusionCmdAppTests" / name;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir.string();
  }

  /** Runs Mitk<app> with the arguments in the working folder; at most 300 s. */
  inline Result Run(const std::string &app, const std::vector<std::string> &args, const std::string &workDir)
  {
    const std::string executable = AppPath(app);
    std::vector<const char *> argv{executable.c_str()};
    for (const auto &arg : args)
      argv.push_back(arg.c_str());
    argv.push_back(nullptr);

    itksysProcess *process = itksysProcess_New();
    itksysProcess_SetCommand(process, argv.data());
    itksysProcess_SetWorkingDirectory(process, workDir.c_str());
    itksysProcess_SetTimeout(process, 300);
    itksysProcess_Execute(process);

    Result result;
    char *data = nullptr;
    int length = 0;
    while (itksysProcess_WaitForData(process, &data, &length, nullptr) != itksysProcess_Pipe_None)
      result.output.append(data, length);
    itksysProcess_WaitForExit(process, nullptr);
    if (itksysProcess_GetState(process) == itksysProcess_State_Exited)
      result.exitCode = itksysProcess_GetExitValue(process);
    itksysProcess_Delete(process);

    MITK_INFO << "Mitk" << app << " exited with " << result.exitCode;
    return result;
  }

  /** Runs the app and fails the test unless it exits with 0; the app output is in the message. */
  inline Result RunOk(const std::string &app, const std::vector<std::string> &args, const std::string &workDir)
  {
    auto result = Run(app, args, workDir);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mitk" + app + " failed:\n" + result.output, 0, result.exitCode);
    return result;
  }

  /** Runs the app with arguments it must refuse: a normal exit with a code other than 0. */
  inline void RunFails(const std::string &app, const std::vector<std::string> &args, const std::string &workDir)
  {
    const auto result = Run(app, args, workDir);
    CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " did not exit normally:\n" + result.output, result.exitCode != -1);
    CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " accepted wrong arguments:\n" + result.output, result.exitCode != 0);
  }

  /**
   * --help and --xml exit with 0; --xml describes the app as a CTK command line module, in
   * well-formed XML: MITK's parser writes the help texts unescaped, so a '<' or '&' in them breaks
   * the description for the Workbench. Parser warnings printed before the XML are skipped.
   */
  inline void CheckHelpAndXml(const std::vector<std::string> &apps)
  {
    const std::string dir = WorkDir("help");
    for (const auto &app : apps)
    {
      const auto help = RunOk(app, {"--help"}, dir);
      CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " --help:\n" + help.output, help.output.find("Use --help") != std::string::npos);
      const auto xml = RunOk(app, {"--xml"}, dir);
      const auto begin = xml.output.find("<executable>");
      const auto end = xml.output.find("</executable>");
      CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " --xml:\n" + xml.output, begin != std::string::npos && end != std::string::npos);
      tinyxml2::XMLDocument document;
      const std::string description = xml.output.substr(begin, end + std::string("</executable>").size() - begin);
      CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " --xml is not well-formed:\n" + xml.output, document.Parse(description.c_str()) == tinyxml2::XML_SUCCESS);
      // tinyxml2 accepts a bare '&' as text, strict XML readers do not; the parser writes no
      // entities, so the description must hold no '&' at all.
      CPPUNIT_ASSERT_MESSAGE("Mitk" + app + " --xml contains '&':\n" + xml.output, description.find('&') == std::string::npos);
    }
  }

  inline std::string Data(const std::string &path)
  {
    return GetDiffusionTestDataFilePath(path);
  }

  inline std::string Path(const std::string &workDir, const std::string &file)
  {
    return (std::filesystem::path(workDir) / file).string();
  }

  /** Loads an output as T, with the named readers preferred (e.g. "Peak Image" for a 4D file). */
  template <class T = mitk::BaseData>
  typename T::Pointer Load(const std::string &workDir, const std::string &file, const std::vector<std::string> &readers = {})
  {
    const std::string path = Path(workDir, file);
    CPPUNIT_ASSERT_MESSAGE("Missing output " + path, std::filesystem::exists(path));
    mitk::PreferenceListReaderOptionsFunctor functor(readers, std::vector<std::string>());
    typename T::Pointer data = dynamic_cast<T *>(mitk::IOUtil::Load(path, &functor).at(0).GetPointer());
    CPPUNIT_ASSERT_MESSAGE("Output " + path + " is no " + T::GetStaticNameOfClass(), data.IsNotNull());
    return data;
  }

  /** The voxel values of a 3D scalar image of any pixel type. */
  inline itk::Image<double, 3>::Pointer Voxels(const mitk::Image *image)
  {
    itk::Image<double, 3>::Pointer voxels;
    mitk::CastToItkImage(image, voxels);
    return voxels;
  }

  inline unsigned int NumFibers(const std::string &workDir, const std::string &file)
  {
    return Load<mitk::FiberBundle>(workDir, file)->GetNumFibers();
  }

  inline void CheckSize(const mitk::Image *image, unsigned int x, unsigned int y, unsigned int z)
  {
    CPPUNIT_ASSERT_EQUAL(x, image->GetDimension(0));
    CPPUNIT_ASSERT_EQUAL(y, image->GetDimension(1));
    CPPUNIT_ASSERT_EQUAL(z, image->GetDimension(2));
  }

  inline std::string ReadText(const std::string &workDir, const std::string &file)
  {
    const std::string path = Path(workDir, file);
    CPPUNIT_ASSERT_MESSAGE("Missing output " + path, std::filesystem::exists(path));
    std::ifstream stream(path);
    return std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
  }
}

#endif
