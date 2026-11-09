/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <itkFlipPeaksFilter.h>
#include <mitkCommandLineParser.h>
#include <mitkDiffusionPeaks.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <map>

namespace
{
  /** Values of --shConvention and the "SH convention" option of the SH reader they select. */
  const std::map<std::string, std::string> SH_CONVENTIONS = {
    {"auto", "auto (file marker, else MITK legacy; vector layout: MITK FSL)"},
    {"MITK", "MITK legacy (MRtrix 0.2)"},
    {"MRtrix", "MRtrix 3"},
    {"MRtrix3", "MRtrix 3"},
    {"FSL", "FSL/Dipy"}};

  /** Negates the selected world components of every peak (FlipPeaksFilter, as the controls view). */
  mitk::PeakImage::Pointer FlipPeaks(mitk::PeakImage *peaks, bool x, bool y, bool z)
  {
    auto flipper = itk::FlipPeaksFilter::New();
    flipper->SetInput(peaks->GetItkPeakImage());
    flipper->SetFlipX(x);
    flipper->SetFlipY(y);
    flipper->SetFlipZ(z);
    flipper->Update();

    auto flipped = mitk::PeakImage::New();
    // GrabItkImageMemory would take over only the first volume of the 4D image.
    mitk::GrabItkImageMemoryChannel(flipper->GetOutput(), nullptr, flipped.GetPointer(), false);
    return flipped;
  }
}

/*!
\brief Extracts the maxima of the ODFs of a spherical harmonics image, as the peak extraction view.
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;
  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input image", "SH coefficient image; the SH order is read from it", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output root", "writes ROOT_PEAKS.nrrd and ROOT_NUM_PEAKS.nrrd", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("mask", "m", mitkCommandLineParser::String, "Mask", "mask image; voxels where it is 0 get no peaks", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("normalization", "n", mitkCommandLineParser::Int, "Normalization", "0=no norm (ODF value), 1=max norm (relative to the largest peak), 2=single vec norm (unit vectors)", 1, true);
  parser.addArgument("numpeaks", "p", mitkCommandLineParser::Int, "Max. number of peaks", "maximum number of extracted peaks per voxel", 2, true);
  parser.addArgument("rel_peakthr", "r", mitkCommandLineParser::Float, "Relative peak threshold", "peaks below this fraction of the largest peak of their voxel are dropped", 0.4, true);
  parser.addArgument("abs_peakthr", "a", mitkCommandLineParser::Float, "Absolute peak threshold", "peaks whose ODF value is below are dropped", 0.03, true);
  parser.addArgument("angular_thr", "", mitkCommandLineParser::Float, "Angular threshold", "in degrees: of two peaks closer than this, the smaller one is dropped", 15, true);
  parser.addArgument("shConvention", "s", mitkCommandLineParser::String, "SH convention", "how the coefficients in the file are defined: auto (file marker, else MITK; MITK's FSL definition for coefficient vectors of a 3D image, as earlier MITK versions saved them), MITK (legacy, MRtrix 0.2), MRtrix3, FSL (dipy's default basis in the FSL bvec frame)", std::string("auto"), true);
  parser.addArgument("flipX", "", mitkCommandLineParser::Bool, "Flip X", "negate the x component (world) of the peaks", us::Any(), true);
  parser.addArgument("flipY", "", mitkCommandLineParser::Bool, "Flip Y", "negate the y component (world) of the peaks", us::Any(), true);
  parser.addArgument("flipZ", "", mitkCommandLineParser::Bool, "Flip Z", "negate the z component (world) of the peaks", us::Any(), true);
  parser.addArgument("scale_by_gfa", "", mitkCommandLineParser::Bool, "Scale by GFA", "scale ODF values and peaks by the GFA", us::Any(), true);
  // Names under which the thresholds were documented before; kept so that existing calls work.
  parser.addArgument("peakthres", "", mitkCommandLineParser::Float, "Relative peak threshold", "older name of --rel_peakthr", us::Any(), true);
  parser.addArgument("abspeakthres", "", mitkCommandLineParser::Float, "Absolute peak threshold", "older name of --abs_peakthr", us::Any(), true);

  parser.setCategory("Preprocessing Tools");
  parser.setTitle("Peak Extraction");
  parser.setDescription("Extract the maxima of the ODFs of a spherical harmonics image (candidates on the ODF directions, each refined to a maximum of the SH function by Newton steps), as the peak extraction view.");
  parser.setContributor("MIC");

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  const std::string imageName = us::any_cast<std::string>(parsedArgs["i"]);
  const std::string outRoot = us::any_cast<std::string>(parsedArgs["o"]);

  std::string maskImageName;
  if (parsedArgs.count("mask"))
    maskImageName = us::any_cast<std::string>(parsedArgs["mask"]);

  mitk::PeakParameters parameters;
  parameters.maxPeaks = 2;
  if (parsedArgs.count("numpeaks"))
  {
    const int numPeaks = us::any_cast<int>(parsedArgs["numpeaks"]);
    if (numPeaks < 1)
    {
      std::cout << "--numpeaks must be at least 1." << std::endl;
      return EXIT_FAILURE;
    }
    parameters.maxPeaks = static_cast<unsigned int>(numPeaks);
  }
  if (parsedArgs.count("peakthres"))
    parameters.relativeThreshold = us::any_cast<float>(parsedArgs["peakthres"]);
  if (parsedArgs.count("rel_peakthr"))
    parameters.relativeThreshold = us::any_cast<float>(parsedArgs["rel_peakthr"]);
  if (parsedArgs.count("abspeakthres"))
    parameters.absoluteThreshold = us::any_cast<float>(parsedArgs["abspeakthres"]);
  if (parsedArgs.count("abs_peakthr"))
    parameters.absoluteThreshold = us::any_cast<float>(parsedArgs["abs_peakthr"]);
  if (parsedArgs.count("angular_thr"))
    parameters.angularThreshold = us::any_cast<float>(parsedArgs["angular_thr"]);
  if (parsedArgs.count("scale_by_gfa"))
    parameters.scaleByGfa = us::any_cast<bool>(parsedArgs["scale_by_gfa"]);
  if (parsedArgs.count("normalization"))
  {
    switch (us::any_cast<int>(parsedArgs["normalization"]))
    {
      case 0:
        parameters.normalization = mitk::PeakNormalization::None;
        break;
      case 1:
        parameters.normalization = mitk::PeakNormalization::MaxVector;
        break;
      case 2:
        parameters.normalization = mitk::PeakNormalization::SingleVector;
        break;
      default:
        std::cout << "--normalization must be 0, 1 or 2." << std::endl;
        return EXIT_FAILURE;
    }
  }

  std::string convention = "auto";
  if (parsedArgs.count("shConvention"))
    convention = us::any_cast<std::string>(parsedArgs["shConvention"]);
  if (SH_CONVENTIONS.count(convention) == 0)
  {
    std::cout << "Unknown --shConvention " << convention << "; use auto, MITK, MRtrix3 or FSL." << std::endl;
    return EXIT_FAILURE;
  }

  const bool flipX = parsedArgs.count("flipX") && us::any_cast<bool>(parsedArgs["flipX"]);
  const bool flipY = parsedArgs.count("flipY") && us::any_cast<bool>(parsedArgs["flipY"]);
  const bool flipZ = parsedArgs.count("flipZ") && us::any_cast<bool>(parsedArgs["flipZ"]);

  try
  {
    // Without the preference a NIfTI SH file is read as a plain 4D image.
    mitk::IFileReader::Options options;
    options["SH convention"] = SH_CONVENTIONS.at(convention);
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, options);
    mitk::ShImage::Pointer shImage = dynamic_cast<mitk::ShImage*>(mitk::IOUtil::Load(imageName, &functor).at(0).GetPointer());
    if (shImage.IsNull())
    {
      std::cout << imageName << " is no spherical harmonics image." << std::endl;
      return EXIT_FAILURE;
    }

    mitk::Image::Pointer mask;
    if (!maskImageName.empty())
      mask = mitk::IOUtil::Load<mitk::Image>(maskImageName);

    std::cout << "SH order " << shImage->ShOrder() << ", up to " << parameters.maxPeaks << " peaks per voxel" << std::endl;
    auto result = mitk::ExtractPeaks(shImage, parameters, mask);

    mitk::PeakImage::Pointer peaks = result.peaks;
    if (flipX || flipY || flipZ)
      peaks = FlipPeaks(peaks, flipX, flipY, flipZ);

    mitk::IOUtil::Save(peaks, outRoot + "_PEAKS.nrrd");
    mitk::IOUtil::Save(result.numberOfPeaks, outRoot + "_NUM_PEAKS.nrrd");
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cout << "ERROR!?!" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
