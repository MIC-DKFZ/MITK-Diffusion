/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkImageCast.h>
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include "mitkCommandLineParser.h"
#include <mitkDiffusionImageMapping.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkModelDirectionRotation.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkShImage.h>
#include <vnl/vnl_inverse.h>

/*!
\brief Copies transformation matrix of one image to another

The directions of a DWI (gradients) and of tensor, SH, ODF and peak images are world vectors, so
they are turned by the rotation between the old and the new geometry.
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("Copy Geometry");
  parser.setCategory("Preprocessing Tools");
  parser.setDescription("Copies transformation matrix of one image to another. The directions of a DWI (gradients) and of tensor, SH, ODF and peak images turn with the image.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input:", "input image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output image", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("ref", "r", mitkCommandLineParser::String, "Reference:", "reference image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("alignCentroid", "a", mitkCommandLineParser::Bool, "align centroids", "align centroids", us::Any(), true);
  parser.addArgument("input_type", "t", mitkCommandLineParser::String, "Input type:", "read the input as 'sh' or 'peak' image; 4D NIfTI files are otherwise read as plain images, whose values do not turn", us::Any(), true);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  // mandatory arguments
  std::string imageName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string refImage = us::any_cast<std::string>(parsedArgs["ref"]);
  std::string outImage = us::any_cast<std::string>(parsedArgs["o"]);

  bool originOnly = false;

  if ( parsedArgs.count("alignCentroid") || parsedArgs.count("a"))
  {
    originOnly = true;
  }

  // A 4D SH or peak NIfTI file cannot be told from other 4D images by its content, so MITK's image
  // reader is its default reader; only the named reader turns its directions.
  std::string inputType;
  if (parsedArgs.count("input_type"))
    inputType = us::any_cast<std::string>(parsedArgs["input_type"]);
  std::vector<std::string> readers;
  if (inputType == "sh")
    readers = {"SH Image"};
  else if (inputType == "peak")
    readers = {"Peak Image"};
  else if (!inputType.empty())
  {
    std::cout << "Unknown input type " << inputType << "; use 'sh' or 'peak'." << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    mitk::Image::Pointer source = mitk::IOUtil::Load<mitk::Image>(refImage);
    mitk::PreferenceListReaderOptionsFunctor functor(readers, std::vector<std::string>());
    mitk::Image::Pointer target = mitk::IOUtil::Load<mitk::Image>(imageName, &functor);
    // The reader list only prefers the named reader; another one (e.g. the tensor reader of a .dti
    // file) can still win, so the loaded type is checked against the requested one.
    const bool typeMatches = inputType.empty() ||
      (inputType == "sh" && dynamic_cast<mitk::ShImage*>(target.GetPointer()) != nullptr) ||
      (inputType == "peak" && dynamic_cast<mitk::PeakImage*>(target.GetPointer()) != nullptr);
    if (!typeMatches)
      mitkThrow() << "The input cannot be read as " << inputType << " image.";
    if (originOnly)
    {
      // Shift the target so that the world centres of both voxel grids coincide.
      const mitk::Vector3D shift = source->GetGeometry()->GetCenter() - target->GetGeometry()->GetCenter();
      const mitk::Point3D newOrigin = target->GetGeometry()->GetOrigin() + shift;
      target->GetGeometry()->SetOrigin(newOrigin);
    }
    else
    {
      mitk::BaseGeometry* s_geom = source->GetGeometry();
      mitk::BaseGeometry* t_geom = target->GetGeometry();
      const mitk::DiffusionImageMapping::MatrixType oldMatrix = t_geom->GetIndexToWorldTransform()->GetMatrix().GetVnlMatrix();

      t_geom->SetIndexToWorldTransform(s_geom->GetIndexToWorldTransform());
      target->SetGeometry(t_geom);

      // The world gradients of a DWI and the world directions of the model images follow the
      // anatomy, so they turn with the image, as in DReg.
      const mitk::DiffusionImageMapping::MatrixType newMatrix = s_geom->GetIndexToWorldTransform()->GetMatrix().GetVnlMatrix();
      const mitk::DiffusionImageMapping::MatrixType change = newMatrix * vnl_inverse(oldMatrix);
      if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(target))
        mitk::DiffusionImageMapping::RotateGradients(target, change);
      else if (mitk::ModelDirectionRotation::HasModelDirections(target))
        mitk::ModelDirectionRotation::Rotate(target, mitk::DiffusionImageMapping::RotationComponent(change));
    }
    mitk::IOUtil::Save(target, outImage);
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cout << "ERROR!?!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
