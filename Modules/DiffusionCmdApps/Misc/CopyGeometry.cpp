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
#include <vnl/vnl_inverse.h>

/*!
\brief Copies transformation matrix of one image to another
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("Copy Geometry");
  parser.setCategory("Preprocessing Tools");
  parser.setDescription("Copies transformation matrix of one image to another");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input:", "input image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output image", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("ref", "r", mitkCommandLineParser::String, "Reference:", "reference image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("alignCentroid", "a", mitkCommandLineParser::Bool, "align centroids", "align centroids", us::Any(), true);

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

  try
  {
    mitk::Image::Pointer source = mitk::IOUtil::Load<mitk::Image>(refImage);
    mitk::Image::Pointer target = mitk::IOUtil::Load<mitk::Image>(imageName);
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

      // The world gradients of a DWI follow the anatomy, so they turn with the image, as in DReg.
      if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(target))
      {
        const mitk::DiffusionImageMapping::MatrixType newMatrix = s_geom->GetIndexToWorldTransform()->GetMatrix().GetVnlMatrix();
        mitk::DiffusionImageMapping::RotateGradients(target, newMatrix * vnl_inverse(oldMatrix));
      }
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
