/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkVectorImage.h>

#include <mitkBaseData.h>
#include <mitkFiberBundle.h>
#include "mitkCommandLineParser.h"
#include <mitkLexicalCast.h>
#include <mitkCoreObjectFactory.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundleIOHelper.h>
#include <mitkImage.h>
#include <mitkRenderingTestHelper.h>
#include <vtkRenderLargeImage.h>
#include <vtkImageWriter.h>
#include <vtkPNGWriter.h>


/*!
\brief Modify input tractogram: fiber resampling, compression, pruning and transformation.
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("Fiber Screenshot");
  parser.setDescription("Take a screenshot of the loaded fiber bundle");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");

  parser.beginGroup("1. Mandatory arguments:");
  parser.addArgument("", "i", mitkCommandLineParser::StringList, "Input:", "input tractograms", us::Any(), false);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "Output png", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("resample", "", mitkCommandLineParser::Float, "", "");
  parser.endGroup();

  parser.addArgument("rotate_x", "", mitkCommandLineParser::Float, "Rotate x-axis:", "Rotate around x-axis (in deg)");
  parser.addArgument("rotate_y", "", mitkCommandLineParser::Float, "Rotate y-axis:", "Rotate around y-axis (in deg)");
  parser.addArgument("rotate_z", "", mitkCommandLineParser::Float, "Rotate z-axis:", "Rotate around z-axis (in deg)");
  parser.addArgument("background", "", mitkCommandLineParser::String, "Background:", "Background color (WHITE, BLACK)");


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;


  mitkCommandLineParser::StringContainerType inFibs = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["i"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  float rotateX = 0;
  if (parsedArgs.count("rotate_x"))
    rotateX = us::any_cast<float>(parsedArgs["rotate_x"]);

  float rotateY = 0;
  if (parsedArgs.count("rotate_y"))
    rotateY = us::any_cast<float>(parsedArgs["rotate_y"]);

  float rotateZ = 0;
  if (parsedArgs.count("rotate_z"))
    rotateZ = us::any_cast<float>(parsedArgs["rotate_z"]);

  std::string color = "BLACK";
  if (parsedArgs.count("background"))
    color = us::any_cast<std::string>(parsedArgs["background"]);

  try
  {
    mitk::RenderingTestHelper renderingHelper(640*2, 480*2);
    std::vector< mitk::FiberBundle::Pointer > tractograms = mitk::FiberBundleIOHelper::load_fibs(inFibs);
    for (auto fib : tractograms)
    {
      fib->RotateAroundAxis(rotateX, rotateY, rotateZ);
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(fib);

      node->SetFloatProperty("shape.tuberadius", 0.5);
      renderingHelper.AddNodeToStorage(node);
    }

    renderingHelper.SetMapperIDToRender3D();

    vtkSmartPointer<vtkRenderer> renderer = renderingHelper.GetVtkRenderer();
    renderer->GetRenderWindow()->DoubleBufferOn();
    double bgcolor[] = {0.0, 0.0, 0.0};
    if (color != "BLACK")
    {
      bgcolor[0] = 1.0;
      bgcolor[1] = 1.0;
      bgcolor[2] = 1.0;
    }

    renderer->SetBackground(bgcolor);

    vtkSmartPointer<vtkRenderLargeImage> magnifier = vtkSmartPointer<vtkRenderLargeImage>::New();
    magnifier->SetInput(renderer);
    magnifier->SetMagnification(1);

    vtkSmartPointer<vtkImageWriter> fileWriter = vtkSmartPointer<vtkPNGWriter>::New();
    fileWriter->SetInputConnection(magnifier->GetOutputPort());
    fileWriter->SetFileName(outFileName.c_str());

    fileWriter->Write();

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
