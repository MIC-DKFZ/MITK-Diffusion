/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkIOUtil.h>
#include <metaCommand.h>
#include "mitkDiffusionCommandLineParser.h"
#include <usAny.h>
#include <mitkIOUtil.h>
#include <mitkLexicalCast.h>
#include <mitkCoreObjectFactory.h>
#include <mitkPlanarFigure.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkFiberBundle.h>
#include <mitkDiffusionDataIOHelper.h>
#include <mitkTractometry.h>
#include <itksys/SystemTools.hxx>
#include <itkTractParcellationFilter.h>

#define _USE_MATH_DEFINES
#include <math.h>

/*!
\brief Join multiple tractograms
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Fiber Statistics");
  parser.setCategory("Fiber Quantification Methods");
  parser.setContributor("MIC");
  parser.setDescription("Calculate mutliple tractometry methods.");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("tracts", "", mitkDiffusionCommandLineParser::StringList, "Input tracts:", "input tracts", us::Any(), false);
  parser.addArgument("images", "", mitkDiffusionCommandLineParser::StringList, "Input images:", "input images", us::Any(), false);
  parser.addArgument("out_file", "", mitkDiffusionCommandLineParser::String, "Output File:", "output file", us::Any(), false);
  parser.addArgument("num_parcels", "", mitkDiffusionCommandLineParser::Int, "Number of parcels:", "", 15);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  mitkDiffusionCommandLineParser::StringContainerType inFibs = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["tracts"]);
  mitkDiffusionCommandLineParser::StringContainerType inImages = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["images"]);
  std::string outfile = us::any_cast<std::string>(parsedArgs["out_file"]);

  int num_parcels = 15;
  if (parsedArgs.count("num_parcels"))
    num_parcels = us::any_cast<int>(parsedArgs["num_parcels"]);

  try
  {

    std::vector< std::string > fib_names, img_names;
    auto input_tracts = mitk::DiffusionDataIOHelper::load_fibs(inFibs, &fib_names);
    auto input_images = mitk::DiffusionDataIOHelper::load_itk_images<itk::Image<float, 3>>(inImages, &img_names);

    bool add_header = !itksys::SystemTools::FileExists(outfile, true);


    std::ofstream statistics_file;
    statistics_file.open (outfile, std::ios_base::app);
    if (add_header)
      statistics_file << "method;tract_file;image_file;parcel;value" << std::endl;

    auto ref = input_tracts.at(0)->GetDeepCopy();
    for (unsigned int i=0; i<inFibs.size(); ++i)
    {
      auto fib = input_tracts.at(i);
      fib->ResampleSpline(1.0);

      {
        auto output = mitk::Tractometry::StaticResamplingTractometry(input_images.at(i), fib->GetDeepCopy(), num_parcels, ref);
        for (unsigned int r=0; r<output.rows(); ++r)
          for (unsigned int c=0; c<output.cols(); ++c)
            statistics_file << "StaticResamplingTractometry;" << fib_names.at(i) << ";" << img_names.at(i) << ";" << boost::lexical_cast<std::string>(r) << ";" << boost::lexical_cast<std::string>(output.get(r,c)) << std::endl;
      }

      {
        auto output = mitk::Tractometry::NearestCentroidPointTractometry(input_images.at(i), fib->GetDeepCopy(), num_parcels, 1, 99, ref);
        for (unsigned int r=0; r<output.size(); ++r)
          for (unsigned int c=0; c<output.at(r).size(); ++c)
            statistics_file << "NearestCentroidPointTractometry;" << fib_names.at(i) << ";" << img_names.at(i) << ";" << boost::lexical_cast<std::string>(r) << ";" << boost::lexical_cast<std::string>(output.at(r).get(c)) << std::endl;
      }
      {
        itk::TractParcellationFilter< itk::Image<unsigned char, 3>, itk::Image<float, 3> >::Pointer parcellator = itk::TractParcellationFilter< itk::Image<unsigned char, 3>, itk::Image<float, 3> >::New();
        parcellator->SetInputImage(input_images.at(i));
        parcellator->SetNumParcels(num_parcels);
        parcellator->SetInputTract(fib->GetDeepCopy());
        parcellator->SetReferenceTract(ref);
        parcellator->Update();
        itk::Image<unsigned char, 3>::Pointer out_image_pp = parcellator->GetOutput(1);

        itk::ImageRegionConstIterator< itk::Image<float, 3> > it(input_images.at(i), input_images.at(i)->GetLargestPossibleRegion());
        itk::ImageRegionConstIterator< itk::Image<unsigned char, 3> > it2(out_image_pp, out_image_pp->GetLargestPossibleRegion());
        while (!it.IsAtEnd()) {
          if (it2.Get()>0)
            statistics_file << "TractParcellationFilter;" << fib_names.at(i) << ";" << img_names.at(i) << ";" << boost::lexical_cast<std::string>(it2.Get()-1) << ";" << boost::lexical_cast<std::string>(it.Get()) << std::endl;
          ++it;
          ++it2;
        }
      }
    }

    statistics_file.close();
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
