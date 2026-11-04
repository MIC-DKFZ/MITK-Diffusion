/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkFiberBundle.h>
#include <mitkFiberBundleIOHelper.h>
#include <mitkImageCast.h>
#include <mitkCommandLineParser.h>
#include <mitkLexicalCast.h>
#include <mitkIOUtil.h>
#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricEuclideanMean.h>
#include <mitkClusteringMetricEuclideanMax.h>
#include <mitkClusteringMetricEuclideanStd.h>
#include <mitkClusteringMetricAnatomic.h>
#include <mitkClusteringMetricScalarMap.h>
#include <mitkClusteringMetricInnerAngles.h>
#include <mitkClusteringMetricLength.h>
#include <itksys/SystemTools.hxx>
#include <mitkDiffusionDataIOHelper.h>
#include <mitkExceptionMacro.h>

/*!
\brief Spatially cluster fibers
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("Fiber Clustering");
  parser.setCategory("Fiber Processing");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input:", "input fiber bundle (.fib; .trk; .tck; .trx)", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output root", us::Any(), false, false, false, mitkCommandLineParser::Output);

  parser.addArgument("cluster_size", "", mitkCommandLineParser::Int, "Cluster size:", "", 10);
  parser.addArgument("fiber_points", "", mitkCommandLineParser::Int, "Fiber points:", "", 12);
  parser.addArgument("min_fibers", "", mitkCommandLineParser::Int, "Min. fibers per cluster:", "", 1);
  parser.addArgument("max_clusters", "", mitkCommandLineParser::Int, "Max. clusters:", "");
  parser.addArgument("merge_clusters", "", mitkCommandLineParser::Float, "Merge clusters:", "Set to 0 to avoid merging and to -1 to use the original cluster size", -1.0);
  parser.addArgument("output_centroids", "", mitkCommandLineParser::Bool, "Output centroids:", "");
  parser.addArgument("only_centroids", "", mitkCommandLineParser::Bool, "Output only centroids:", "");
  parser.addArgument("merge_centroids", "", mitkCommandLineParser::Bool, "Merge centroids:", "");
  parser.addArgument("metrics", "", mitkCommandLineParser::StringList, "Metrics:", "EU_MEAN; EU_STD; EU_MAX; ANGLES; ANAT; MAP; LENGTH", std::string("EU_MEAN"));
  parser.addArgument("metric_weights", "", mitkCommandLineParser::StringList, "Metric weights:", "add one float weight for each used metric");
  parser.addArgument("input_centroids", "", mitkCommandLineParser::String, "Input centroids:", "", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("scalar_map", "", mitkCommandLineParser::String, "Scalar map:", "", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("parcellation", "", mitkCommandLineParser::String, "Parcellation:", "", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("file_ending", "", mitkCommandLineParser::String, "File ending:", "");

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string out_root = us::any_cast<std::string>(parsedArgs["o"]);

  bool only_centroids = false;
  if (parsedArgs.count("only_centroids"))
    only_centroids = us::any_cast<bool>(parsedArgs["only_centroids"]);

  bool merge_centroids = false;
  if (parsedArgs.count("merge_centroids"))
    merge_centroids = us::any_cast<bool>(parsedArgs["merge_centroids"]);

  int cluster_size = 10;
  if (parsedArgs.count("cluster_size"))
    cluster_size = us::any_cast<int>(parsedArgs["cluster_size"]);

  int fiber_points = 12;
  if (parsedArgs.count("fiber_points"))
    fiber_points = us::any_cast<int>(parsedArgs["fiber_points"]);

  int min_fibers = 1;
  if (parsedArgs.count("min_fibers"))
    min_fibers = us::any_cast<int>(parsedArgs["min_fibers"]);

  int max_clusters = 0;
  if (parsedArgs.count("max_clusters"))
    max_clusters = us::any_cast<int>(parsedArgs["max_clusters"]);

  float merge_clusters = -1.0;
  if (parsedArgs.count("merge_clusters"))
    merge_clusters = us::any_cast<float>(parsedArgs["merge_clusters"]);

  bool output_centroids = false;
  if (parsedArgs.count("output_centroids"))
    output_centroids = us::any_cast<bool>(parsedArgs["output_centroids"]);

  std::vector< std::string > metric_strings = {"EU_MEAN"};
  if (parsedArgs.count("metrics"))
    metric_strings = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["metrics"]);

  std::vector< std::string > metric_weights = {"1.0"};
  if (parsedArgs.count("metric_weights"))
    metric_weights = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["metric_weights"]);

  std::string input_centroids = "";
  if (parsedArgs.count("input_centroids"))
    input_centroids = us::any_cast<std::string>(parsedArgs["input_centroids"]);

  std::string scalar_map = "";
  if (parsedArgs.count("scalar_map"))
    scalar_map = us::any_cast<std::string>(parsedArgs["scalar_map"]);

  std::string parcellation = "";
  if (parsedArgs.count("parcellation"))
    parcellation = us::any_cast<std::string>(parsedArgs["parcellation"]);

  std::string file_ending = ".fib";
  if (parsedArgs.count("file_ending"))
    file_ending = us::any_cast<std::string>(parsedArgs["file_ending"]);

  if (metric_strings.size()!=metric_weights.size())
  {
    MITK_INFO << "Each metric needs an associated metric weight!";
    return EXIT_FAILURE;
  }

  // A cluster size of 0 would never reach the largest cluster distance, and the metrics compare
  // fibers point by point.
  if (cluster_size<=0 || fiber_points<2)
  {
    MITK_INFO << "--cluster_size must be positive and --fiber_points at least 2.";
    return EXIT_FAILURE;
  }

  try
  {
    typedef itk::Image< float, 3 > FloatImageType;
    typedef itk::Image< short, 3 > ShortImageType;

    mitk::FiberBundle::Pointer fib = mitk::FiberBundleIOHelper::load_fib(inFileName);

    float max_d = 0;
    int i=1;
    std::vector< float > distances;
    while (max_d < fib->GetGeometry()->GetDiagonalLength()/2)
    {
      distances.push_back(cluster_size*i);
      max_d = cluster_size*i;
      ++i;
    }

    std::shared_ptr< mitk::TractClusteringFilter > clusterer = std::make_shared<mitk::TractClusteringFilter>();
    clusterer->SetDistances(distances);
    clusterer->SetTractogram(fib);

    if (input_centroids!="")
    {

      mitk::FiberBundle::Pointer in_centroids = mitk::FiberBundleIOHelper::load_fib(input_centroids);
      clusterer->SetInCentroids(in_centroids);
    }

    // Each metric name must create its metric, so that its weight is not given to another one.
    int mc = 0;
    for (auto m : metric_strings)
    {
      float w = boost::lexical_cast<float>(metric_weights.at(mc));
      MITK_INFO << "Metric: " << m << " (w=" << w << ")";
      std::unique_ptr<mitk::ClusteringMetric> metric;
      if (m=="EU_MEAN")
        metric = std::make_unique<mitk::ClusteringMetricEuclideanMean>();
      else if (m=="EU_STD")
        metric = std::make_unique<mitk::ClusteringMetricEuclideanStd>();
      else if (m=="EU_MAX")
        metric = std::make_unique<mitk::ClusteringMetricEuclideanMax>();
      else if (m=="ANGLES")
        metric = std::make_unique<mitk::ClusteringMetricInnerAngles>();
      else if (m=="LENGTH")
        metric = std::make_unique<mitk::ClusteringMetricLength>();
      else if (m=="MAP" || m=="ANAT")
      {
        const std::string& image_file = m=="MAP" ? scalar_map : parcellation;
        if (image_file.empty())
          mitkThrow() << "Metric " << m << " needs --" << (m=="MAP" ? "scalar_map" : "parcellation") << ".";
        mitk::Image::Pointer mitk_map = mitk::IOUtil::Load<mitk::Image>(image_file);
        if (mitk_map->GetDimension()!=3)
          mitkThrow() << "Metric " << m << " needs a 3D image: " << image_file;

        if (m=="MAP")
        {
          FloatImageType::Pointer itk_map = FloatImageType::New();
          mitk::CastToItkImage(mitk_map, itk_map);
          auto map_metric = std::make_unique<mitk::ClusteringMetricScalarMap>();
          map_metric->SetImages({itk_map});
          metric = std::move(map_metric);
        }
        else
        {
          ShortImageType::Pointer itk_map = ShortImageType::New();
          mitk::CastToItkImage(mitk_map, itk_map);
          auto anat_metric = std::make_unique<mitk::ClusteringMetricAnatomic>();
          anat_metric->SetParcellations({itk_map});
          metric = std::move(anat_metric);
        }
      }
      else
        mitkThrow() << "Unknown metric " << m << ".";

      metric->SetScale(w);
      clusterer->AddMetric(std::move(metric));
      mc++;
    }

    if (mc==0)
    {
      MITK_INFO << "No metric selected!";
      return EXIT_FAILURE;
    }

    clusterer->SetMergeDuplicateThreshold(merge_clusters);
    clusterer->SetNumPoints(fiber_points);
    clusterer->SetMaxClusters(max_clusters);
    clusterer->SetMinClusterSize(min_fibers);
    clusterer->Update();
    std::vector<mitk::FiberBundle::Pointer> tracts = clusterer->GetOutTractograms();
    std::vector<mitk::FiberBundle::Pointer> centroids = clusterer->GetOutCentroids();

    MITK_INFO << "Saving clusters";
    mitk::CoutSilencer silencer;

    if (!only_centroids)
    {
      for (unsigned int i=0; i<tracts.size(); ++i)
        mitk::IOUtil::Save(tracts.at(i), out_root + "Cluster_" + boost::lexical_cast<std::string>(i) + file_ending);
      if (clusterer->GetOutUnassigned()->GetNumFibers()>0)
        mitk::IOUtil::Save(clusterer->GetOutUnassigned(), out_root + "Cluster_unassigned" + file_ending);
    }


    if (output_centroids && !merge_centroids)
    {
      for (unsigned int i=0; i<centroids.size(); ++i)
        mitk::IOUtil::Save(centroids.at(i), out_root + "Centroid_" + boost::lexical_cast<std::string>(i) + file_ending);
    }
    else if (output_centroids)
    {
      mitk::FiberBundle::Pointer centroid = mitk::FiberBundle::New();
      centroid = centroid->AddBundles(centroids);
      mitk::IOUtil::Save(centroid, out_root + mitk::DiffusionDataIOHelper::GetNameWithoutExtension(inFileName) + "_Centroids" + file_ending);
    }
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
    std::cout << "ERROR!?!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
