/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkPeakImage.h>
#include <mitkTestFixture.h>

#include <vtkCell.h>

using namespace mitk::CmdAppTest;

/** Smoke tests of the apps in Modules/DiffusionCmdApps/FiberProcessing, on fiberBundleX.fib (196 fibers). */
class mitkDiffusionCmdAppsFiberProcessingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsFiberProcessingTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(FiberModification_LengthThenResampling);
  MITK_TEST(FiberModification_NegativeTranslation);
  MITK_TEST(FiberModification_MissingInput_Fails);
  MITK_TEST(FiberJoin_TwiceTheFibers);
  MITK_TEST(FiberExtraction_PlanarFigure);
  MITK_TEST(FiberExtractionRoi_DottedOutputName);
  MITK_TEST(FiberClustering_ClustersAndCentroids);
  MITK_TEST(FiberColoring_AllFibers);
  MITK_TEST(FiberDirectionExtraction_Peaks);
  MITK_TEST(TractDensityFilter_Subset);
  MITK_TEST(Sift2WeightCopy_Weights);
  MITK_TEST(FitFibersToImage_FittedBundles);
  CPPUNIT_TEST_SUITE_END();

  static std::string Fibers() { return Data("fiberBundleX.fib"); }
  static constexpr unsigned int NumInputFibers = 196;

public:
  void HelpAndXml()
  {
    CheckHelpAndXml("FiberProcessing", {"Sift2WeightCopy", "FiberExtraction", "FiberExtractionRoi", "FiberModification",
                                        "FitFibersToImage", "FiberDirectionExtraction", "FiberJoin", "FiberClustering",
                                        "TractDensityFilter", "FiberColoring"});
  }

  // The length filter runs before the resampling, which would change the lengths.
  void FiberModification_LengthThenResampling()
  {
    const auto dir = WorkDir("FiberModification");
    RunOk("FiberModification", {"-i", Fibers(), "-o", "modified.fib", "--min_length", "20", "--num_resampling", "10"}, dir);
    auto fibers = Load<mitk::FiberBundle>(dir, "modified.fib");
    CPPUNIT_ASSERT(fibers->GetNumFibers() > 0 && fibers->GetNumFibers() < NumInputFibers);
    for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
    {
      CPPUNIT_ASSERT_EQUAL(vtkIdType(10), fibers->GetFiberPolyData()->GetCell(i)->GetNumberOfPoints());
      CPPUNIT_ASSERT(fibers->GetFiberLength(i) >= 20 * 0.9);
    }
  }

  void FiberModification_NegativeTranslation()
  {
    const auto dir = WorkDir("FiberModificationTranslate");
    RunOk("FiberModification", {"-i", Fibers(), "-o", "moved.fib", "--translate_x", "-5"}, dir);
    auto input = mitk::IOUtil::Load<mitk::FiberBundle>(Fibers());
    auto moved = Load<mitk::FiberBundle>(dir, "moved.fib");
    CPPUNIT_ASSERT_EQUAL(input->GetNumFibers(), moved->GetNumFibers());
    const double *before = input->GetFiberPolyData()->GetPoint(0);
    const double *after = moved->GetFiberPolyData()->GetPoint(0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.0, after[0] - before[0], 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, after[1] - before[1], 1e-4);
  }

  void FiberModification_MissingInput_Fails()
  {
    const auto dir = WorkDir("FiberModificationFails");
    RunFails("FiberModification", {"-i", Path(dir, "missing.fib"), "-o", "out.fib", "--min_length", "20"}, dir);
  }

  void FiberJoin_TwiceTheFibers()
  {
    const auto dir = WorkDir("FiberJoin");
    RunOk("FiberJoin", {"-i", Fibers(), Fibers(), "-o", "joined.fib"}, dir);
    CPPUNIT_ASSERT_EQUAL(2 * NumInputFibers, NumFibers(dir, "joined.fib"));
  }

  void FiberExtraction_PlanarFigure()
  {
    const auto dir = WorkDir("FiberExtraction");
    RunOk("FiberExtraction", {"-i", Fibers(), "-o", "extracted.fib", "-pf1", Data("ROI1.pf")}, dir);
    const auto count = NumFibers(dir, "extracted.fib");
    CPPUNIT_ASSERT(count > 0 && count < NumInputFibers);
  }

  // Without labels, the output is <root><extension> (#43: the root keeps its dot).
  void FiberExtractionRoi_DottedOutputName()
  {
    const auto dir = WorkDir("FiberExtractionRoi");
    RunOk("FiberExtractionRoi", {"-i", Fibers(), "-o", "roi.01", "--rois", Data("ROIIMAGE.nrrd"), "--overlap_fraction", "0.1"}, dir);
    const auto count = NumFibers(dir, "roi.01.trk");
    CPPUNIT_ASSERT(count > 0 && count <= NumInputFibers);
  }

  void FiberClustering_ClustersAndCentroids()
  {
    const auto dir = WorkDir("FiberClustering");
    RunOk("FiberClustering", {"-i", Fibers(), "-o", "clusters_", "--cluster_size", "15", "--output_centroids", "--merge_centroids"}, dir);
    unsigned int clustered = 0;
    unsigned int clusters = 0;
    for (; std::filesystem::exists(Path(dir, "clusters_Cluster_" + std::to_string(clusters) + ".fib")); ++clusters)
      clustered += NumFibers(dir, "clusters_Cluster_" + std::to_string(clusters) + ".fib");
    if (std::filesystem::exists(Path(dir, "clusters_Cluster_unassigned.fib")))
      clustered += NumFibers(dir, "clusters_Cluster_unassigned.fib");
    CPPUNIT_ASSERT(clusters > 1);
    CPPUNIT_ASSERT_EQUAL(NumInputFibers, clustered);
    CPPUNIT_ASSERT_EQUAL(clusters, NumFibers(dir, "clusters_fiberBundleX_Centroids.fib"));
  }

  void FiberColoring_AllFibers()
  {
    const auto dir = WorkDir("FiberColoring");
    RunOk("FiberColoring", {"-i", Fibers(), "-o", "colored.fib"}, dir);
    CPPUNIT_ASSERT_EQUAL(NumInputFibers, NumFibers(dir, "colored.fib"));
  }

  void FiberDirectionExtraction_Peaks()
  {
    const auto dir = WorkDir("FiberDirectionExtraction");
    RunOk("FiberDirectionExtraction", {"-i", Fibers(), "-o", "fibers"}, dir);
    auto peaks = Load<mitk::PeakImage>(dir, "fibers_DIRECTIONS.nii.gz", {"Peak Image"});
    CPPUNIT_ASSERT_EQUAL(9u, peaks->GetDimension(3));
  }

  void TractDensityFilter_Subset()
  {
    const auto dir = WorkDir("TractDensityFilter");
    RunOk("TractDensityFilter", {"-i", Fibers(), "-o", "dense.fib"}, dir);
    const auto count = NumFibers(dir, "dense.fib");
    CPPUNIT_ASSERT(count > 0 && count <= NumInputFibers);
  }

  void Sift2WeightCopy_Weights()
  {
    const auto dir = WorkDir("Sift2WeightCopy");
    {
      std::ofstream weights(Path(dir, "weights.txt"));
      for (unsigned int i = 0; i < NumInputFibers; ++i)
        weights << (i % 2 == 0 ? 0.5 : 2.0) << " ";
    }
    RunOk("Sift2WeightCopy", {"-i", Fibers(), "-o", "weighted.fib", "--weights", "weights.txt"}, dir);
    auto fibers = Load<mitk::FiberBundle>(dir, "weighted.fib");
    CPPUNIT_ASSERT_EQUAL(NumInputFibers, fibers->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, fibers->GetFiberWeight(0), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, fibers->GetFiberWeight(1), 1e-6);
  }

  void FitFibersToImage_FittedBundles()
  {
    const auto dir = WorkDir("FitFibersToImage");
    RunOk("FitFibersToImage", {"-i1", Data("FiberFit/Cluster_0.fib"), Data("FiberFit/Cluster_1.fib"), "-i2",
                               Data("FiberFit/csd_peak_image.nii.gz"), "-o", "fit_", "--max_iter", "2"}, dir);
    CPPUNIT_ASSERT(NumFibers(dir, "fit_Cluster_0_fitted.fib") > 0);
    CPPUNIT_ASSERT(NumFibers(dir, "fit_Cluster_1_fitted.fib") > 0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsFiberProcessing)
