/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestTestUtil.h"

#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkLog.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionIteratorWithIndex.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>

#include <cstdlib>
#include <filesystem>
#include <functional>

/**
 * \brief The extraction route keeps the fibers that pass, end in or avoid synthetic box ROIs as
 * the fiber extraction filter and MitkFiberExtractionRoi do; the metrics route counts fibers,
 * lengths, volume and the fibers per ROI.
 */
class mitkDiffusionRestExtractionTestSuite : public DiffusionRestTestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestExtractionTestSuite);
  MITK_TEST(Extract_PassEndNot_KnownFibers);
  MITK_TEST(Extract_SameAsTheApp);
  MITK_TEST(Metrics_FibersLengthsVolumeRois);
  MITK_TEST(WrongInput_400);
  CPPUNIT_TEST_SUITE_END();

  static constexpr unsigned int Size = 40;

  mitk::FiberBundle::Pointer m_Bundle;
  mitk::DataNode::Pointer m_BundleNode;
  mitk::DataNode::Pointer m_StartBox;   ///< x in [0, 5]: the x fibers start here
  mitk::DataNode::Pointer m_MiddleBox;  ///< x in [18, 22]: the x fibers pass, the y fibers (at x = 25) do not
  mitk::DataNode::Pointer m_YStartBox;  ///< y in [0, 5]: the y fibers start here

  /** Nine fibers along x (35 mm) and four along y (35 mm), one point per millimeter. */
  static mitk::FiberBundle::Pointer MakeBundle()
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto cells = vtkSmartPointer<vtkCellArray>::New();

    const auto addFiber = [&](const std::function<std::array<double, 3>(int)>& point) {
      auto line = vtkSmartPointer<vtkPolyLine>::New();
      for (int i = 0; i <= 35; ++i)
      {
        const auto p = point(i);
        line->GetPointIds()->InsertNextId(points->InsertNextPoint(p[0], p[1], p[2]));
      }
      cells->InsertNextCell(line);
    };

    for (double y : { 10.0, 15.0, 20.0 })
      for (double z : { 10.0, 15.0, 20.0 })
        addFiber([y, z](int i) { return std::array<double, 3>{ 2.0 + i, y, z }; });

    for (double z : { 10.0, 15.0, 20.0, 25.0 })
      addFiber([z](int i) { return std::array<double, 3>{ 25.0, 2.0 + i, z }; });

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);
    return mitk::FiberBundle::New(polyData);
  }

  /** An unsigned char image of Size^3 voxels at 1 mm, 1 inside the box. */
  static mitk::Image::Pointer MakeBox(const std::function<bool(const itk::Index<3>&)>& inside)
  {
    using ImageType = itk::Image<unsigned char, 3>;
    auto image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize({{Size, Size, Size}});
    image->SetRegions(region);
    image->Allocate();
    itk::ImageRegionIteratorWithIndex<ImageType> it(image, region);
    for (; !it.IsAtEnd(); ++it)
      it.Set(inside(it.GetIndex()) ? 1 : 0);
    return mitk::GrabItkImageMemory(image.GetPointer());
  }

  static mitk::FiberBundle* Fibers(mitk::DataNode* node)
  {
    auto* fibers = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    CPPUNIT_ASSERT(fibers != nullptr);
    return fibers;
  }

public:
  void setUp() override
  {
    SetUpServer();
    m_Bundle = MakeBundle();
    CPPUNIT_ASSERT_EQUAL(13u, m_Bundle->GetNumFibers());
    m_BundleNode = AddNode(m_Bundle, "fibers");
    m_StartBox = AddNode(MakeBox([](const itk::Index<3>& i) { return i[0] <= 5; }), "start");
    m_MiddleBox = AddNode(MakeBox([](const itk::Index<3>& i) { return i[0] >= 18 && i[0] <= 22; }), "middle");
    m_YStartBox = AddNode(MakeBox([](const itk::Index<3>& i) { return i[1] <= 5; }), "ystart");
  }

  void tearDown() override
  {
    m_BundleNode = m_StartBox = m_MiddleBox = m_YStartBox = nullptr;
    m_Bundle = nullptr;
    TearDownServer();
  }

  unsigned int ExtractCount(const nlohmann::json& rois, std::string* nodeName = nullptr)
  {
    auto body = rois;
    body["tractogram"] = Uid(m_BundleNode);
    const auto result = JobResult("POST", "/diffusion/tractography/extract", body);
    CPPUNIT_ASSERT_EQUAL(13u, result.at("fibers_in").get<unsigned int>());
    auto node = NodeOf(result.at("tractogram"));
    CPPUNIT_ASSERT(m_DataStorage->GetSources(node)->front() == m_BundleNode);
    if (nodeName != nullptr)
      *nodeName = node->GetName();
    CPPUNIT_ASSERT_EQUAL(result.at("fibers").get<unsigned int>(), Fibers(node)->GetNumFibers());
    return result.at("fibers").get<unsigned int>();
  }

  void Extract_PassEndNot_KnownFibers()
  {
    std::string name;
    CPPUNIT_ASSERT_EQUAL(9u, ExtractCount({ { "pass", { Uid(m_MiddleBox) } } }, &name));
    CPPUNIT_ASSERT_EQUAL(std::string("fibers_extracted"), name);
    CPPUNIT_ASSERT_EQUAL(9u, ExtractCount({ { "end", { Uid(m_StartBox) } } }));
    CPPUNIT_ASSERT_EQUAL(0u, ExtractCount({ { "end", { Uid(m_StartBox) } }, { "both_ends", true } }));
    CPPUNIT_ASSERT_EQUAL(0u, ExtractCount({ { "end", { Uid(m_StartBox), Uid(m_YStartBox) } } }));
    CPPUNIT_ASSERT_EQUAL(4u, ExtractCount({ { "not", { Uid(m_MiddleBox) } } }));
    CPPUNIT_ASSERT_EQUAL(0u, ExtractCount({ { "pass", { Uid(m_MiddleBox) } }, { "not", { Uid(m_StartBox) } } }));
    CPPUNIT_ASSERT_EQUAL(4u, ExtractCount({ { "end", { { { "mask", Uid(m_YStartBox) }, { "name", "y start" } } } } }));
    // The middle box holds 5 of 35 mm of each x fiber: a larger required overlap keeps none.
    CPPUNIT_ASSERT_EQUAL(0u, ExtractCount({ { "pass", { Uid(m_MiddleBox) } }, { "overlap_fraction", 0.5 } }));
    CPPUNIT_ASSERT_EQUAL(0u, ExtractCount({ { "pass", { Uid(m_MiddleBox) } }, { "min_fibers", 10 } }));
    CPPUNIT_ASSERT_EQUAL(std::string("CST"), NodeOf(JobResult("POST", "/diffusion/tractography/extract",
      { { "tractogram", Uid(m_BundleNode) }, { "pass", { Uid(m_MiddleBox) } }, { "name", "CST" } }).at("tractogram"))->GetName());
  }

  void Extract_SameAsTheApp()
  {
    const char* app = std::getenv("MITK_DIFFUSION_EXTRACTION_APP");

    if (app == nullptr || !std::filesystem::exists(app))
    {
      MITK_WARN << "MitkFiberExtractionRoi is not built (MITK_DIFFUSION_EXTRACTION_APP); the comparison with the app is left out.";
      return;
    }

    const auto dir = std::filesystem::temp_directory_path() / "mitkDiffusionRestExtractionTest";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const auto input = (dir / "fibers.fib").string();
    const auto middle = (dir / "middle.nii.gz").string();
    const auto start = (dir / "start.nii.gz").string();
    mitk::IOUtil::Save(m_Bundle, input);
    mitk::IOUtil::Save(m_MiddleBox->GetData(), middle);
    mitk::IOUtil::Save(m_StartBox->GetData(), start);

    const auto run = [&](const std::string& output, const std::string& arguments) {
      const auto command = "\"" + std::string(app) + "\" -i \"" + input + "\" -o \"" + output + "\" " + arguments + " > /dev/null 2>&1";
      CPPUNIT_ASSERT_EQUAL_MESSAGE(command, 0, std::system(command.c_str()));
      return mitk::IOUtil::Load<mitk::FiberBundle>(output);
    };

    // Passing fibers: any part inside (overlap fraction 0) as the route's "pass".
    const auto appPass = run((dir / "pass.fib").string(), "--rois \"" + middle + "\" --overlap_fraction 0");
    const auto routePass = JobResult("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) }, { "pass", { Uid(m_MiddleBox) } } });
    CPPUNIT_ASSERT(appPass->Equals(Fibers(NodeOf(routePass.at("tractogram")))));

    // Ending fibers: one end inside, as the route's "end".
    const auto appEnd = run((dir / "end.fib").string(), "--rois \"" + start + "\"");
    const auto routeEnd = JobResult("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) }, { "end", { Uid(m_StartBox) } } });
    CPPUNIT_ASSERT(appEnd->Equals(Fibers(NodeOf(routeEnd.at("tractogram")))));
    CPPUNIT_ASSERT_EQUAL(9u, appEnd->GetNumFibers());

    std::filesystem::remove_all(dir);
  }

  void Metrics_FibersLengthsVolumeRois()
  {
    const auto result = JobResult("POST", "/diffusion/tractography/metrics", {
      { "tractogram", Uid(m_BundleNode) },
      { "rois", { Uid(m_MiddleBox), { { "mask", Uid(m_StartBox) }, { "name", "start box" } } } } });

    CPPUNIT_ASSERT_EQUAL(13u, result.at("fibers").get<unsigned int>());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, result.at("length_mm").at("min").get<double>(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, result.at("length_mm").at("max").get<double>(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0, result.at("length_mm").at("mean").get<double>(), 1e-3);
    // 13 fibers of 35 mm on a 1 mm grid, one voxel per millimeter of fiber, two crossing voxels shared.
    const double volume = result.at("volume_ml").get<double>();
    CPPUNIT_ASSERT(volume > 0.40 && volume < 0.48);

    const auto rois = result.at("rois");
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), rois.size());
    CPPUNIT_ASSERT_EQUAL(Uid(m_MiddleBox), rois[0].at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(9u, rois[0].at("fibers_passing").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(0u, rois[0].at("fibers_ending").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(std::string("start box"), rois[1].at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(9u, rois[1].at("fibers_passing").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(9u, rois[1].at("fibers_ending").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(0u, rois[1].at("fibers_both_ends").get<unsigned int>());
    // Nothing was added.
    CPPUNIT_ASSERT_EQUAL(0u, m_DataStorage->GetDerivations(m_BundleNode)->Size());
  }

  void WrongInput_400()
  {
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_StartBox) }, { "pass", { Uid(m_MiddleBox) } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) }, { "pass", Uid(m_MiddleBox) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) }, { "pass", { Uid(m_BundleNode) } } }));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("POST", "/diffusion/tractography/extract", { { "tractogram", Uid(m_BundleNode) }, { "pass", { "node_999999" } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/metrics", { { "tractogram", Uid(m_StartBox) } }));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("POST", "/diffusion/tractography/metrics", { { "tractogram", "node_999999" } }));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestExtraction)
