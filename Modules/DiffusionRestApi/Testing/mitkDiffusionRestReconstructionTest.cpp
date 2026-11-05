/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestTestUtil.h"

#include <mitkDiffusionIndices.h>
#include <mitkDiffusionPeaks.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkDiffusionTestData.h>
#include <mitkIOUtil.h>
#include <mitkLevelWindowProperty.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkTestingMacros.h>

/**
 * \brief The reconstruction routes give the results of the reconstruction functions, as nodes
 * below their input, and answer wrong input with 400.
 */
class mitkDiffusionRestReconstructionTestSuite : public DiffusionRestTestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestReconstructionTestSuite);
  MITK_TEST(Tensors_SameAsDirectCall);
  MITK_TEST(TensorIndices_SameAsDirectCall);
  MITK_TEST(ColorFa_RgbImageOnTheTensorGrid);
  MITK_TEST(Qball_SameAsDirectCall_WithAutomaticOrder);
  MITK_TEST(Peaks_SameAsDirectCall);
  MITK_TEST(WrongInput_400_UnknownUid_404);
  CPPUNIT_TEST_SUITE_END();

  mitk::Image::Pointer m_Dwi;
  mitk::DataNode::Pointer m_DwiNode;

public:
  void setUp() override
  {
    SetUpServer();
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"), &functor);
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(m_Dwi));
    m_DwiNode = AddNode(m_Dwi, "dwi");
  }

  void tearDown() override
  {
    m_DwiNode = nullptr;
    m_Dwi = nullptr;
    TearDownServer();
  }

  mitk::DataNode::Pointer ReconstructTensorsByRoute()
  {
    const auto result = JobResult("POST", "/diffusion/reconstruction/tensors", { { "dwi", Uid(m_DwiNode) } });
    return NodeOf(result.at("tensors"));
  }

  void Tensors_SameAsDirectCall()
  {
    auto node = ReconstructTensorsByRoute();
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_LinearLeastSquares_DT"), node->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(node)->front() == m_DwiNode);

    const auto expected = mitk::ReconstructTensors(m_Dwi, mitk::TensorFit::LinearLeastSquares, 0, false);
    CPPUNIT_ASSERT(dynamic_cast<mitk::TensorImage*>(node->GetData()) != nullptr);
    CPPUNIT_ASSERT(SameImage(expected, dynamic_cast<mitk::Image*>(node->GetData())));

    // The other fit with its options.
    const auto corrected = JobResult("POST", "/diffusion/reconstruction/tensors",
      { { "dwi", Uid(m_DwiNode) }, { "fit", "eigenvalue_correction" }, { "b0_threshold", 10 } });
    auto correctedNode = NodeOf(corrected.at("tensors"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_EigenvalueCorrected_DT"), correctedNode->GetName());
    const auto expectedCorrected = mitk::ReconstructTensors(m_Dwi, mitk::TensorFit::EigenvalueCorrection, 10, true);
    CPPUNIT_ASSERT(SameImage(expectedCorrected, dynamic_cast<mitk::Image*>(correctedNode->GetData())));
  }

  void TensorIndices_SameAsDirectCall()
  {
    auto tensorNode = ReconstructTensorsByRoute();
    auto* tensors = dynamic_cast<mitk::TensorImage*>(tensorNode->GetData());

    const auto result = JobResult("POST", "/diffusion/reconstruction/tensor-indices", { { "tensors", Uid(tensorNode) }, { "indices", { "FA", "MD" } } });
    const auto images = result.at("images");
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), images.size());
    CPPUNIT_ASSERT_EQUAL(std::string("FA"), images[0].at("index").get<std::string>());

    auto faNode = NodeOf(images[0]);
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_LinearLeastSquares_DT_FA"), faNode->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(faNode)->front() == tensorNode);
    CPPUNIT_ASSERT(SameImage(mitk::TensorIndexImage(tensors, mitk::TensorIndex::FA, 1.0), dynamic_cast<mitk::Image*>(faNode->GetData())));
    CPPUNIT_ASSERT(SameImage(mitk::TensorIndexImage(tensors, mitk::TensorIndex::MD, 1.0), dynamic_cast<mitk::Image*>(NodeOf(images[1])->GetData())));

    mitk::LevelWindow levelWindow;
    CPPUNIT_ASSERT(faNode->GetLevelWindow(levelWindow));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, levelWindow.GetLevel(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, levelWindow.GetWindow(), 1e-6);

    // The default is the FA alone.
    const auto fa = JobResult("POST", "/diffusion/reconstruction/tensor-indices", { { "tensors", Uid(tensorNode) } });
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), fa.at("images").size());
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/tensor-indices", { { "tensors", Uid(tensorNode) }, { "indices", { "XY" } } }));
  }

  void ColorFa_RgbImageOnTheTensorGrid()
  {
    auto tensorNode = ReconstructTensorsByRoute();
    const auto result = JobResult("POST", "/diffusion/reconstruction/color-fa", { { "tensors", Uid(tensorNode) } });
    auto node = NodeOf(result.at("image"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_LinearLeastSquares_DT_ColorFA"), node->GetName());

    auto* image = dynamic_cast<mitk::Image*>(node->GetData());
    CPPUNIT_ASSERT(image != nullptr);
    CPPUNIT_ASSERT_EQUAL(4u, static_cast<unsigned int>(image->GetPixelType().GetNumberOfComponents()));
    CPPUNIT_ASSERT_EQUAL(3u, image->GetDimension());
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_EQUAL(m_Dwi->GetDimension(i), image->GetDimension(i));
  }

  void Qball_SameAsDirectCall_WithAutomaticOrder()
  {
    const auto result = JobResult("POST", "/diffusion/reconstruction/qball", { { "dwi", Uid(m_DwiNode) }, { "sampled_odf", true } });
    auto shNode = NodeOf(result.at("sh"));
    auto odfNode = NodeOf(result.at("odf"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_CSA_Qball_SH"), shNode->GetName());
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_CSA_Qball_SH_Sampled"), odfNode->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(odfNode)->front() == shNode);

    // The order follows the number of diffusion-weighted volumes: the largest up to 8 whose
    // coefficients they determine.
    const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(m_Dwi);
    const auto weighted = gradients->Size() - mitk::DiffusionPropertyHelper::GetB0Indices(gradients, mitk::DiffusionPropertyHelper::GetReferenceBValue(m_Dwi)).size();
    unsigned int expectedOrder = 2;
    for (unsigned int order = 8; order >= 2; order -= 2)
    {
      if ((order + 1) * (order + 2) / 2 <= weighted)
      {
        expectedOrder = order;
        break;
      }
    }
    CPPUNIT_ASSERT_EQUAL(expectedOrder, result.at("sh_order").get<unsigned int>());

    mitk::QballParameters parameters;
    parameters.method = mitk::QballMethod::SolidAngle;
    parameters.shOrder = expectedOrder;
    parameters.sampledOdf = true;
    const auto expected = mitk::ReconstructQball(m_Dwi, parameters);
    CPPUNIT_ASSERT(SameImage(expected.sh, dynamic_cast<mitk::Image*>(shNode->GetData())));
    CPPUNIT_ASSERT(SameImage(expected.odf, dynamic_cast<mitk::Image*>(odfNode->GetData())));

    // Without sampled ODFs, and with an explicit order.
    const auto plain = JobResult("POST", "/diffusion/reconstruction/qball", { { "dwi", Uid(m_DwiNode) }, { "method", "standard" }, { "sh_order", 4 } });
    CPPUNIT_ASSERT(plain.at("odf").is_null());
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_Qball_SH"), NodeOf(plain.at("sh"))->GetName());
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/qball", { { "dwi", Uid(m_DwiNode) }, { "sh_order", 5 } }));
  }

  void Peaks_SameAsDirectCall()
  {
    const auto qball = JobResult("POST", "/diffusion/reconstruction/qball", { { "dwi", Uid(m_DwiNode) } });
    auto shNode = NodeOf(qball.at("sh"));
    auto* sh = dynamic_cast<mitk::ShImage*>(shNode->GetData());

    const auto result = JobResult("POST", "/diffusion/reconstruction/peaks", { { "image", Uid(shNode) } });
    auto peakNode = NodeOf(result.at("peaks"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_CSA_Qball_SH_PEAKS"), peakNode->GetName());
    CPPUNIT_ASSERT(dynamic_cast<mitk::PeakImage*>(peakNode->GetData()) != nullptr);

    const auto expected = mitk::ExtractPeaks(sh, mitk::PeakParameters());
    CPPUNIT_ASSERT(SameImage(expected.peaks, dynamic_cast<mitk::Image*>(peakNode->GetData())));

    // Tensor peaks: the principal directions.
    auto tensorNode = ReconstructTensorsByRoute();
    const auto tensorPeaks = JobResult("POST", "/diffusion/reconstruction/peaks", { { "image", Uid(tensorNode) }, { "fa_threshold", 0.1 } });
    auto principal = NodeOf(tensorPeaks.at("peaks"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_LinearLeastSquares_DT_PrincipalDirection"), principal->GetName());
    const auto expectedPrincipal = mitk::ExtractTensorPeaks(dynamic_cast<mitk::TensorImage*>(tensorNode->GetData()), 0.1, true);
    CPPUNIT_ASSERT(SameImage(expectedPrincipal.peaks, dynamic_cast<mitk::Image*>(principal->GetData())));
  }

  void WrongInput_400_UnknownUid_404()
  {
    auto tensorNode = ReconstructTensorsByRoute();

    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/tensors", { { "dwi", Uid(tensorNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/tensors", nlohmann::json::object()));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("POST", "/diffusion/reconstruction/tensors", { { "dwi", "node_999999" } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/tensors", { { "dwi", Uid(m_DwiNode) }, { "fit", "magic" } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/tensor-indices", { { "tensors", Uid(m_DwiNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/color-fa", { { "tensors", Uid(m_DwiNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/qball", { { "dwi", Uid(tensorNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/peaks", { { "image", Uid(m_DwiNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/reconstruction/peaks", { { "image", Uid(tensorNode) }, { "mask", Uid(m_DwiNode) } }));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestReconstruction)
