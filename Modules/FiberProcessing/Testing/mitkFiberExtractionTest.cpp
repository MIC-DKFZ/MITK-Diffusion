/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <itkFiberExtractionFilter.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

/** Fiber extraction with planar figures and ROI images on fiberBundleX.fib, subtraction and
 *  addition of bundles, against the MITK-Data references. */
class mitkFiberExtractionTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberExtractionTestSuite);
  MITK_TEST(PlanarFigures_ExtractAsReference);
  MITK_TEST(SubtractAndAdd_AsReference);
  MITK_TEST(RoiImage_PassingMask);
  MITK_TEST(RoiImage_EndingInMask);
  CPPUNIT_TEST_SUITE_END();

  using UcharImageType = itk::Image<unsigned char, 3>;
  using FilterType = itk::FiberExtractionFilter<unsigned char>;

  mitk::FiberBundle::Pointer Load(const std::string& file)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("" + file));
  }

  /** ROI1 OR (ROI2 AND ROI3) */
  mitk::FiberBundle::Pointer ExtractWithPlanarFigures(const mitk::FiberBundle* fib)
  {
    auto storage = mitk::StandaloneDataStorage::New();
    auto figure = [this](const std::string& file) {
      auto node = mitk::DataNode::New();
      node->SetData(mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath("" + file)).at(0));
      return node;
    };

    auto orComposite = mitk::PlanarFigureComposite::New();
    orComposite->setOperationType(mitk::PlanarFigureComposite::OR);
    auto orNode = mitk::DataNode::New();
    orNode->SetData(orComposite);
    auto orParent = mitk::DataStorage::SetOfObjects::New();
    orParent->push_back(orNode);

    auto andComposite = mitk::PlanarFigureComposite::New();
    andComposite->setOperationType(mitk::PlanarFigureComposite::AND);
    auto andNode = mitk::DataNode::New();
    andNode->SetData(andComposite);
    auto andParent = mitk::DataStorage::SetOfObjects::New();
    andParent->push_back(andNode);

    storage->Add(orNode);
    storage->Add(figure("ROI1.pf"), orParent);
    storage->Add(andNode, orParent);
    storage->Add(figure("ROI2.pf"), andParent);
    storage->Add(figure("ROI3.pf"), andParent);

    return fib->ExtractFiberSubset(orNode, storage);
  }

  UcharImageType::Pointer LoadRoiImage()
  {
    auto image = UcharImageType::New();
    mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ROIIMAGE.nrrd")), image);
    return image;
  }

  static void AssertEqualsReference(mitk::FiberBundle* result, mitk::FiberBundle* reference, const std::string& name)
  {
    const bool equal = result->Equals(reference);
    if (!equal)
      mitk::IOUtil::Save(result, mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + name + ".fib");
    CPPUNIT_ASSERT_MESSAGE(name + " equals the reference", equal);
  }

public:
  void PlanarFigures_ExtractAsReference()
  {
    auto extracted = ExtractWithPlanarFigures(Load("fiberBundleX.fib"));
    AssertEqualsReference(extracted, Load("fiberBundleX_extracted.fib"), "extracted");
  }

  void SubtractAndAdd_AsReference()
  {
    auto fib = Load("fiberBundleX.fib");
    auto extracted = ExtractWithPlanarFigures(fib);
    auto notExtracted = fib->SubtractBundle(extracted);
    AssertEqualsReference(notExtracted, Load("fiberBundleX_subtracted.fib"), "subtracted");
    AssertEqualsReference(extracted->AddBundle(notExtracted), Load("fiberBundleX_added.fib"), "added");
  }

  void RoiImage_PassingMask()
  {
    auto roi = LoadRoiImage();
    auto extractor = FilterType::New();
    extractor->SetInputFiberBundle(Load("fiberBundleX.fib"));
    extractor->SetRoiImages({roi});
    extractor->SetOverlapFraction(0.0);
    extractor->SetMode(FilterType::MODE::OVERLAP);
    extractor->Update();
    AssertEqualsReference(extractor->GetPositives().at(0), Load("fiberBundleX_passing-mask_v2.fib"), "passing_mask");
  }

  void RoiImage_EndingInMask()
  {
    auto roi = LoadRoiImage();
    auto extractor = FilterType::New();
    extractor->SetInputFiberBundle(Load("fiberBundleX.fib"));
    extractor->SetRoiImages({roi});
    extractor->SetOverlapFraction(0.0);
    extractor->SetBothEnds(true);
    extractor->SetMode(FilterType::MODE::ENDPOINTS);
    extractor->Update();
    AssertEqualsReference(extractor->GetPositives().at(0), Load("fiberBundleX_ending-in-mask.fib"), "ending_in_mask");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberExtraction)
