/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkNodePredicateIsDWI.h>

class mitkDiffusionPropertyHelperTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionPropertyHelperTestSuite);
  MITK_TEST(GetMeasurementFrame_Missing_ReturnsIdentity);
  MITK_TEST(InitializeImage_NoMeasurementFrame_KeepsGradients);
  MITK_TEST(Getters_NoDiffusionProperties_Throw);
  MITK_TEST(AverageRedundantGradients_NegativeValues_AveragesSigned);
  MITK_TEST(InitializeImage_RotationAndFrame_DerivesOriginalGradientsExactly);
  MITK_TEST(InitializeImage_Flags_SelectRotationAndFrame);
  MITK_TEST(InitializeImage_Gradients_BValueMapByRoundedBValue);
  MITK_TEST(IsDiffusionWeightedImage_IncompleteImages_False);
  MITK_TEST(CopyProperties_IgnoreOriginal_RemovesOriginalGradients);
  MITK_TEST(CalcAveragedDirectionSet_AntipodalAndClose_Merged);
  MITK_TEST(NodePredicateIsDWI_Nodes_ChecksData);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;

  mitk::Image::Pointer m_Dwi;

  static PropHelper::GradientDirectionType Vec(double x, double y, double z)
  {
    PropHelper::GradientDirectionType v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
  }

  /** One voxel, one component per gradient. */
  static mitk::Image::Pointer CreateDwi(const std::vector<short>& values, const std::vector<PropHelper::GradientDirectionType>& gradients)
  {
    auto itkImage = PropHelper::ImageType::New();
    PropHelper::ImageType::RegionType region;
    region.SetSize({1, 1, 1});
    itkImage->SetRegions(region);
    itkImage->SetVectorLength(values.size());
    itkImage->Allocate();
    PropHelper::ImageType::PixelType pixel(values.size());
    for (unsigned int i = 0; i < values.size(); ++i)
      pixel[i] = values[i];
    itkImage->SetPixel({0, 0, 0}, pixel);

    auto image = mitk::GrabItkImageMemory(itkImage);
    auto container = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : gradients)
      container->push_back(g);
    PropHelper::SetOriginalGradientContainer(image, container);
    PropHelper::SetReferenceBValue(image, 1000);
    PropHelper::InitializeImage(image);
    return image;
  }

  /** Replaces the rotation part of the image matrix, keeping spacing and origin. */
  static void SetRotation(mitk::Image* image, const vnl_matrix_fixed<double, 3, 3>& rotation)
  {
    auto geometry = image->GetGeometry();
    mitk::AffineTransform3D::MatrixType matrix;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        matrix[r][c] = rotation[r][c] * geometry->GetSpacing()[c];
    auto transform = mitk::AffineTransform3D::New();
    transform->SetMatrix(matrix);
    transform->SetOffset(geometry->GetIndexToWorldTransform()->GetOffset());
    geometry->SetIndexToWorldTransform(transform);
  }

public:
  void setUp() override
  {
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test.dwi"));
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void GetMeasurementFrame_Missing_ReturnsIdentity()
  {
    m_Dwi->GetPropertyList()->RemoveProperty(PropHelper::GetMeasurementFramePropertyName());
    PropHelper::MeasurementFrameType identity;
    identity.set_identity();
    CPPUNIT_ASSERT(PropHelper::GetMeasurementFrame(m_Dwi) == identity);
  }

  void InitializeImage_NoMeasurementFrame_KeepsGradients()
  {
    // Filters set only the modified gradients and no measurement frame before InitializeImage;
    // it then derives the original gradients while reading the (missing) frame.
    auto expected = PropHelper::GetGradientContainer(m_Dwi);
    auto image = m_Dwi->Clone();
    image->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    image->GetPropertyList()->RemoveProperty(PropHelper::GetMeasurementFramePropertyName());
    PropHelper::InitializeImage(image);

    auto original = PropHelper::GetOriginalGradientContainer(image);
    CPPUNIT_ASSERT_EQUAL(expected->Size(), original->Size());
    for (unsigned int i = 0; i < expected->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (expected->ElementAt(i) - original->ElementAt(i)).two_norm(), 1e-12);
  }

  void Getters_NoDiffusionProperties_Throw()
  {
    auto image = mitk::Image::New();
    image->Initialize(m_Dwi);
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT_THROW(PropHelper::GetGradientContainer(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetOriginalGradientContainer(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetReferenceBValue(image), mitk::Exception);
    CPPUNIT_ASSERT_THROW(PropHelper::GetBValueMap(image), mitk::Exception);
  }

  void AverageRedundantGradients_NegativeValues_AveragesSigned()
  {
    auto image = CreateDwi({-10, -20, 5}, {Vec(1, 0, 0), Vec(1, 0, 0), Vec(0, 1, 0)});
    PropHelper::AverageRedundantGradients(image, 0.001);

    CPPUNIT_ASSERT_EQUAL(size_t(2), image->GetPixelType().GetNumberOfComponents());
    CPPUNIT_ASSERT_EQUAL(2u, PropHelper::GetGradientContainer(image)->Size());
    auto itkImage = PropHelper::GetItkVectorImage(image);
    auto pixel = itkImage->GetPixel({0, 0, 0});
    CPPUNIT_ASSERT_EQUAL(short(-15), pixel[0]);
    CPPUNIT_ASSERT_EQUAL(short(5), pixel[1]);
  }

  void InitializeImage_RotationAndFrame_DerivesOriginalGradientsExactly()
  {
    // Rotation and frame that do not commute: deriving original from modified gradients must undo
    // them in reverse order.
    auto image = m_Dwi->Clone();
    SetRotation(image, mitk::imv::GetRotationMatrixVnl<double>(30, 0, 0));
    PropHelper::SetMeasurementFrame(image, mitk::imv::GetRotationMatrixVnl<double>(0, 0, 40));
    PropHelper::SetApplyMatrixToGradients(image, true);
    PropHelper::InitializeImage(image);
    auto original = PropHelper::GetOriginalGradientContainer(image);

    image->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    PropHelper::InitializeImage(image);
    auto derived = PropHelper::GetOriginalGradientContainer(image);
    CPPUNIT_ASSERT_EQUAL(original->Size(), derived->Size());
    for (unsigned int i = 0; i < original->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (original->ElementAt(i) - derived->ElementAt(i)).two_norm(), 1e-12);
  }

  void InitializeImage_Flags_SelectRotationAndFrame()
  {
    // Row vectors: modified = original * R * MF, each factor only if its flag is set.
    const auto rotation = mitk::imv::GetRotationMatrixVnl<double>(30, 0, 0);
    const auto frame = mitk::imv::GetRotationMatrixVnl<double>(0, 0, 40);
    PropHelper::MeasurementFrameType identity;
    identity.set_identity();

    for (bool applyRotation : {false, true})
    {
      for (bool applyFrame : {false, true})
      {
        auto image = m_Dwi->Clone();
        SetRotation(image, rotation);
        PropHelper::SetMeasurementFrame(image, frame);
        PropHelper::SetApplyMatrixToGradients(image, applyRotation);
        PropHelper::SetApplyMfToGradients(image, applyFrame);
        PropHelper::InitializeImage(image);

        const auto original = PropHelper::GetOriginalGradientContainer(image);
        const auto modified = PropHelper::GetGradientContainer(image);
        for (unsigned int i = 0; i < original->Size(); ++i)
        {
          const auto expected = original->ElementAt(i) * (applyRotation ? rotation : identity) * (applyFrame ? frame : identity);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (expected - modified->ElementAt(i)).two_norm(), 1e-12);
        }
      }
    }
  }

  void InitializeImage_Gradients_BValueMapByRoundedBValue()
  {
    auto image = CreateDwi({1, 2, 3, 4, 5}, {Vec(0, 0, 0), Vec(1, 0, 0), Vec(0, std::sqrt(0.4999), 0), Vec(0, 0, std::sqrt(0.5)), Vec(0, -1, 0)});
    const auto& map = PropHelper::GetBValueMap(image);
    CPPUNIT_ASSERT_EQUAL(size_t(3), map.size());
    CPPUNIT_ASSERT(map.at(0) == std::vector<unsigned int>({0}));
    // 1000 * 0.4999 rounds to 500, like 1000 * 0.5.
    CPPUNIT_ASSERT(map.at(500) == std::vector<unsigned int>({2, 3}));
    CPPUNIT_ASSERT(map.at(1000) == std::vector<unsigned int>({1, 4}));
  }

  void IsDiffusionWeightedImage_IncompleteImages_False()
  {
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(m_Dwi));
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(static_cast<const mitk::Image*>(nullptr)));

    auto noBValue = m_Dwi->Clone();
    noBValue->GetPropertyList()->RemoveProperty("DWMRI.ReferenceBValue");
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(noBValue));

    auto tooFewGradients = m_Dwi->Clone();
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    gradients->push_back(Vec(1, 0, 0));
    PropHelper::SetGradientContainer(tooFewGradients, gradients);
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(tooFewGradients));
  }

  void CopyProperties_IgnoreOriginal_RemovesOriginalGradients()
  {
    auto all = mitk::Image::New();
    all->Initialize(m_Dwi);
    PropHelper::CopyProperties(m_Dwi, all);
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(all));
    CPPUNIT_ASSERT(all->GetProperty(PropHelper::GetOriginalGradientContainerPropertyName().c_str()).IsNotNull());

    auto withoutOriginal = mitk::Image::New();
    withoutOriginal->Initialize(m_Dwi);
    PropHelper::CopyProperties(m_Dwi, withoutOriginal, true);
    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(withoutOriginal));
    CPPUNIT_ASSERT(withoutOriginal->GetProperty(PropHelper::GetOriginalGradientContainerPropertyName().c_str()).IsNull());
  }

  void CalcAveragedDirectionSet_AntipodalAndClose_Merged()
  {
    auto directions = PropHelper::GradientDirectionsContainerType::New();
    directions->push_back(Vec(1, 0, 0));
    directions->push_back(Vec(-1, 0, 0));
    directions->push_back(Vec(0.9999, 0.0001, 0));
    directions->push_back(Vec(0, 1, 0));
    auto merged = PropHelper::CalcAveragedDirectionSet(0.001, directions.GetPointer());
    CPPUNIT_ASSERT_EQUAL(2u, merged->Size());
    CPPUNIT_ASSERT(merged->ElementAt(0) == Vec(1, 0, 0));
    CPPUNIT_ASSERT(merged->ElementAt(1) == Vec(0, 1, 0));
  }

  void NodePredicateIsDWI_Nodes_ChecksData()
  {
    auto predicate = mitk::NodePredicateIsDWI::New();
    auto node = mitk::DataNode::New();
    CPPUNIT_ASSERT(!predicate->CheckNode(node));
    node->SetData(m_Dwi);
    CPPUNIT_ASSERT(predicate->CheckNode(node));
    auto plain = mitk::Image::New();
    plain->Initialize(m_Dwi);
    node->SetData(plain);
    CPPUNIT_ASSERT(!predicate->CheckNode(node));
    CPPUNIT_ASSERT_THROW(predicate->CheckNode(nullptr), std::invalid_argument);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionPropertyHelper)
