/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

class mitkDiffusionImageEditingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageEditingTestSuite);
  MITK_TEST(FlipAxes_AxisAligned_MirrorsDataAndNegatesComponent);
  MITK_TEST(FlipAxes_RotatedImage_MirrorsGradientsAlongImageAxis);
  MITK_TEST(FlipAxes_ScalarImage_MirrorsData);
  MITK_TEST(FlipAxes_MultiComponentImage_Throws);
  MITK_TEST(ApplyHeader_NewDirection_GradientsTurnWithImageOrStay);
  MITK_TEST(ApplyHeader_FromRotatedLegacyImage_GradientsTurnBack);
  MITK_TEST(ApplyHeader_NewFrame_AppliedToOriginalGradients);
  MITK_TEST(ApplyHeader_OnlyOrigin_KeepsEverythingElse);
  MITK_TEST(ApplyHeader_InvalidValues_Throw);
  MITK_TEST(ClearGradientRotation_RotatedLegacyImage_RemovesRotation);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using MatrixType = mitk::DiffusionImageEditing::MatrixType;

  mitk::Image::Pointer m_Dwi;

  /** Replaces the rotation part of the image matrix, keeping spacing and origin. */
  static void SetRotation(mitk::Image* image, const MatrixType& rotation)
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

  /** 90 degrees about z: image axis x points along world y, image axis y along world -x. */
  static MatrixType RotationZ()
  {
    MatrixType rotation;
    rotation.fill(0.0);
    rotation[1][0] = 1;
    rotation[0][1] = -1;
    rotation[2][2] = 1;
    return rotation;
  }

  static void AssertGradientsEqual(PropHelper::GradientDirectionsContainerType::ConstPointer expected, PropHelper::GradientDirectionsContainerType::ConstPointer actual, double tolerance)
  {
    CPPUNIT_ASSERT_EQUAL(expected->Size(), actual->Size());
    for (unsigned int i = 0; i < expected->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (expected->ElementAt(i) - actual->ElementAt(i)).two_norm(), tolerance);
  }

  /** Voxel (i, j, k) of the output equals voxel (n - 1 - i, j, k) of the input. */
  static void AssertMirroredAlongX(mitk::Image* input, mitk::Image* output)
  {
    auto in = PropHelper::GetItkVectorImage(input);
    auto out = PropHelper::GetItkVectorImage(output);
    const auto n = static_cast<itk::IndexValueType>(in->GetLargestPossibleRegion().GetSize(0));
    itk::ImageRegionConstIteratorWithIndex<DwiType> it(out, out->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      auto index = it.GetIndex();
      index[0] = n - 1 - index[0];
      CPPUNIT_ASSERT(it.Get() == in->GetPixel(index));
    }
  }

  /** The DWI with the image-rotation flag set as given and its gradients derived again. */
  mitk::Image::Pointer WithImageRotationFlag(bool applyImageRotation)
  {
    auto image = m_Dwi->Clone();
    PropHelper::SetApplyMatrixToGradients(image, applyImageRotation);
    PropHelper::InitializeImage(image);
    return image;
  }

  static PropHelper::GradientDirectionsContainerType::Pointer Transformed(PropHelper::GradientDirectionsContainerType::ConstPointer gradients, const MatrixType& matrix)
  {
    auto result = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : gradients->CastToSTLConstContainer())
      result->push_back(matrix * g);
    return result;
  }

  /** The gradients stay as they are when InitializeImage derives them again, as later edits do. */
  static void AssertStableUnderInitialize(mitk::Image* image)
  {
    auto expected = PropHelper::GetGradientContainer(image);
    auto copy = image->Clone();
    PropHelper::InitializeImage(copy);
    AssertGradientsEqual(expected, PropHelper::GetGradientContainer(copy), 1e-12);
  }

public:
  void setUp() override
  {
    m_Dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
  }

  void tearDown() override
  {
    m_Dwi = nullptr;
  }

  void FlipAxes_AxisAligned_MirrorsDataAndNegatesComponent()
  {
    auto flipped = mitk::DiffusionImageEditing::FlipAxes(m_Dwi, {true, false, false});

    AssertMirroredAlongX(m_Dwi, flipped);
    CPPUNIT_ASSERT(mitk::Equal(*m_Dwi->GetGeometry(), *flipped->GetGeometry(), mitk::eps, true));
    auto expected = PropHelper::GradientDirectionsContainerType::New();
    for (auto g : PropHelper::GetGradientContainer(m_Dwi)->CastToSTLConstContainer())
    {
      g[0] = -g[0];
      expected->push_back(g);
    }
    // Exactly as before the gradient convention fix for axis-aligned images.
    AssertGradientsEqual(expected.GetPointer(), PropHelper::GetGradientContainer(flipped), 0.0);
  }

  void FlipAxes_RotatedImage_MirrorsGradientsAlongImageAxis()
  {
    SetRotation(m_Dwi, RotationZ());
    auto flipped = mitk::DiffusionImageEditing::FlipAxes(m_Dwi, {true, false, false});

    // Image axis x points along world y: the world y component is negated.
    AssertMirroredAlongX(m_Dwi, flipped);
    CPPUNIT_ASSERT(mitk::Equal(*m_Dwi->GetGeometry(), *flipped->GetGeometry(), mitk::eps, true));
    auto expected = PropHelper::GradientDirectionsContainerType::New();
    for (auto g : PropHelper::GetGradientContainer(m_Dwi)->CastToSTLConstContainer())
    {
      g[1] = -g[1];
      expected->push_back(g);
    }
    AssertGradientsEqual(expected.GetPointer(), PropHelper::GetGradientContainer(flipped), 1e-12);

    auto flippedBack = mitk::DiffusionImageEditing::FlipAxes(flipped, {true, false, false});
    AssertGradientsEqual(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(flippedBack), 1e-12);
  }

  void FlipAxes_ScalarImage_MirrorsData()
  {
    using ScalarType = itk::Image<float, 3>;
    auto itkImage = ScalarType::New();
    ScalarType::RegionType region;
    region.SetSize({3, 2, 1});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    for (unsigned int x = 0; x < 3; ++x)
      for (unsigned int y = 0; y < 2; ++y)
        itkImage->SetPixel({x, y, 0}, 10.0f * x + y);
    ScalarType::PointType origin;
    origin[0] = 5;
    origin[1] = -3;
    origin[2] = 2;
    itkImage->SetOrigin(origin);
    auto image = mitk::GrabItkImageMemory(itkImage.GetPointer());

    auto flipped = mitk::DiffusionImageEditing::FlipAxes(image, {false, true, false});

    CPPUNIT_ASSERT(mitk::Equal(*image->GetGeometry(), *flipped->GetGeometry(), mitk::eps, true));
    ScalarType::Pointer itkFlipped;
    mitk::CastToItkImage(flipped, itkFlipped);
    for (unsigned int x = 0; x < 3; ++x)
      for (unsigned int y = 0; y < 2; ++y)
        CPPUNIT_ASSERT_EQUAL(10.0f * x + (1 - y), itkFlipped->GetPixel({x, y, 0}));
  }

  void FlipAxes_MultiComponentImage_Throws()
  {
    // A vector image without diffusion properties.
    auto itkImage = DwiType::New();
    DwiType::RegionType region;
    region.SetSize({2, 2, 2});
    itkImage->SetRegions(region);
    itkImage->SetVectorLength(2);
    itkImage->Allocate();
    auto image = mitk::GrabItkImageMemory(itkImage.GetPointer());
    CPPUNIT_ASSERT(!PropHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::FlipAxes(image, {true, false, false}), mitk::Exception);
  }
  void ApplyHeader_NewDirection_GradientsTurnWithImageOrStay()
  {
    // The result must not depend on the image-rotation flag of the input.
    for (bool flag : {true, false})
    {
      auto image = WithImageRotationFlag(flag);
      const auto gradients = PropHelper::GetGradientContainer(image);
      for (bool rotate : {true, false})
      {
        mitk::DiffusionImageEditing::HeaderChange change;
        change.direction = RotationZ();
        change.rotateGradientsWithImage = rotate;
        auto result = mitk::DiffusionImageEditing::ApplyHeader(image, change);

        // test.dwi has an identity direction matrix, so turning with the image applies RotationZ.
        MatrixType identity;
        identity.set_identity();
        AssertGradientsEqual(Transformed(gradients, rotate ? RotationZ() : identity).GetPointer(), PropHelper::GetGradientContainer(result), 1e-12);
        CPPUNIT_ASSERT(!mitk::DiffusionImageEditing::AppliesImageRotationToGradients(result));
        AssertStableUnderInitialize(result);
        CPPUNIT_ASSERT(PropHelper::GetOriginalGradientContainer(result)->Size() == gradients->Size());
      }
    }
  }

  void ApplyHeader_FromRotatedLegacyImage_GradientsTurnBack()
  {
    // A legacy image applies its rotation to the gradients; resetting the direction matrix with
    // "rotate with the image" turns them back by the inverse rotation.
    SetRotation(m_Dwi, RotationZ());
    auto image = WithImageRotationFlag(true);
    const auto gradients = PropHelper::GetGradientContainer(image);

    mitk::DiffusionImageEditing::HeaderChange change;
    MatrixType identity;
    identity.set_identity();
    change.direction = identity;
    auto result = mitk::DiffusionImageEditing::ApplyHeader(image, change);

    AssertGradientsEqual(Transformed(gradients, RotationZ().transpose()).GetPointer(), PropHelper::GetGradientContainer(result), 1e-12);
    AssertStableUnderInitialize(result);
  }

  void ApplyHeader_NewFrame_AppliedToOriginalGradients()
  {
    auto image = WithImageRotationFlag(false);
    const auto original = PropHelper::GetOriginalGradientContainer(image);

    // Frame axes as rows: modified = frame^T * original.
    mitk::DiffusionImageEditing::HeaderChange change;
    change.measurementFrame = RotationZ();
    auto result = mitk::DiffusionImageEditing::ApplyHeader(image, change);

    AssertGradientsEqual(Transformed(original, RotationZ().transpose()).GetPointer(), PropHelper::GetGradientContainer(result), 1e-12);
    AssertGradientsEqual(original, PropHelper::GetOriginalGradientContainer(result), 0.0);
    AssertStableUnderInitialize(result);
  }

  void ApplyHeader_OnlyOrigin_KeepsEverythingElse()
  {
    mitk::DiffusionImageEditing::HeaderChange change;
    mitk::Point3D origin;
    origin[0] = 1.5;
    origin[1] = -2;
    origin[2] = 7;
    change.origin = origin;
    auto result = mitk::DiffusionImageEditing::ApplyHeader(m_Dwi, change);

    CPPUNIT_ASSERT(mitk::Equal(origin, result->GetGeometry()->GetOrigin(), 1e-12, true));
    CPPUNIT_ASSERT(mitk::Equal(m_Dwi->GetGeometry()->GetSpacing(), result->GetGeometry()->GetSpacing(), 1e-12, true));
    CPPUNIT_ASSERT(mitk::MatrixEqualElementWise(m_Dwi->GetGeometry()->GetIndexToWorldTransform()->GetMatrix(), result->GetGeometry()->GetIndexToWorldTransform()->GetMatrix(), 1e-12));
    AssertGradientsEqual(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(result), 0.0);
  }

  void ApplyHeader_InvalidValues_Throw()
  {
    mitk::DiffusionImageEditing::HeaderChange zeroSpacing;
    mitk::Vector3D spacing;
    spacing.Fill(1.0);
    spacing[1] = 0;
    zeroSpacing.spacing = spacing;
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::ApplyHeader(m_Dwi, zeroSpacing), mitk::Exception);

    mitk::DiffusionImageEditing::HeaderChange skewed;
    MatrixType direction;
    direction.set_identity();
    direction[0][1] = 0.5;
    skewed.direction = direction;
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::ApplyHeader(m_Dwi, skewed), mitk::Exception);
  }
  void ClearGradientRotation_RotatedLegacyImage_RemovesRotation()
  {
    SetRotation(m_Dwi, RotationZ());
    auto image = WithImageRotationFlag(true);
    CPPUNIT_ASSERT(mitk::DiffusionImageEditing::AppliesImageRotationToGradients(image));
    const auto original = PropHelper::GetOriginalGradientContainer(image);

    auto cleared = mitk::DiffusionImageEditing::ClearGradientRotation(image);

    // The modified gradients change at once: the frame applied to the originals, no rotation.
    CPPUNIT_ASSERT(!mitk::DiffusionImageEditing::AppliesImageRotationToGradients(cleared));
    AssertGradientsEqual(Transformed(original, PropHelper::GetMeasurementFrame(image).transpose()).GetPointer(), PropHelper::GetGradientContainer(cleared), 1e-12);
    AssertStableUnderInitialize(cleared);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageEditing)
