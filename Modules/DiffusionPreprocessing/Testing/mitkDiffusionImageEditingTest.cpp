/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

#include <algorithm>
#include <limits>

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
  MITK_TEST(ApplyHeader_TimeSeries_ChangesEveryTimeStep);
  MITK_TEST(ClearGradientRotation_RotatedLegacyImage_RemovesRotation);
  MITK_TEST(MergeDwis_RotatedLegacyImages_WorldGradientsWithoutRotation);
  MITK_TEST(MergeDwis_NotADwiOrOtherGeometry_Throws);
  MITK_TEST(Crop_DwiAndVolume_KeepsWorldPositionOfVoxels);
  MITK_TEST(Crop_NothingLeft_Throws);
  MITK_TEST(Resample_DwiAndVolume_SameGridAndValues);
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

  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;

  static VolumeType::Pointer ToItkVolume(const mitk::Image* image)
  {
    VolumeType::Pointer itkImage;
    mitk::CastToItkImage(image, itkImage);
    return itkImage;
  }

  /** Same geometry and the same values as channel c of the DWI. */
  static void AssertVolumeIsChannel(mitk::Image* volume, mitk::Image* dwi, unsigned int c)
  {
    CPPUNIT_ASSERT(mitk::Equal(*dwi->GetGeometry(), *volume->GetGeometry(), 1e-9, 1e-9, true));
    auto itkVolume = ToItkVolume(volume);
    auto itkDwi = PropHelper::GetItkVectorImage(dwi);
    itk::ImageRegionConstIterator<VolumeType> v(itkVolume, itkVolume->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> d(itkDwi, itkDwi->GetLargestPossibleRegion());
    for (; !v.IsAtEnd(); ++v, ++d)
      CPPUNIT_ASSERT_EQUAL(d.Get()[c], v.Get());
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
    // For an axis-aligned image the world component of the flipped axis is negated exactly.
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

  void ApplyHeader_TimeSeries_ChangesEveryTimeStep()
  {
    // A 3D+t image like the "_B0_ALL" result of the preprocessing view, which "Align to" accepts.
    const unsigned int dimensions[] = {3, 2, 2, 3};
    auto image = mitk::Image::New();
    image->Initialize(mitk::MakeScalarPixelType<float>(), 4, dimensions);
    const auto inputGeometry = image->GetGeometry(0)->Clone();

    mitk::DiffusionImageEditing::HeaderChange change;
    mitk::Vector3D spacing;
    spacing[0] = 0.5;
    spacing[1] = 2;
    spacing[2] = 3;
    change.spacing = spacing;
    mitk::Point3D origin;
    origin[0] = 1.5;
    origin[1] = -2;
    origin[2] = 7;
    change.origin = origin;
    change.direction = RotationZ();
    auto result = mitk::DiffusionImageEditing::ApplyHeader(image, change);

    mitk::AffineTransform3D::MatrixType expectedMatrix;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        expectedMatrix[r][c] = RotationZ()[r][c] * spacing[c];
    CPPUNIT_ASSERT_EQUAL(3u, result->GetTimeSteps());
    for (unsigned int t = 0; t < 3; ++t)
    {
      const auto geometry = result->GetGeometry(t);
      CPPUNIT_ASSERT(mitk::Equal(origin, geometry->GetOrigin(), 1e-12, true));
      CPPUNIT_ASSERT(mitk::Equal(spacing, geometry->GetSpacing(), 1e-12, true));
      CPPUNIT_ASSERT(mitk::MatrixEqualElementWise(expectedMatrix, geometry->GetIndexToWorldTransform()->GetMatrix(), 1e-12));
      // The input stays as it was.
      CPPUNIT_ASSERT(mitk::Equal(*inputGeometry, *image->GetGeometry(t), 1e-12, 1e-12, true));
    }

    // The bounding box over all time steps is that of the moved geometry.
    mitk::Point3D lower, upper;
    lower.Fill(std::numeric_limits<double>::max());
    upper.Fill(std::numeric_limits<double>::lowest());
    for (int i = 0; i < 8; ++i)
    {
      const auto corner = result->GetGeometry(0)->GetCornerPoint(i);
      for (unsigned int d = 0; d < 3; ++d)
      {
        lower[d] = std::min(lower[d], corner[d]);
        upper[d] = std::max(upper[d], corner[d]);
      }
    }
    const auto bounds = result->GetTimeGeometry()->GetBoundsInWorld();
    for (unsigned int d = 0; d < 3; ++d)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(lower[d], bounds[2 * d], 1e-9);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(upper[d], bounds[2 * d + 1], 1e-9);
    }
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
  void MergeDwis_RotatedLegacyImages_WorldGradientsWithoutRotation()
  {
    SetRotation(m_Dwi, RotationZ());
    auto image = WithImageRotationFlag(true);
    const auto gradients = PropHelper::GetGradientContainer(image);

    auto merged = mitk::DiffusionImageEditing::MergeDwis({image, image});

    auto expected = PropHelper::GradientDirectionsContainerType::New();
    for (unsigned int copy = 0; copy < 2; ++copy)
      for (const auto& g : gradients->CastToSTLConstContainer())
        expected->push_back(g);
    AssertGradientsEqual(expected.GetPointer(), PropHelper::GetGradientContainer(merged), 1e-9);
    CPPUNIT_ASSERT(!mitk::DiffusionImageEditing::AppliesImageRotationToGradients(merged));
    AssertStableUnderInitialize(merged);
    CPPUNIT_ASSERT(mitk::Equal(*image->GetGeometry(), *merged->GetGeometry(), mitk::eps, true));
  }

  void MergeDwis_NotADwiOrOtherGeometry_Throws()
  {
    auto scalar = mitk::Image::New();
    scalar->Initialize(mitk::MakeScalarPixelType<short>(), *m_Dwi->GetGeometry());
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::MergeDwis({m_Dwi, scalar}), mitk::Exception);

    mitk::DiffusionImageEditing::HeaderChange change;
    mitk::Point3D origin;
    origin.Fill(10);
    change.origin = origin;
    auto moved = mitk::DiffusionImageEditing::ApplyHeader(m_Dwi, change);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::MergeDwis({m_Dwi, moved}), itk::ExceptionObject);
  }

  void Crop_DwiAndVolume_KeepsWorldPositionOfVoxels()
  {
    // test.dwi: 11 x 11 x 3 voxels; turned so that the shift of the origin runs through the
    // direction matrix.
    SetRotation(m_Dwi, RotationZ());
    PropHelper::InitializeImage(m_Dwi);
    const std::array<unsigned int, 3> lower = {1, 2, 0};
    const std::array<unsigned int, 3> upper = {2, 1, 1};
    auto cropped = mitk::DiffusionImageEditing::Crop(m_Dwi, lower, upper);

    CPPUNIT_ASSERT_EQUAL(8u, cropped->GetDimension(0));
    CPPUNIT_ASSERT_EQUAL(8u, cropped->GetDimension(1));
    CPPUNIT_ASSERT_EQUAL(2u, cropped->GetDimension(2));
    AssertGradientsEqual(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(cropped), 1e-12);

    // Voxel (i, j, k) of the result is voxel (i + 1, j + 2, k) of the input, at the same place.
    mitk::Point3D first;
    m_Dwi->GetGeometry()->IndexToWorld(mitk::Point3D(std::array<double, 3>{1, 2, 0}.data()), first);
    for (unsigned int d = 0; d < 3; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(first[d], cropped->GetGeometry()->GetOrigin()[d], 1e-9);
    auto in = PropHelper::GetItkVectorImage(m_Dwi);
    auto out = PropHelper::GetItkVectorImage(cropped);
    itk::ImageRegionConstIteratorWithIndex<DwiType> it(out, out->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      auto index = it.GetIndex();
      index[0] += 1;
      index[1] += 2;
      CPPUNIT_ASSERT(it.Get() == in->GetPixel(index));
    }

    // A volume of the DWI is cropped alike.
    auto volume = mitk::DiffusionImageEditing::Crop(mitk::DiffusionImageEditing::ExtractVolume(m_Dwi, 3), lower, upper);
    AssertVolumeIsChannel(volume, cropped, 3);
  }

  void Crop_NothingLeft_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::Crop(m_Dwi, {5, 0, 0}, {6, 0, 0}), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::Crop(m_Dwi, {0, 0, 3}, {0, 0, 0}), mitk::Exception);
  }

  void Resample_DwiAndVolume_SameGridAndValues()
  {
    // A flipped y axis, as ITK reads most clinical NIfTI files.
    MatrixType flip;
    flip.set_identity();
    flip[1][1] = -1;
    SetRotation(m_Dwi, flip);
    PropHelper::InitializeImage(m_Dwi);
    auto volume = mitk::DiffusionImageEditing::ExtractVolume(m_Dwi, 7);

    using mitk::DiffusionImageEditing::Interpolation;
    using mitk::DiffusionImageEditing::ResampleGrid;
    const std::vector<std::pair<ResampleGrid, mitk::Vector3D>> grids = {
      {ResampleGrid::SamplingFactor, mitk::Vector3D(std::array<double, 3>{1.5, 0.7, 2.0}.data())},
      {ResampleGrid::Spacing, mitk::Vector3D(std::array<double, 3>{3.0, 1.0, 2.5}.data())},
      {ResampleGrid::Size, mitk::Vector3D(std::array<double, 3>{5, 20, 3}.data())}};
    for (const auto& [grid, value] : grids)
    {
      for (auto interpolation : {Interpolation::NearestNeighbour, Interpolation::Linear, Interpolation::BSpline, Interpolation::WindowedSinc})
      {
        auto resampledDwi = mitk::DiffusionImageEditing::Resample(m_Dwi, grid, value, interpolation);
        AssertGradientsEqual(PropHelper::GetGradientContainer(m_Dwi), PropHelper::GetGradientContainer(resampledDwi), 1e-12);
        AssertVolumeIsChannel(mitk::DiffusionImageEditing::Resample(volume, grid, value, interpolation), resampledDwi, 7);
      }
    }
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::Resample(m_Dwi, ResampleGrid::Spacing, mitk::Vector3D(0.0), Interpolation::Linear), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageEditing)
