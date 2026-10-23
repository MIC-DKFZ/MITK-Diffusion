/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkAffineTransform.h>
#include <itkImageRegionConstIterator.h>
#include <mapPreCachedRegistrationKernel.h>
#include <mapRegistration.h>
#include <mapRegistrationManipulator.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionImageMapping.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkMAPAlgorithmHelper.h>

class mitkDiffusionImageMappingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageMappingTestSuite);
  MITK_TEST(RotateGradients_OneAffine_RotatesAllWeightedGradients);
  MITK_TEST(RotateGradients_OneMatrixNoB0_RotatesFirstVolumeToo);
  MITK_TEST(RotateGradients_OneAffinePerVolume_RotatesEachVolumeByItsOwn);
  MITK_TEST(RotationComponent_SingularMatrix_Throws);
  MITK_TEST(MapImage_IdentityRegistration_KeepsVolumesAndGradients);
  MITK_TEST(MapImage_Rotation_MapsEachVolumeAndRotatesGradients);
  MITK_TEST(RefineGeometry_Rotation_TurnsGeometryAndGradients);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;
  using MatrixType = mitk::DiffusionImageMapping::MatrixType;

  mitk::Image::Pointer m_Image;
  MatrixType m_Rotation;

  /** A registration whose direct kernel turns about the image centre, the inverse kernel back. */
  mitk::MAPRegistrationWrapper::Pointer RotationRegistration() const
  {
    using RegistrationType = map::core::Registration<3, 3>;
    using KernelType = map::core::PreCachedRegistrationKernel<3, 3>;
    using TransformType = itk::AffineTransform<double, 3>;

    auto direct = TransformType::New();
    direct->SetCenter(m_Image->GetGeometry()->GetCenter());
    TransformType::MatrixType matrix;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        matrix[r][c] = m_Rotation[r][c];
    direct->SetMatrix(matrix);
    auto inverse = TransformType::New();
    direct->GetInverse(inverse);

    auto directKernel = KernelType::New();
    directKernel->setTransformModel(direct);
    auto inverseKernel = KernelType::New();
    inverseKernel->setTransformModel(inverse);

    auto registration = RegistrationType::New();
    map::core::RegistrationManipulator<RegistrationType> manipulator(registration);
    manipulator.setDirectMapping(directKernel);
    manipulator.setInverseMapping(inverseKernel);
    return mitk::MAPRegistrationWrapper::New(registration);
  }

  static void CheckGradients(const PropHelper::GradientDirectionsContainerType* expected, const mitk::Image* image)
  {
    auto output = PropHelper::GetGradientContainer(image);
    CPPUNIT_ASSERT_EQUAL(expected->Size(), output->Size());
    for (unsigned int i = 0; i < expected->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (output->ElementAt(i) - expected->ElementAt(i)).two_norm(), 1e-9);
  }

  /** The gradients of m_Image, each multiplied by its matrix. */
  PropHelper::GradientDirectionsContainerType::Pointer Rotated(const std::vector<MatrixType>& rotations) const
  {
    auto input = PropHelper::GetGradientContainer(m_Image);
    auto result = PropHelper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < input->Size(); ++i)
      result->push_back(rotations[i] * input->ElementAt(i));
    return result;
  }

  static VolumeType::Pointer ToItk(const mitk::Image* image)
  {
    VolumeType::Pointer itkImage;
    mitk::CastToItkImage(image, itkImage);
    return itkImage;
  }

  static void CheckEqualVolumes(const VolumeType* expected, const DwiType* dwi, unsigned int channel)
  {
    itk::ImageRegionConstIterator<VolumeType> e(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> d(dwi, dwi->GetLargestPossibleRegion());
    for (; !e.IsAtEnd(); ++e, ++d)
      CPPUNIT_ASSERT_EQUAL(e.Get(), d.Get()[channel]);
  }

public:
  void setUp() override
  {
    m_Image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"));
    // 90 degrees about z
    m_Rotation.fill(0.0);
    m_Rotation[0][1] = -1;
    m_Rotation[1][0] = 1;
    m_Rotation[2][2] = 1;
  }

  void tearDown() override
  {
    m_Image = nullptr;
  }

  void RotateGradients_OneAffine_RotatesAllWeightedGradients()
  {
    // The rotation part of rotation * scaling is the rotation.
    MatrixType scaling;
    scaling.set_identity();
    scaling[0][0] = 2;
    auto expected = Rotated(std::vector<MatrixType>(PropHelper::GetGradientContainer(m_Image)->Size(), m_Rotation));
    mitk::DiffusionImageMapping::RotateGradients(m_Image, m_Rotation * scaling);
    CheckGradients(expected, m_Image);
  }

  void RotateGradients_OneMatrixNoB0_RotatesFirstVolumeToo()
  {
    // Volume 0, the only b0 volume of test.dwi, becomes a weighted one: the lowest shell is then
    // b = 3000 and starts with volume 0.
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : PropHelper::GetGradientContainer(m_Image)->CastToSTLConstContainer())
      gradients->push_back(g);
    gradients->ElementAt(0).fill(0.0);
    gradients->ElementAt(0)[0] = 1.0;
    m_Image->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    PropHelper::SetGradientContainer(m_Image, gradients);
    PropHelper::InitializeImage(m_Image);
    CPPUNIT_ASSERT_EQUAL(size_t(1), PropHelper::GetBValueMap(m_Image).size());

    auto expected = Rotated(std::vector<MatrixType>(gradients->Size(), m_Rotation));
    mitk::DiffusionImageMapping::RotateGradients(m_Image, m_Rotation);
    CheckGradients(expected, m_Image);
  }

  void RotateGradients_OneAffinePerVolume_RotatesEachVolumeByItsOwn()
  {
    const unsigned int n = PropHelper::GetGradientContainer(m_Image)->Size();
    MatrixType identity;
    identity.set_identity();
    std::vector<MatrixType> rotations(n, identity);
    rotations[1] = m_Rotation;
    rotations[2] = m_Rotation.transpose();
    auto expected = Rotated(rotations);
    mitk::DiffusionImageMapping::RotateGradients(m_Image, rotations);
    CheckGradients(expected, m_Image);

    rotations.pop_back();
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageMapping::RotateGradients(m_Image, rotations), mitk::Exception);
  }

  void RotationComponent_SingularMatrix_Throws()
  {
    MatrixType singular = m_Rotation;
    singular[2][2] = 0;
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageMapping::RotationComponent(singular), mitk::Exception);
  }

  void MapImage_IdentityRegistration_KeepsVolumesAndGradients()
  {
    auto result = mitk::DiffusionImageMapping::MapImage(m_Image, mitk::GenerateIdentityRegistration3D(), m_Image->GetGeometry(), mitk::ImageMappingInterpolator::NearestNeighbor);

    CPPUNIT_ASSERT(PropHelper::IsDiffusionWeightedImage(result));
    CheckGradients(PropHelper::GetGradientContainer(m_Image), result);
    CPPUNIT_ASSERT(PropHelper::GetBValueMap(m_Image) == PropHelper::GetBValueMap(result));
    CPPUNIT_ASSERT_EQUAL(PropHelper::GetReferenceBValue(m_Image), PropHelper::GetReferenceBValue(result));

    auto input = PropHelper::GetItkVectorImage(m_Image);
    auto output = PropHelper::GetItkVectorImage(result);
    itk::ImageRegionConstIterator<DwiType> in(input, input->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> out(output, output->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
      CPPUNIT_ASSERT(in.Get() == out.Get());
  }

  void MapImage_Rotation_MapsEachVolumeAndRotatesGradients()
  {
    auto registration = RotationRegistration();
    auto result = mitk::DiffusionImageMapping::MapImage(m_Image, registration, m_Image->GetGeometry(), mitk::ImageMappingInterpolator::Linear);

    const unsigned int n = PropHelper::GetGradientContainer(m_Image)->Size();
    CheckGradients(Rotated(std::vector<MatrixType>(n, m_Rotation)), result);

    // Each volume is mapped as a scalar image would be.
    auto output = PropHelper::GetItkVectorImage(result);
    CPPUNIT_ASSERT_EQUAL(n, output->GetVectorLength());
    for (unsigned int i : {0u, 5u, n - 1})
    {
      auto expected = mitk::ImageMappingHelper::map(mitk::DiffusionImageEditing::ExtractVolume(m_Image, i), registration.GetPointer(), false, 0, m_Image->GetGeometry(), false, 0, mitk::ImageMappingInterpolator::Linear);
      CheckEqualVolumes(ToItk(expected), output, i);
    }
  }

  void RefineGeometry_Rotation_TurnsGeometryAndGradients()
  {
    auto result = mitk::DiffusionImageMapping::RefineGeometry(m_Image, RotationRegistration());

    const unsigned int n = PropHelper::GetGradientContainer(m_Image)->Size();
    CheckGradients(Rotated(std::vector<MatrixType>(n, m_Rotation)), result);

    // The voxels stay as they are; the image axes turn with the registration.
    for (unsigned int c = 0; c < 3; ++c)
    {
      const auto column = result->GetGeometry()->GetMatrixColumn(c);
      const auto expected = m_Rotation * m_Image->GetGeometry()->GetMatrixColumn(c);
      for (unsigned int r = 0; r < 3; ++r)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[r], column[r], 1e-9);
    }
    auto input = PropHelper::GetItkVectorImage(m_Image);
    auto output = PropHelper::GetItkVectorImage(result);
    itk::ImageRegionConstIterator<DwiType> in(input, input->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<DwiType> out(output, output->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
      CPPUNIT_ASSERT(in.Get() == out.Get());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageMapping)
