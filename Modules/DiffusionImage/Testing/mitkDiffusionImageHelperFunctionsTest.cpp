/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkIOUtil.h>

#include <itksys/SystemTools.hxx>

#include <algorithm>
#include <format>
#include <fstream>

class mitkDiffusionImageHelperFunctionsTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionImageHelperFunctionsTestSuite);
  MITK_TEST(ReadBvalsBvecs_ThreeRows_ScalesByBValue);
  MITK_TEST(ReadBvalsBvecs_ManyDigits_ReadsDoublePrecision);
  MITK_TEST(ReadBvalsBvecs_WrongShape_Throws);
  MITK_TEST(ReadBvalsBvecs_InvalidEntry_Throws);
  MITK_TEST(WriteBvalsBvecs_ReadBack_SameGradients);
  MITK_TEST(UniqueDirections_ShellsWithSameDirections_CountedOnce);
  MITK_TEST(FslBvecsMeasurementFrame_Determinant_FlipsXOnlyIfPositive);
  MITK_TEST(IntersectImage_Segment_LengthPerVoxel);
  MITK_TEST(RotationMatrices_VnlAndItk_Equal);
  MITK_TEST(Cart2Sph_Sph2Cart_RoundTrip);
  MITK_TEST(VectorImageToVolumes_VolumesToVectorImage_TransposeAndKeepGeometry);
  CPPUNIT_TEST_SUITE_END();

  using GradientType = mitk::DiffusionPropertyHelper::GradientDirectionType;
  using ContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

  std::string m_TempDir;

  static GradientType Vec(double x, double y, double z)
  {
    GradientType v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
  }

  void WriteFiles(const std::string& bvals, const std::string& bvecs)
  {
    std::ofstream(m_TempDir + "/test.bvals") << bvals;
    std::ofstream(m_TempDir + "/test.bvecs") << bvecs;
  }

  ContainerType::ConstPointer Read(double& bref)
  {
    return mitk::gradients::ReadBvalsBvecs(m_TempDir + "/test.bvals", m_TempDir + "/test.bvecs", bref);
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("DiffusionImageHelperFunctionsTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void ReadBvalsBvecs_ThreeRows_ScalesByBValue()
  {
    WriteFiles("0 1000 2000\n", "0 1 0\n0 0 0.6\n0 0 0.8\n");
    double bref = 0;
    auto gradients = Read(bref);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2000.0, bref, 0);
    CPPUNIT_ASSERT_EQUAL(3u, gradients->Size());
    CPPUNIT_ASSERT(gradients->ElementAt(0) == Vec(0, 0, 0));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (gradients->ElementAt(1) - Vec(std::sqrt(0.5), 0, 0)).two_norm(), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (gradients->ElementAt(2) - Vec(0, 0.6, 0.8)).two_norm(), 1e-15);
  }

  void ReadBvalsBvecs_ManyDigits_ReadsDoublePrecision()
  {
    const double x = 0.12345678901234;
    const double y = std::sqrt(1 - x * x);
    WriteFiles("1000\n", std::format("{}\n{}\n0\n", x, y));
    double bref = 0;
    auto gradients = Read(bref);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(x, gradients->ElementAt(0)[0], 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(y, gradients->ElementAt(0)[1], 1e-15);
  }

  void ReadBvalsBvecs_WrongShape_Throws()
  {
    double bref = 0;
    // One gradient per line instead of three rows.
    WriteFiles("0 1000 1000 1000\n", "0 0 0\n1 0 0\n0 1 0\n0 0 1\n");
    CPPUNIT_ASSERT_THROW(Read(bref), mitk::Exception);
    // One b-value more than gradients.
    WriteFiles("0 1000 1000\n", "0 1\n0 0\n0 0\n");
    CPPUNIT_ASSERT_THROW(Read(bref), mitk::Exception);
    // Missing file.
    CPPUNIT_ASSERT_THROW(mitk::gradients::ReadBvalsBvecs(m_TempDir + "/test.bvals", m_TempDir + "/missing.bvecs", bref), mitk::Exception);
  }

  void ReadBvalsBvecs_InvalidEntry_Throws()
  {
    double bref = 0;
    WriteFiles("0 1000\n", "0 1\n0 zero\n0 0\n");
    CPPUNIT_ASSERT_THROW(Read(bref), mitk::Exception);
  }

  void WriteBvalsBvecs_ReadBack_SameGradients()
  {
    auto gradients = ContainerType::New();
    gradients->push_back(Vec(0, 0, 0));
    gradients->push_back(Vec(0.36, 0.48, 0.8));
    gradients->push_back(Vec(0.6, -0.8, 0) * std::sqrt(0.5));
    mitk::gradients::WriteBvalsBvecs(m_TempDir + "/test.bvals", m_TempDir + "/test.bvecs", gradients.GetPointer(), 2000);

    double bref = 0;
    auto read = Read(bref);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2000.0, bref, 0);
    CPPUNIT_ASSERT_EQUAL(gradients->Size(), read->Size());
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (gradients->ElementAt(i) - read->ElementAt(i)).two_norm(), 1e-15);
  }

  void UniqueDirections_ShellsWithSameDirections_CountedOnce()
  {
    // b0, then the same three directions at b=1000 and at the reference b=2000.
    auto gradients = ContainerType::New();
    gradients->push_back(Vec(0, 0, 0));
    for (double scale : {std::sqrt(0.5), 1.0})
    {
      gradients->push_back(Vec(1, 0, 0) * scale);
      gradients->push_back(Vec(0, 1, 0) * scale);
      gradients->push_back(Vec(0, 0, 1) * scale);
    }
    mitk::BValueMapProperty::BValueMap shells;
    shells[0] = {0};
    shells[1000] = {1, 2, 3};
    shells[2000] = {4, 5, 6};

    CPPUNIT_ASSERT_EQUAL(size_t(3), mitk::gradients::GetAllUniqueDirections(shells, gradients.GetPointer()).size());
    // Plus the b0 entry.
    CPPUNIT_ASSERT_EQUAL(4u, mitk::gradients::CreateNormalizedUniqueGradientDirectionContainer(shells, gradients)->Size());
  }

  void FslBvecsMeasurementFrame_Determinant_FlipsXOnlyIfPositive()
  {
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    auto frame = mitk::gradients::FslBvecsMeasurementFrame(identity);
    mitk::DiffusionPropertyHelper::MeasurementFrameType expected;
    expected.set_identity();
    expected(0, 0) = -1;
    CPPUNIT_ASSERT(frame == expected);

    // LAS storage, as dcm2niix usually writes: diag(1,-1,1) in LPS, determinant -1.
    itk::Matrix<double, 3, 3> las;
    las.SetIdentity();
    las(1, 1) = -1;
    frame = mitk::gradients::FslBvecsMeasurementFrame(las);
    CPPUNIT_ASSERT(frame == las.GetVnlMatrix().transpose());
  }

  void IntersectImage_Segment_LengthPerVoxel()
  {
    itk::Vector<double, 3> spacing;
    spacing[0] = 2;
    spacing[1] = 1;
    spacing[2] = 1;

    // Along x from 0.2 to 1.6 (in index units): 0.3 in voxel 0, 1.0 in voxel 1, 0.1 in voxel 2.
    itk::ContinuousIndex<float, 3> sf, ef;
    sf[0] = 0.2f; sf[1] = 0; sf[2] = 0;
    ef[0] = 1.6f; ef[1] = 0; ef[2] = 0;
    itk::Index<3> si = {{0, 0, 0}}, ei = {{2, 0, 0}};
    auto segments = mitk::imv::IntersectImage(spacing, si, ei, sf, ef);
    CPPUNIT_ASSERT_EQUAL(size_t(3), segments.size());
    const double expected[3] = {0.6, 2.0, 0.2};
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(itk::IndexValueType(i), segments[i].first[0]);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], segments[i].second, 1e-6);
    }

    // Backwards and within one voxel.
    itk::Index<3> backStart = {{2, 0, 0}}, backEnd = {{0, 0, 0}};
    double total = 0;
    for (const auto& segment : mitk::imv::IntersectImage(spacing, backStart, backEnd, ef, sf))
      total += segment.second;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.8, total, 1e-6);

    itk::ContinuousIndex<float, 3> inside = sf;
    inside[1] = 0.3f;
    itk::Index<3> same = {{0, 0, 0}}, same2 = {{0, 0, 0}};
    auto single = mitk::imv::IntersectImage(spacing, same, same2, sf, inside);
    CPPUNIT_ASSERT_EQUAL(size_t(1), single.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, single[0].second, 1e-6);
  }

  void RotationMatrices_VnlAndItk_Equal()
  {
    const auto vnlRotation = mitk::imv::GetRotationMatrixVnl<double>(10, -20, 30);
    const auto itkRotation = mitk::imv::GetRotationMatrixItk<double>(10, -20, 30);
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(vnlRotation(r, c), itkRotation(r, c), 1e-12);
    // z * y * x, so a pure z rotation turns x into y.
    const auto rz = mitk::imv::GetRotationMatrixVnl<double>(0, 0, 90);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (rz * Vec(1, 0, 0) - Vec(0, 1, 0)).two_norm(), 1e-12);
  }

  void Cart2Sph_Sph2Cart_RoundTrip()
  {
    const auto v = Vec(0.3, -0.4, 1.2);
    double spherical[3];
    mitk::gradients::Cart2Sph(v[0], v[1], v[2], spherical);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.two_norm(), spherical[2], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (mitk::gradients::Sph2Cart(spherical[1], spherical[0], spherical[2]) - v).two_norm(), 1e-12);
  }

  void VectorImageToVolumes_VolumesToVectorImage_TransposeAndKeepGeometry()
  {
    // Component t of voxel (x, y, z) is volume t of the 4-D image; an odd size and an oblique
    // geometry catch swapped strides and lost geometry.
    using VectorImageType = itk::VectorImage<short, 3>;
    const unsigned int numberOfVolumes = 5;
    const itk::Size<3> size = {{3, 4, 2}};
    auto value = [](itk::IndexValueType x, itk::IndexValueType y, itk::IndexValueType z, itk::IndexValueType t)
    { return static_cast<short>(x + 10 * y + 100 * z + 1000 * t); };

    auto image = VectorImageType::New();
    image->SetRegions(VectorImageType::RegionType(size));
    image->SetVectorLength(numberOfVolumes);
    image->SetSpacing(itk::MakeVector(0.5, 1.5, 2.0));
    image->SetOrigin(itk::MakePoint(-10.0, 20.0, 5.0));
    image->SetDirection(VectorImageType::DirectionType(mitk::imv::GetRotationMatrixVnl<double>(30, 0, 40)));
    image->Allocate();
    for (itk::IndexValueType z = 0; z < 2; ++z)
      for (itk::IndexValueType y = 0; y < 4; ++y)
        for (itk::IndexValueType x = 0; x < 3; ++x)
        {
          VectorImageType::PixelType pixel(numberOfVolumes);
          for (unsigned int t = 0; t < numberOfVolumes; ++t)
            pixel[t] = value(x, y, z, t);
          image->SetPixel({{x, y, z}}, pixel);
        }

    auto volumes = mitk::imv::VectorImageToVolumes(image.GetPointer());
    const auto size4 = volumes->GetLargestPossibleRegion().GetSize();
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(numberOfVolumes), size4[3]);
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(size[i], size4[i]);
      CPPUNIT_ASSERT_EQUAL(image->GetSpacing()[i], volumes->GetSpacing()[i]);
      CPPUNIT_ASSERT_EQUAL(image->GetOrigin()[i], volumes->GetOrigin()[i]);
      for (unsigned int j = 0; j < 3; ++j)
        CPPUNIT_ASSERT_EQUAL(image->GetDirection()[i][j], volumes->GetDirection()[i][j]);
      CPPUNIT_ASSERT_EQUAL(0.0, volumes->GetDirection()[i][3]);
      CPPUNIT_ASSERT_EQUAL(0.0, volumes->GetDirection()[3][i]);
    }
    CPPUNIT_ASSERT_EQUAL(1.0, volumes->GetSpacing()[3]);
    CPPUNIT_ASSERT_EQUAL(0.0, volumes->GetOrigin()[3]);
    CPPUNIT_ASSERT_EQUAL(1.0, volumes->GetDirection()[3][3]);
    for (itk::IndexValueType t = 0; t < numberOfVolumes; ++t)
      for (itk::IndexValueType z = 0; z < 2; ++z)
        for (itk::IndexValueType y = 0; y < 4; ++y)
          for (itk::IndexValueType x = 0; x < 3; ++x)
            CPPUNIT_ASSERT_EQUAL(value(x, y, z, t), volumes->GetPixel({{x, y, z, t}}));

    auto back = mitk::imv::VolumesToVectorImage(volumes.GetPointer());
    CPPUNIT_ASSERT_EQUAL(numberOfVolumes, back->GetVectorLength());
    CPPUNIT_ASSERT(back->GetLargestPossibleRegion() == image->GetLargestPossibleRegion());
    CPPUNIT_ASSERT(back->GetSpacing() == image->GetSpacing());
    CPPUNIT_ASSERT(back->GetOrigin() == image->GetOrigin());
    CPPUNIT_ASSERT(back->GetDirection() == image->GetDirection());
    const auto numberOfValues = size[0] * size[1] * size[2] * numberOfVolumes;
    CPPUNIT_ASSERT(std::equal(image->GetBufferPointer(), image->GetBufferPointer() + numberOfValues, back->GetBufferPointer()));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageHelperFunctions)
