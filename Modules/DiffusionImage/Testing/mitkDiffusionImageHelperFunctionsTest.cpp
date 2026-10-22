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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionImageHelperFunctions)
