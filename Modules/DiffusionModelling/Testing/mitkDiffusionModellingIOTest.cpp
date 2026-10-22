/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkIOUtil.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itksys/SystemTools.hxx>
#include <vnl/vnl_matrix_fixed.h>

#include <cmath>
#include <fstream>
#include <sstream>

class mitkDiffusionModellingIOTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingIOTestSuite);
  MITK_TEST(TensorReader_StandardFrame_RotatedToWorld);
  MITK_TEST(TensorReader_MitkFrameKey_RotatedToWorld);
  MITK_TEST(TensorReader_PermutationFrame_RotatedToWorld);
  MITK_TEST(TensorReader_NoFrame_Unchanged);
  CPPUNIT_TEST_SUITE_END();

  using Matrix3 = vnl_matrix_fixed<double, 3, 3>;

  std::string m_TempDir;

  /** Writes a one-voxel ASCII NRRD tensor file (elements xx, xy, xz, yy, yz, zz) with the given header lines. */
  std::string WriteTensorNrrd(const std::string &name, const double (&elements)[6], const std::string &frameLines)
  {
    const std::string path = m_TempDir + "/" + name;
    std::ofstream file(path);
    file << "NRRD0004\n"
         << "type: float\n"
         << "dimension: 4\n"
         << "space: left-posterior-superior\n"
         << "sizes: 6 1 1 1\n"
         << "space directions: none (1,0,0) (0,1,0) (0,0,1)\n"
         << "kinds: 3D-symmetric-matrix domain domain domain\n"
         << "endian: little\n"
         << "encoding: ascii\n"
         << "space origin: (0,0,0)\n"
         << frameLines << (frameLines.empty() ? "" : "\n") << "\n";
    for (double e : elements)
      file << e << " ";
    file << "\n";
    return path;
  }

  Matrix3 ReadTensor(const std::string &path)
  {
    auto image = mitk::IOUtil::Load<mitk::TensorImage>(path);
    mitk::ImagePixelReadAccessor<mitk::TensorImage::PixelType, 3> access(image.GetPointer());
    itk::Index<3> idx = {{0, 0, 0}};
    auto tensor = access.GetPixelByIndex(idx);
    Matrix3 m;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        m(r, c) = tensor(r, c);
    return m;
  }

  static Matrix3 ToMatrix(const double (&e)[6])
  {
    Matrix3 m;
    m(0, 0) = e[0]; m(0, 1) = e[1]; m(0, 2) = e[2];
    m(1, 0) = e[1]; m(1, 1) = e[3]; m(1, 2) = e[4];
    m(2, 0) = e[2]; m(2, 1) = e[4]; m(2, 2) = e[5];
    return m;
  }

  static void AssertEqual(const Matrix3 &expected, const Matrix3 &actual)
  {
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected(r, c), actual(r, c), 1e-5);
  }

  /** The frame vectors are the rows of F; the world tensor is F^T T F, as for DWI gradients. */
  void CheckFrame(const std::string &name, const Matrix3 &frameRows, bool mitkKey)
  {
    const double elements[6] = {3, 0.2, 0, 1, 0, 0.5};
    std::ostringstream frame;
    frame.precision(17);
    if (mitkKey)
      frame << "measurement frame:=";
    else
      frame << "measurement frame: ";
    for (unsigned int r = 0; r < 3; ++r)
      frame << (mitkKey ? " ( " : "(") << frameRows(r, 0) << (mitkKey ? " , " : ",") << frameRows(r, 1) << (mitkKey ? " , " : ",")
            << frameRows(r, 2) << (mitkKey ? " )" : ") ");

    const Matrix3 expected = frameRows.transpose() * ToMatrix(elements) * frameRows;
    AssertEqual(expected, ReadTensor(WriteTensorNrrd(name, elements, frame.str())));
  }

  static Matrix3 RotationZ30()
  {
    const double c = std::cos(itk::Math::pi / 6);
    const double s = std::sin(itk::Math::pi / 6);
    Matrix3 f;
    f(0, 0) = c;  f(0, 1) = s; f(0, 2) = 0;
    f(1, 0) = -s; f(1, 1) = c; f(1, 2) = 0;
    f(2, 0) = 0;  f(2, 1) = 0; f(2, 2) = 1;
    return f;
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("mitkDiffusionModellingIOTest_XXXXXX");
  }

  void tearDown() override
  {
    itksys::SystemTools::RemoveADirectory(m_TempDir);
  }

  void TensorReader_StandardFrame_RotatedToWorld()
  {
    CheckFrame("standard.dti", RotationZ30(), false);
  }

  void TensorReader_MitkFrameKey_RotatedToWorld()
  {
    CheckFrame("mitk.dti", RotationZ30(), true);
  }

  void TensorReader_PermutationFrame_RotatedToWorld()
  {
    Matrix3 f;
    f.fill(0);
    f(0, 1) = 1;
    f(1, 2) = 1;
    f(2, 0) = 1;
    CheckFrame("permutation.dti", f, false);
  }

  void TensorReader_NoFrame_Unchanged()
  {
    const double elements[6] = {3, 0.2, 0, 1, 0, 0.5};
    AssertEqual(ToMatrix(elements), ReadTensor(WriteTensorNrrd("noframe.dti", elements, "")));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingIO)
