/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkFiberWeightsFile.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkExceptionMacro.h>

#include "mitkFiberProcessingTestHelper.h"

#include <sstream>

class mitkFiberWeightsFileTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberWeightsFileTestSuite);
  MITK_TEST(Read_MRtrixHeader_SkipsComments);
  MITK_TEST(Read_OneWeightPerLine);
  MITK_TEST(Read_InvalidValue_Throws);
  MITK_TEST(Read_MissingFile_Throws);
  MITK_TEST(Apply_MatchingCount_SetsWeightsInOrder);
  MITK_TEST(Apply_CountMismatch_ThrowsAndKeepsWeights);
  CPPUNIT_TEST_SUITE_END();

  mitk::FiberBundle::Pointer ThreeFibers()
  {
    using namespace mitk::FiberProcessingTestHelper;
    return CreateBundle({Line({0, 0, 0}, {10, 0, 0}, 5), Line({0, 1, 0}, {10, 1, 0}, 5), Line({0, 2, 0}, {10, 2, 0}, 5)});
  }

public:
  void Read_MRtrixHeader_SkipsComments()
  {
    // tcksift2 3.0 writes a command history comment line before the weights.
    std::istringstream stream("# command_history: tcksift2 in.tck fod.mif weights.txt  (version=3.0.8)\n"
                              "0.01009899323 1.5 2e-3\n");
    const auto weights = mitk::FiberWeightsFile::Read(stream);
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), weights.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.01009899323, weights[0], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, weights[1], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.002, weights[2], 1e-9);
  }

  void Read_OneWeightPerLine()
  {
    std::istringstream stream("0.5\r\n\r\n  0.25\r\n# comment\n4\n");
    const auto weights = mitk::FiberWeightsFile::Read(stream);
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), weights.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, weights[0], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, weights[1], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, weights[2], 1e-9);
  }

  void Read_InvalidValue_Throws()
  {
    std::istringstream text("0.5 abc 1\n");
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Read(text), mitk::Exception);
    std::istringstream trailing("0.5 1.0x\n");
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Read(trailing), mitk::Exception);
    std::istringstream notFinite("0.5 nan\n");
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Read(notFinite), mitk::Exception);
  }

  void Read_MissingFile_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Read(std::string("/nonexistent/sift2_weights.txt")), mitk::Exception);
  }

  void Apply_MatchingCount_SetsWeightsInOrder()
  {
    auto fib = ThreeFibers();
    mitk::FiberWeightsFile::Apply({0.1f, 0.2f, 0.3f}, fib);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, fib->GetFiberWeight(0), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2, fib->GetFiberWeight(1), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, fib->GetFiberWeight(2), 1e-6);
  }

  void Apply_CountMismatch_ThrowsAndKeepsWeights()
  {
    auto fib = ThreeFibers();
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Apply({0.1f, 0.2f}, fib), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::FiberWeightsFile::Apply({0.1f, 0.2f, 0.3f, 0.4f}, fib), mitk::Exception);
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fib->GetFiberWeight(i), 1e-6);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberWeightsFile)
