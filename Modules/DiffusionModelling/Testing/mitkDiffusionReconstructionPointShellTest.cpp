/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkPointShell.h>

#include <algorithm>
#include <cmath>
#include <vector>

/**
 * Guards the sampling tables of the sphere. Their values and their order are part of the ODF file
 * format (one value per direction) and of Fiberfox's stick directions, so they must not change.
 */
class mitkDiffusionReconstructionPointShellTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionPointShellTestSuite);
  MITK_TEST(Tables_MatchChecksums);
  MITK_TEST(Tables_UnitLengthAndAntipodal);
  MITK_TEST(Directions_SharedAndFixedSizeEqual);
  MITK_TEST(Directions_UnsupportedSize_Throws);
  CPPUNIT_TEST_SUITE_END();

  /** Sum of all values, sum weighted with the position (order sensitive), first and last direction. */
  struct Checksum
  {
    unsigned int n;
    double sum;
    double weightedSum;
    double first[3];
    double last[3];
  };

  static const std::vector<Checksum> &Checksums()
  {
    // Computed from the tables as they are in use (2026-10-23).
    static const std::vector<Checksum> checksums = {
      {12, 1.9999999995301554e-06, -39.006871000000011, {0.25025399999999998, -0.91616200000000003, -0.313083}, {-0.97567599999999999, 0.177817, 0.12820400000000001}},
      {42, 9.3901971478516089e-06, -798.1471725346737, {0.97567675495554884, -0.17781613758991291, -0.1282040992012935}, {-0.89710163704047841, -0.34103744532956171, 0.28089520056843209}},
      {92, 2.0415103920701227e-05, -3303.9986165318837, {0.97567675495554873, -0.17781613758991291, -0.1282040992012935}, {-0.77447189788264303, -0.21063795660617651, 0.59651062909794961}},
      {162, 3.5868395638627426e-05, -9436.6593833276183, {0.97567675495554884, -0.17781613758991291, -0.1282040992012935}, {-0.74067547630238095, -0.38941630235214841, 0.5474986596025796}},
      {252, 5.5739721668945919e-05, -21721.136635228082, {0.97567675495554884, -0.17781613758991291, -0.12820409920129339}, {-0.70943413100008901, -0.48577126053830838, 0.51061678018566092}},
      {362, 8.0028000297938107e-05, -43369.766797419587, {0.97567675495554873, -0.17781613758991291, -0.1282040992012935}, {-0.68535824374160914, -0.54394172278153985, 0.4841605931432954}},
      {492, 0.0001087328319805958, -78282.248778216628, {0.97567675495554873, -0.17781613758991291, -0.12820409920129339}, {-0.66703583440727177, -0.58229800262959064, 0.4647496441635906}},
      {642, 0.00014185405326988532, -131045.64102464732, {0.97567675495554884, -0.17781613758991291, -0.1282040992012935}, {-0.65285131290297793, -0.60929303830661863, 0.45005239329654428}},
      {812, 0.00017939158758745188, -206934.36060032054, {0.97567675495554873, -0.17781613758991291, -0.12820409920129339}, {-0.6416277774997261, -0.62924386513043207, 0.4385954324163403}},
      {1002, 0.00022134539412632659, -311910.18274749577, {0.97567675495554884, -0.17781613758991291, -0.12820409920129339}, {-0.63256094633541793, -0.64455332568219015, 0.42943877272932518}}};
    return checksums;
  }

public:
  void Tables_MatchChecksums()
  {
    for (const auto &expected : Checksums())
    {
      const vnl_matrix<double> &table = mitk::pointshell::Directions(expected.n);
      CPPUNIT_ASSERT_EQUAL(3u, table.rows());
      CPPUNIT_ASSERT_EQUAL(expected.n, table.cols());
      double sum = 0;
      double weightedSum = 0;
      for (unsigned int i = 0; i < table.cols(); ++i)
      {
        for (unsigned int r = 0; r < 3; ++r)
        {
          sum += table(r, i);
          weightedSum += table(r, i) * (3 * i + r + 1);
        }
      }
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.sum, sum, 1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.weightedSum, weightedSum, 1e-9);
      for (unsigned int r = 0; r < 3; ++r)
      {
        CPPUNIT_ASSERT_EQUAL(expected.first[r], table(r, 0));
        CPPUNIT_ASSERT_EQUAL(expected.last[r], table(r, expected.n - 1));
      }
    }
  }

  // The tables are unit vectors in antipodal pairs (to the precision they were written with).
  void Tables_UnitLengthAndAntipodal()
  {
    for (const auto &expected : Checksums())
    {
      const vnl_matrix<double> &table = mitk::pointshell::Directions(expected.n);
      for (unsigned int i = 0; i < table.cols(); ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, table.get_column(i).two_norm(), 2e-6);
        double nearestAntipode = 10;
        for (unsigned int j = 0; j < table.cols(); ++j)
          nearestAntipode = std::min(nearestAntipode, (table.get_column(i) + table.get_column(j)).inf_norm());
        CPPUNIT_ASSERT(nearestAntipode < 2e-6);
      }
    }
  }

  // Each table is built once: every call returns the same object, and the fixed-size variant has
  // the same values.
  void Directions_SharedAndFixedSizeEqual()
  {
    CPPUNIT_ASSERT(&mitk::pointshell::Directions(252) == &mitk::pointshell::Directions(252));
    CPPUNIT_ASSERT(&mitk::pointshell::Directions<252>() == &mitk::pointshell::Directions<252>());
    CPPUNIT_ASSERT(mitk::pointshell::Directions<42>().as_matrix() == mitk::pointshell::Directions(42));
    CPPUNIT_ASSERT_EQUAL(Checksums().size(), mitk::pointshell::SupportedSizes().size());
    for (unsigned int i = 0; i < Checksums().size(); ++i)
      CPPUNIT_ASSERT_EQUAL(Checksums()[i].n, mitk::pointshell::SupportedSizes()[i]);
  }

  void Directions_UnsupportedSize_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::pointshell::Directions(100), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionPointShell)
