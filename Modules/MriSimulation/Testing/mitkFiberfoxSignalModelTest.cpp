/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <mitkTestFixture.h>
#include <mitkIOUtil.h>
#include <mitkFiberfoxParameters.h>
#include <mitkRawShModel.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itksys/SystemTools.hxx>
#include <boost/property_tree/xml_parser.hpp>
#include <vnl/algo/vnl_svd.h>

/**
 * Fiberfox signal models. The prototype (raw-SH) model gets kernels whose signal a + b (u.p)^2 is
 * an exact order-2 expansion, so the turned kernel has a closed form: a + b (g.f)^2 for fiber f.
 * MITK's SH basis is evaluated in float precision, hence the tolerance of 1e-5 (the basis the
 * model used before deviated by about 0.5).
 */
class mitkFiberfoxSignalModelTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxSignalModelTestSuite);
  MITK_TEST(PrototypeFollowsFiber);
  MITK_TEST(PrototypeDirectionFromCoefficients);
  MITK_TEST(PrototypeSurvivesCopy);
  MITK_TEST(PrototypeSaveLoad);
  CPPUNIT_TEST_SUITE_END();

  using ModelType = mitk::RawShModel<double>;
  using GradientType = ModelType::GradientType;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;
  ModelType::GradientListType m_Gradients;

  GradientType RandomDirection()
  {
    GradientType v;
    for (int i = 0; i < 3; ++i)
      v[i] = m_RandGen->GetNormalVariate();
    v.Normalize();
    return v;
  }

  static GradientType Direction(double x, double y, double z)
  {
    GradientType v;
    v[0] = x; v[1] = y; v[2] = z;
    v.Normalize();
    return v;
  }

  /** Coefficients of a + b (u.p)^2, least squares on many directions (exact: the function is of order 2). */
  vnl_vector<double> KernelCoefficients(double a, double b, const GradientType& p)
  {
    const unsigned int n = 300;
    vnl_matrix<double> spherical(3, n);
    vnl_vector<double> values(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      const GradientType u = RandomDirection();
      double s[3];
      mitk::gradients::Cart2Sph(u[0], u[1], u[2], s);
      for (int r = 0; r < 3; ++r)
        spherical(r, i) = s[r];
      values[i] = a + b * (u * p) * (u * p);
    }
    const vnl_matrix<double> basis = mitk::sh::ComputeSphericalHarmonicsBasis(spherical, 2);
    return vnl_svd<double>(basis).solve(values);
  }

  static double Expected(double a, double b, const GradientType& g, const GradientType& f)
  {
    GradientType u = g;
    u.Normalize();
    return a + b * (u * f) * (u * f);
  }

public:

  void setUp() override
  {
    m_RandGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    m_RandGen->SetSeed(42);
    m_Gradients.clear();
    GradientType b0;
    b0.Fill(0.0);
    m_Gradients.push_back(b0);
    for (int i = 0; i < 40; ++i)
      m_Gradients.push_back(RandomDirection());
  }

  void tearDown() override { m_RandGen = nullptr; }

  void PrototypeFollowsFiber()
  {
    const GradientType p = Direction(1, 2, 3);
    ModelType model;
    model.SetGradientList(m_Gradients);
    model.AddKernel(KernelCoefficients(0.8, -0.5, p), 0.9, p);

    std::vector<GradientType> fibers = { p, p * -1.0, Direction(1, 0, 0), Direction(-3, 0, 1) };
    for (int i = 0; i < 50; ++i)
      fibers.push_back(RandomDirection());
    double maxError = 0;
    for (auto f : fibers)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.9, model.SimulateMeasurement(0, f), 1e-12);
      for (unsigned int g = 1; g < m_Gradients.size(); ++g)
        maxError = std::max(maxError, std::fabs(model.SimulateMeasurement(g, f) - Expected(0.8, -0.5, m_Gradients[g], f)));
    }
    MITK_INFO << "Largest deviation of the turned kernel: " << maxError;
    CPPUNIT_ASSERT_MESSAGE("the turned kernel equals the closed form", maxError < 1e-5);
  }

  void PrototypeDirectionFromCoefficients()
  {
    for (auto p : { Direction(1, 2, 3), Direction(0, 0, 1), Direction(1, 0, 0), RandomDirection() })
    {
      ModelType model;
      model.AddKernel(KernelCoefficients(0.8, -0.5, p), 1.0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::fabs(model.GetPrototypeDirection(0) * p), 1e-9);
    }
  }

  /** Parameters with a prototype fiber compartment of two kernels. */
  mitk::FiberfoxParameters PrototypeParameters()
  {
    mitk::FiberfoxParameters parameters;
    auto model = new ModelType();
    model->AddKernel(KernelCoefficients(0.8, -0.5, Direction(1, 2, 3)), 0.9, Direction(1, 2, 3));
    model->AddKernel(KernelCoefficients(0.7, -0.3, Direction(0, -1, 1)), 0.95, Direction(0, -1, 1));
    model->m_CompartmentId = 1;
    delete parameters.m_FiberModelList[0];
    parameters.m_FiberModelList[0] = model;
    parameters.SetGradienDirections(m_Gradients);
    return parameters;
  }

  void PrototypeSurvivesCopy()
  {
    mitk::FiberfoxParameters copy;
    {
      const mitk::FiberfoxParameters parameters = PrototypeParameters();
      copy = mitk::FiberfoxParameters(parameters);
    }
    auto model = dynamic_cast<ModelType*>(copy.m_FiberModelList[0]);
    CPPUNIT_ASSERT(model != nullptr);
    CPPUNIT_ASSERT_EQUAL(2u, model->GetNumberOfKernels());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, model->GetPrototypeDirection(1) * Direction(0, -1, 1), 1e-12);
    GradientType f = Direction(1, 2, 3);
    model->SetSeed(0);
    const double value = model->SimulateMeasurement(3, f);
    CPPUNIT_ASSERT_MESSAGE("the copy simulates with one of its kernels",
                           std::fabs(value - Expected(0.8, -0.5, m_Gradients[3], f)) < 1e-5 || std::fabs(value - Expected(0.7, -0.3, m_Gradients[3], f)) < 1e-5);
  }

  void PrototypeSaveLoad()
  {
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_prototype_XXXXXX");
    const std::string file = folder + "/prototype.ffp";
    PrototypeParameters().SaveParameters(file);

    for (bool withDirections : { true, false })
    {
      if (!withDirections)
      {
        // .ffp files of earlier MITK versions store no kernel direction
        boost::property_tree::ptree tree;
        boost::property_tree::xml_parser::read_xml(file, tree);
        for (auto& kernel : tree.get_child("fiberfox.image.compartments.c0.kernels"))
          kernel.second.erase("dir");
        boost::property_tree::xml_parser::write_xml(file, tree);
      }
      mitk::FiberfoxParameters loaded;
      loaded.LoadParameters(file, true);
      auto model = dynamic_cast<ModelType*>(loaded.m_FiberModelList.at(0));
      CPPUNIT_ASSERT(model != nullptr);
      CPPUNIT_ASSERT_EQUAL(2u, model->GetNumberOfKernels());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::fabs(model->GetPrototypeDirection(0) * Direction(1, 2, 3)), 1e-9);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::fabs(model->GetPrototypeDirection(1) * Direction(0, -1, 1)), 1e-9);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.95, model->GetBaselineSignal(1), 1e-12);
    }
    itksys::SystemTools::RemoveADirectory(folder);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxSignalModel)
