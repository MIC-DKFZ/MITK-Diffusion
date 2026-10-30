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
#include <mitkAstroStickModel.h>
#include <mitkNoncentralChiNoiseModel.h>
#include <cmath>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
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
  MITK_TEST(SampleKernels_B0WithShortVector_TensorFitTakesItAsB0);
  MITK_TEST(SampleKernels_FirstB0ByTheRule_NormalizesKernels);
  MITK_TEST(AstroSticksVector);
  MITK_TEST(NoncentralChiNoiseMoments);
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

  /**
   * One-voxel DWI, reference b-value 1000, of a tensor with eigenvalues 1.7e-3, 0.3e-3, 0.3e-3
   * along x (FA 0.80): S_i = s0[i] exp(-1000 g_i^T D g_i), rounded to short.
   */
  static mitk::Image::Pointer TensorDwi(const std::vector<mitk::DiffusionPropertyHelper::GradientDirectionType>& gradients, const std::vector<double>& s0)
  {
    auto container = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
    auto dwi = mitk::DiffusionPropertyHelper::ImageType::New();
    mitk::DiffusionPropertyHelper::ImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    dwi->SetRegions(region);
    dwi->SetVectorLength(gradients.size());
    dwi->Allocate();
    mitk::DiffusionPropertyHelper::ImageType::PixelType pixel(gradients.size());
    for (unsigned int i = 0; i < gradients.size(); ++i)
    {
      const auto& g = gradients[i];
      container->push_back(g);
      const double gDg = 1.7e-3 * g[0] * g[0] + 0.3e-3 * (g[1] * g[1] + g[2] * g[2]);
      pixel[i] = static_cast<short>(std::round(s0[i] * std::exp(-1000 * gDg)));
    }
    dwi->SetPixel({{0, 0, 0}}, pixel);
    auto image = mitk::GrabItkImageMemory(dwi.GetPointer());
    mitk::DiffusionPropertyHelper::MeasurementFrameType identity;
    identity.set_identity();
    mitk::DiffusionPropertyHelper::InitializeImage(image, container.GetPointer(), 1000, identity, false);
    return image;
  }

  /** Six unit directions, x and y and z first. */
  static std::vector<mitk::DiffusionPropertyHelper::GradientDirectionType> SixDirections()
  {
    const double s = std::sqrt(0.5);
    return { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {s, s, 0}, {s, 0, s}, {0, s, s} };
  }

  /** A b0 volume by MITK's rule: b = 1000 * |g|^2 = 0.3 rounds to 0. Along x, the fiber direction. */
  static mitk::DiffusionPropertyHelper::GradientDirectionType ShortB0Vector()
  {
    return { std::sqrt(0.0003), 0, 0 };
  }

  /** Samples the kernels of the DWI with FA and ADC ranges around the tensor's (FA 0.80, ADC 0.77e-3). */
  static unsigned int SampledKernels(mitk::Image::Pointer dwi, ModelType& model)
  {
    model.SetFaRange(0.7, 0.9);
    model.SetAdcRange(0.0005, 0.001);
    model.SampleKernels(dwi, nullptr);
    return model.GetNumberOfKernels();
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

  /**
   * The tensor that selects prototype voxels by FA treats the b = 0.3 volume as b0, as the b0
   * normalization, the Q-ball and the ADC fit do. ITK's tensor fit would take the short vector
   * as a volume at the full b-value with no attenuation along x, and the FA would leave the range.
   */
  void SampleKernels_B0WithShortVector_TensorFitTakesItAsB0()
  {
    auto gradients = SixDirections();
    gradients.insert(gradients.begin(), ShortB0Vector());
    gradients.insert(gradients.begin(), mitk::DiffusionPropertyHelper::GradientDirectionType(0.0));
    ModelType model;
    CPPUNIT_ASSERT_EQUAL(1u, SampledKernels(TensorDwi(gradients, std::vector<double>(gradients.size(), 1000)), model));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::fabs(model.GetPrototypeDirection(0) * Direction(1, 0, 0)), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, model.GetBaselineSignal(0), 1e-12);
  }

  /**
   * The first b0 volume by MITK's rule (here the b = 0.3 volume, index 0) normalizes the kernels; a
   * voxel with a volume above it is skipped. The exact zero vector at index 1 carries a lower
   * baseline (noise), so taking it as the reference would skip the only voxel.
   */
  void SampleKernels_FirstB0ByTheRule_NormalizesKernels()
  {
    auto gradients = SixDirections();
    gradients.insert(gradients.begin(), mitk::DiffusionPropertyHelper::GradientDirectionType(0.0));
    gradients.insert(gradients.begin(), ShortB0Vector());
    std::vector<double> s0(gradients.size(), 1000);
    s0[1] = 990;
    ModelType model;
    CPPUNIT_ASSERT_EQUAL(1u, SampledKernels(TensorDwi(gradients, s0), model));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, model.GetBaselineSignal(0), 1e-12);
  }

  /** The vector overload of the astrosticks model equals the per-direction one (it summed into uninitialized memory). */
  void AstroSticksVector()
  {
    mitk::AstroStickModel<double> model;
    model.SetGradientList(m_Gradients);
    GradientType f = Direction(1, 0, 0);
    const auto signal = model.SimulateMeasurement(f);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(m_Gradients.size()), signal.Size());
    for (unsigned int g = 0; g < m_Gradients.size(); ++g)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(model.SimulateMeasurement(g, f), signal[g], 1e-12);
  }

  /**
   * Magnitude noise of n coils combined by sum of squares, variance s2 per real and imaginary
   * part, follows a noncentral chi distribution: E[M^2] = 2 n s2 + S^2, at S = 0
   * E[M] = sqrt(2 s2) Gamma(n + 1/2) / Gamma(n), and at high signal E[M] = sqrt(S^2 + (2n - 1) s2)
   * with variance s2. One coil is Rician.
   */
  void NoncentralChiNoiseMoments()
  {
    const double variance = 100;
    const unsigned int samples = 200000;
    for (unsigned int coils : { 1u, 4u, 8u })
      for (double signal : { 0.0, 1000.0 })
      {
        mitk::NoncentralChiNoiseModel<double> model;
        model.SetSeed(1);
        model.SetNoiseVariance(variance);
        model.SetNumberOfCoils(coils);
        mitk::NoncentralChiNoiseModel<double>::PixelType pixel(samples);
        pixel.Fill(signal);
        model.AddNoise(pixel);
        double sum = 0, sumOfSquares = 0;
        for (unsigned int i = 0; i < samples; ++i)
        {
          sum += pixel[i];
          sumOfSquares += pixel[i] * pixel[i];
        }
        const double mean = sum / samples;
        const double meanSquare = sumOfSquares / samples;
        const double expectedMeanSquare = 2 * coils * variance + signal * signal;
        const double expectedMean = signal == 0
          ? std::sqrt(2 * variance) * std::exp(std::lgamma(coils + 0.5) - std::lgamma(coils))
          : std::sqrt(signal * signal + (2.0 * coils - 1) * variance);
        const double expectedVariance = expectedMeanSquare - expectedMean * expectedMean;
        // M^2 / s2 is noncentral chi-square with 2n degrees of freedom and noncentrality S^2 / s2
        const double varianceOfSquare = 4 * coils * variance * variance + 4 * signal * signal * variance;
        MITK_INFO << "Noncentral chi noise, " << coils << " coils, signal " << signal << ": mean " << mean
                  << " (expected " << expectedMean << "), variance " << meanSquare - mean * mean
                  << " (expected " << expectedVariance << "), mean square " << meanSquare << " (expected " << expectedMeanSquare << ")";
        // five standard errors of the sample mean and of the sample mean square
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedMean, mean, 5 * std::sqrt(expectedVariance / samples));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedMeanSquare, meanSquare, 5 * std::sqrt(varianceOfSquare / samples));
      }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxSignalModel)
