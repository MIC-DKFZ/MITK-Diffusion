/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <mitkTestFixture.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <set>
#include <string>
#include <vector>
#include <mitkFiberfoxParameters.h>
#include <mitkBallModel.h>
#include <itkTractsToDWIImageFilter.h>
#include <itkKspaceImageFilter.h>
#include <itkFieldmapGeneratorFilter.h>
#include <mitkSingleShotEpi.h>
#include <mitkConventionalSpinEcho.h>
#include <mitkFastSpinEcho.h>

/** The k-space filter with the count of time segments of its last update readable (0: exact sum). */
class KspaceSegmentsProbe : public itk::KspaceImageFilter<float>
{
public:
  using Self = KspaceSegmentsProbe;
  using Pointer = itk::SmartPointer<Self>;
  itkFactorylessNewMacro(Self)
  unsigned int NumSegments() const { return m_NumSegments; }
};

/**
 * k-space acquisition of Fiberfox on synthetic input images (the input-image mode of the
 * simulation: no signal model, no relaxation, signal scale 1): each artifact against its
 * definition, the noise statistics, the readout schemes and the field map generator. Needs no
 * mitk-diffusion-test-data. The definitions are those of the method check of review #13.
 */
class mitkFiberfoxKspaceTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxKspaceTestSuite);
  MITK_TEST(RoundTrip);
  MITK_TEST(GibbsRingingIsIdealLowPass);
  MITK_TEST(NyquistGhost);
  MITK_TEST(Aliasing);
  MITK_TEST(OffResonanceShift);
  MITK_TEST(EddyCurrentsEqualFieldMap);
  MITK_TEST(TimeSegmentsApproximateExactSum);
  MITK_TEST(TimeSegmentsHoldForStrongPhase);
  MITK_TEST(NoiseFullFourier);
  MITK_TEST(NoisePartialFourier);
  MITK_TEST(NoiseMultiCoil);
  MITK_TEST(NoKspaceNoiseWithNoiseModel);
  MITK_TEST(ReadoutSchemes);
  MITK_TEST(FieldmapGenerator);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::TractsToDWIImageFilter<short>;
  using DwiType = FilterType::OutputImageType;
  using DoubleDwiType = FilterType::DoubleDwiType;
  using FloatImageType = itk::Image<float, 3>;
  using Complex = std::complex<double>;
  using Pattern = std::function<double(int x, int y)>;

  /** Input DWI: every slice and volume holds pattern(x, y). */
  static DwiType::Pointer MakeImage(int nx, int ny, int nz, unsigned int numVolumes, const Pattern& pattern)
  {
    auto image = DwiType::New();
    itk::ImageRegion<3> region;
    region.SetSize(0, nx);
    region.SetSize(1, ny);
    region.SetSize(2, nz);
    image->SetRegions(region);
    itk::Vector<double, 3> spacing;
    spacing.Fill(2.0);
    image->SetSpacing(spacing);
    image->SetVectorLength(numVolumes);
    image->Allocate();
    for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
        {
          DwiType::PixelType pixel(numVolumes);
          pixel.Fill(static_cast<short>(pattern(x, y)));
          image->SetPixel({{x, y, z}}, pixel);
        }
    return image;
  }

  /** Acquisition of the input image without artifacts; the gradients are b = 0 except where given. */
  static mitk::FiberfoxParameters MakeParameters(const DwiType* image, mitk::SignalGenerationParameters::GradientListType gradients = {})
  {
    mitk::FiberfoxParameters parameters;
    parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    parameters.m_SignalGen.m_SignalScale = 1;
    parameters.m_SignalGen.m_ImageRegion = image->GetLargestPossibleRegion();
    parameters.m_SignalGen.m_ImageSpacing = image->GetSpacing();
    parameters.m_SignalGen.m_ImageOrigin = image->GetOrigin();
    parameters.m_SignalGen.m_ImageDirection = image->GetDirection();
    parameters.m_Misc.m_OutputPath = "";  // no log file
    // the exact sum of the phase of distortions and eddy currents, the model the cases below check;
    // TimeSegmentsApproximateExactSum compares the default, the time-segmented phase, with it
    parameters.m_SignalGen.m_KspaceTimeSegments = 0;
    // the k-space simulation reads T1/T2 of one model per compartment image
    parameters.m_NonFiberModelList.push_back(new mitk::BallModel<double>());
    gradients.resize(image->GetVectorLength());
    parameters.SetGradienDirections(gradients);
    return parameters;
  }

  static FilterType::Pointer Acquire(const DwiType::Pointer& image, const mitk::FiberfoxParameters& parameters)
  {
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(true);
    filter->SetParameters(parameters);
    filter->SetInputImage(image);
    filter->Update();
    return filter;
  }

  /** Complex image of coil c, volume g, slice z. */
  static Complex Value(const FilterType* filter, int x, int y, int z = 0, unsigned int g = 0, unsigned int c = 0)
  {
    const itk::Index<3> index = {{x, y, z}};
    return Complex(filter->GetOutputImagesReal().at(c)->GetPixel(index)[g], filter->GetOutputImagesImag().at(c)->GetPixel(index)[g]);
  }

  /**
   * Fiberfox's centred DFT pair (KspaceImageFilter forward with 1/N, DftImageFilter inverse; index
   * n of N at n - N/2), with the k-space samples for which keep(kx index, ky index) is false set to 0.
   */
  static std::vector<Complex> KeepFrequencies(const Pattern& f, int nx, int ny, const std::function<bool(int, int)>& keep)
  {
    const double pi2 = 2 * itk::Math::pi;
    std::vector<Complex> k(nx * ny, 0.0);
    for (int ky = 0; ky < ny; ++ky)
      for (int kx = 0; kx < nx; ++kx)
      {
        if (!keep(kx, ky))
          continue;
        Complex s = 0;
        for (int y = 0; y < ny; ++y)
          for (int x = 0; x < nx; ++x)
            s += f(x, y) * std::exp(Complex(0, pi2 * (double(kx - nx/2) * (x - nx/2) / nx + double(ky - ny/2) * (y - ny/2) / ny)));
        k[kx + nx * ky] = s / double(nx * ny);
      }
    std::vector<Complex> image(nx * ny, 0.0);
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
      {
        Complex s = 0;
        for (int ky = 0; ky < ny; ++ky)
          for (int kx = 0; kx < nx; ++kx)
            s += k[kx + nx * ky] * std::exp(Complex(0, -pi2 * (double(kx - nx/2) * (x - nx/2) / nx + double(ky - ny/2) * (y - ny/2) / ny)));
        image[x + nx * y] = s;
      }
    return image;
  }

  /** A box with a gradient inside, so that every voxel value is distinct. */
  static Pattern Box(int x0, int x1, int y0, int y1)
  {
    return [=](int x, int y) { return (x >= x0 && x < x1 && y >= y0 && y < y1) ? 500.0 + 7 * x + 3 * y : 0.0; };
  }

  /** Centre of mass of the magnitude in x and y, in the lines y0 to y1 - 1. */
  static std::array<double, 2> CentreOfMass(const FilterType* filter, int nx, int y0, int y1)
  {
    double sum = 0, cx = 0, cy = 0;
    for (int y = y0; y < y1; ++y)
      for (int x = 0; x < nx; ++x)
      {
        const double m = std::abs(Value(filter, x, y));
        sum += m;
        cx += m * x;
        cy += m * y;
      }
    return {cx / sum, cy / sum};
  }

  /** Variance of the real and imaginary parts and mean magnitude of all voxels, volumes and slices of coil 0. */
  static void NoiseStatistics(const FilterType* filter, double& varReal, double& varImag, double& meanMagnitude)
  {
    auto real = filter->GetOutputImagesReal().at(0);
    auto imag = filter->GetOutputImagesImag().at(0);
    const std::size_t n = real->GetLargestPossibleRegion().GetNumberOfPixels() * real->GetNumberOfComponentsPerPixel();
    double sumR = 0, sumI = 0, sumR2 = 0, sumI2 = 0, sumM = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      const double r = real->GetBufferPointer()[i];
      const double m = imag->GetBufferPointer()[i];
      sumR += r; sumI += m; sumR2 += r * r; sumI2 += m * m;
      sumM += std::sqrt(r * r + m * m);
    }
    varReal = sumR2 / n - (sumR / n) * (sumR / n);
    varImag = sumI2 / n - (sumI / n) * (sumI / n);
    meanMagnitude = sumM / n;
  }

public:

  /** Without artifacts the acquisition returns the input: real part = input, imaginary part 0. */
  void RoundTrip()
  {
    for (auto size : { std::array<int, 2>{32, 32}, std::array<int, 2>{31, 29} })
    {
      const int nx = size[0], ny = size[1];
      const Pattern f = [](int x, int y) { return 100.0 + 10 * x + 3 * y + (x * y) % 7; };
      auto image = MakeImage(nx, ny, 1, 2, f);
      auto filter = Acquire(image, MakeParameters(image));
      double maxError = 0;
      for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
        {
          maxError = std::max(maxError, std::abs(Value(filter, x, y, 0, 1) - Complex(f(x, y), 0)));
          CPPUNIT_ASSERT_EQUAL(static_cast<short>(f(x, y)), filter->GetOutput()->GetPixel({{x, y, 0}})[1]);
        }
      MITK_INFO << nx << "x" << ny << " round trip: largest deviation " << maxError;
      CPPUNIT_ASSERT_MESSAGE("input reproduced", maxError < 0.05);
    }
  }

  /** Gibbs ringing by zeroing p % of the frequencies is the ideal low-pass: the outer ceil(N/2 p/100) lines on each side are 0. */
  void GibbsRingingIsIdealLowPass()
  {
    const int nx = 32, ny = 30, percent = 50;
    const Pattern f = Box(8, 22, 10, 21);
    auto image = MakeImage(nx, ny, 1, 1, f);
    auto parameters = MakeParameters(image);
    parameters.m_SignalGen.m_DoAddGibbsRinging = true;
    parameters.m_SignalGen.m_ZeroRinging = percent;
    auto filter = Acquire(image, parameters);

    const int rx = static_cast<int>(std::ceil(nx / 2.0 * percent / 100.0));
    const int ry = static_cast<int>(std::ceil(ny / 2.0 * percent / 100.0));
    const auto expected = KeepFrequencies(f, nx, ny, [=](int kx, int ky) { return kx >= rx && kx < nx - rx && ky >= ry && ky < ny - ry; });
    double maxError = 0, overshoot = 0;
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
      {
        maxError = std::max(maxError, std::abs(Value(filter, x, y) - expected[x + nx * y]));
        overshoot = std::max(overshoot, std::abs(expected[x + nx * y]) - f(x, y));
      }
    MITK_INFO << "Gibbs ringing: largest deviation from the ideal low-pass " << maxError << ", overshoot " << overshoot;
    CPPUNIT_ASSERT_MESSAGE("ringing present", overshoot > 20);
    CPPUNIT_ASSERT_MESSAGE("ideal low-pass", maxError < 0.05);
  }

  /**
   * N/2 ghost: readout lines shifted by +d and -d in alternation give the image f cos(t) plus the
   * image shifted by N/2 in phase direction times i sin(t), t = 2 pi d x/Nx (x from the centre), so
   * |image|^2 = (f(x,y) cos t)^2 + (f(x,y+N/2) sin t)^2.
   */
  void NyquistGhost()
  {
    const int nx = 32, ny = 32;
    const double d = 0.25;
    const Pattern f = Box(6, 26, 3, 14);
    auto image = MakeImage(nx, ny, 1, 1, f);
    auto parameters = MakeParameters(image);
    parameters.m_Misc.m_DoAddGhosts = true;
    parameters.m_SignalGen.m_KspaceLineOffset = d;
    auto filter = Acquire(image, parameters);

    double maxError = 0, maxGhost = 0;
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
      {
        const double t = 2 * itk::Math::pi * d * (x - nx/2) / nx;
        const double main = f(x, y) * std::cos(t);
        const double ghost = f(x, (y + ny/2) % ny) * std::sin(t);
        maxError = std::max(maxError, std::fabs(std::abs(Value(filter, x, y)) - std::sqrt(main * main + ghost * ghost)));
        maxGhost = std::max(maxGhost, std::fabs(ghost));
      }
    MITK_INFO << "N/2 ghost: largest deviation " << maxError << ", ghost up to " << maxGhost;
    CPPUNIT_ASSERT(maxGhost > 100);
    CPPUNIT_ASSERT_MESSAGE("ghost at y + N/2 with amplitude sin(2 pi d x/N)", maxError < 0.05);
  }

  /**
   * Aliasing: the field of view in phase direction is reduced to M = N * factor lines around the
   * centre, and every input line folds onto the output line congruent to it modulo M.
   */
  void Aliasing()
  {
    const int nx = 16, ny = 32;
    const Pattern f = [](int x, int y) { return 100.0 + 5 * x + 11 * y; };
    for (float factor : { 0.5f, 0.75f })
    {
      auto image = MakeImage(nx, ny, 1, 1, f);
      auto parameters = MakeParameters(image);
      parameters.m_Misc.m_DoAddAliasing = true;
      parameters.m_SignalGen.m_CroppingFactor = factor;
      auto filter = Acquire(image, parameters);

      const int m = static_cast<int>(ny * factor);
      CPPUNIT_ASSERT_EQUAL(static_cast<itk::SizeValueType>(m), filter->GetOutput()->GetLargestPossibleRegion().GetSize(1));
      std::vector<double> expected(nx * m, 0.0);
      for (int y = 0; y < ny; ++y)
      {
        const int folded = (((y - ny/2 + m/2) % m) + m) % m;
        for (int x = 0; x < nx; ++x)
          expected[x + nx * folded] += f(x, y);
      }
      double maxError = 0;
      for (int y = 0; y < m; ++y)
        for (int x = 0; x < nx; ++x)
          maxError = std::max(maxError, std::abs(Value(filter, x, y) - Complex(expected[x + nx * y], 0)));
      MITK_INFO << "Aliasing " << factor << ": largest deviation " << maxError;
      CPPUNIT_ASSERT_MESSAGE("lines outside the field of view folded back", maxError < 0.05);
    }
  }

  /**
   * A constant off-resonance f shifts an EPI image in phase direction by f * tLine * Ny pixels
   * (2 pixels here), towards lower y with the normal and towards higher y with the reversed phase
   * direction; nothing moves in readout direction.
   */
  void OffResonanceShift()
  {
    const int nx = 32, ny = 32;
    const Pattern f = Box(10, 22, 8, 20);
    auto image = MakeImage(nx, ny, 1, 1, f);
    auto input = Acquire(image, MakeParameters(image));
    // the lines of the object (8 to 19) and of its shifted copies; the alternating readout direction
    // adds a faint N/2 ghost 16 lines away, outside these lines
    const auto centre = CentreOfMass(input, nx, 6, 22);

    const double frequency = 62.5;
    for (bool reverse : { false, true })
    {
      auto parameters = MakeParameters(image);
      parameters.m_SignalGen.m_ReversePhase = reverse;
      parameters.m_Misc.m_DoAddDistortions = true;
      auto fmap = FloatImageType::New();
      fmap->SetRegions(image->GetLargestPossibleRegion());
      fmap->SetSpacing(image->GetSpacing());
      fmap->Allocate();
      fmap->FillBuffer(frequency);
      parameters.m_SignalGen.m_FrequencyMap = fmap;
      auto filter = Acquire(image, parameters);

      const auto moved = CentreOfMass(filter, nx, 6, 22);
      const double expected = (reverse ? 1 : -1) * frequency * parameters.m_SignalGen.m_tLine / 1000 * ny;
      MITK_INFO << "Off-resonance, " << (reverse ? "reversed" : "normal") << " phase direction: shift " << moved[1] - centre[1] << " (expected " << expected << "), readout " << moved[0] - centre[0];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, moved[1] - centre[1], 0.05);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, moved[0] - centre[0], 0.05);
    }
  }

  /**
   * Eddy currents of strength G (mT/m, b = 1000, without decay) along x act like the field map
   * gamma_bar G x with gamma_bar = 42.576 MHz/T.
   */
  void EddyCurrentsEqualFieldMap()
  {
    const int nx = 32, ny = 32;
    const Pattern f = Box(6, 26, 8, 24);
    mitk::SignalGenerationParameters::GradientListType gradients(2);
    gradients[1][0] = 1;
    auto image = MakeImage(nx, ny, 1, 2, f);
    const double strength = 0.02;

    auto eddyParameters = MakeParameters(image, gradients);
    eddyParameters.m_Misc.m_DoAddEddyCurrents = true;
    eddyParameters.m_SignalGen.m_EddyStrength = strength;
    eddyParameters.m_SignalGen.m_Tau = 1e9;
    auto eddy = Acquire(image, eddyParameters);

    auto fmapParameters = MakeParameters(image, gradients);
    fmapParameters.m_Misc.m_DoAddDistortions = true;
    auto fmap = FloatImageType::New();
    fmap->SetRegions(image->GetLargestPossibleRegion());
    fmap->SetSpacing(image->GetSpacing());
    fmap->Allocate();
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
        fmap->SetPixel({{x, y, 0}}, 42576000.0 * strength / 1000 * (x - nx/2) * image->GetSpacing()[0] / 1000);
    fmapParameters.m_SignalGen.m_FrequencyMap = fmap;
    auto distorted = Acquire(image, fmapParameters);

    double maxError = 0, maxDistortion = 0;
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
      {
        maxError = std::max(maxError, std::abs(Value(eddy, x, y, 0, 1) - Value(distorted, x, y, 0, 1)));
        maxDistortion = std::max(maxDistortion, std::abs(Value(eddy, x, y, 0, 1) - Complex(f(x, y), 0)));
        // the b = 0 volume has no eddy currents
        CPPUNIT_ASSERT_DOUBLES_EQUAL(f(x, y), Value(eddy, x, y, 0, 0).real(), 0.05);
      }
    MITK_INFO << "Eddy currents: largest deviation from the field map " << maxError << ", distortion up to " << maxDistortion;
    CPPUNIT_ASSERT(maxDistortion > 100);
    CPPUNIT_ASSERT_MESSAGE("eddy currents = equivalent field map", maxError < 0.05);
  }

  /** Relative RMS difference of the complex images of all volumes of coil 0. */
  static double RelativeDifference(const FilterType* test, const FilterType* reference, int nx, int ny, unsigned int numVolumes)
  {
    double difference = 0, norm = 0;
    for (unsigned int g = 0; g < numVolumes; ++g)
      for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
        {
          difference += std::norm(Value(test, x, y, 0, g) - Value(reference, x, y, 0, g));
          norm += std::norm(Value(reference, x, y, 0, g));
        }
    return std::sqrt(difference / norm);
  }

  /**
   * The time-segmented phase (SignalGenerationParameters::m_KspaceTimeSegments) approximates the exact
   * sum of a field map and eddy currents: with the automatic segment count (1.25 phase spread + 8) the
   * complex image differs by less than 0.1 % (relative RMS), with too few segments by far more. Single
   * shot EPI in both phase directions with ghosts and aliasing, and conventional spin echo.
   */
  void TimeSegmentsApproximateExactSum()
  {
    const int nx = 40, ny = 40;
    const Pattern f = [](int x, int y)
    {
      const double dx = (x - 20) / 15.0, dy = (y - 21) / 17.0;
      return dx * dx + dy * dy < 1 ? 400.0 + 6 * x - 4 * y + 150 * std::cos(x / 3.0) : 0.0;
    };
    mitk::SignalGenerationParameters::GradientListType gradients(2);
    gradients[1][0] = 0.6;
    gradients[1][1] = 0.8;
    auto image = MakeImage(nx, ny, 1, 2, f);
    // a field map of about 150 Hz range (6 cycles over the EPI readout of 40 lines) and strong eddy currents
    auto fmap = FloatImageType::New();
    fmap->SetRegions(image->GetLargestPossibleRegion());
    fmap->SetSpacing(image->GetSpacing());
    fmap->Allocate();
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
        fmap->SetPixel({{x, y, 0}}, 120 * std::exp(-((x - 26) * (x - 26) + (y - 14) * (y - 14)) / 60.0) - 30 + 0.8 * y);

    struct Case { std::string name; mitk::SignalGenerationParameters::AcquisitionType type; bool reverse; bool ghostsAndAliasing; };
    for (const auto& c : { Case{"EPI", mitk::SignalGenerationParameters::SingleShotEpi, false, false},
                           Case{"EPI reversed, ghosts, aliasing", mitk::SignalGenerationParameters::SingleShotEpi, true, true},
                           Case{"spin echo", mitk::SignalGenerationParameters::ConventionalSpinEcho, false, false} })
    {
      auto parameters = MakeParameters(image, gradients);
      parameters.m_SignalGen.m_AcquisitionType = c.type;
      parameters.m_SignalGen.m_ReversePhase = c.reverse;
      parameters.m_Misc.m_DoAddDistortions = true;
      parameters.m_SignalGen.m_FrequencyMap = fmap;
      parameters.m_Misc.m_DoAddEddyCurrents = true;
      parameters.m_SignalGen.m_EddyStrength = 0.05;
      parameters.m_Misc.m_DoAddGhosts = c.ghostsAndAliasing;
      parameters.m_SignalGen.m_KspaceLineOffset = 0.2;
      parameters.m_Misc.m_DoAddAliasing = c.ghostsAndAliasing;
      parameters.m_SignalGen.m_CroppingFactor = 0.9;
      auto exact = Acquire(image, parameters);
      const int my = static_cast<int>(exact->GetOutput()->GetLargestPossibleRegion().GetSize(1));

      double distortion = 0, norm = 0;
      for (int y = 0; y < std::min(ny, my); ++y)
        for (int x = 0; x < nx; ++x)
        {
          distortion += std::norm(Value(exact, x, y, 0, 1) - Complex(f(x, y), 0));
          norm += f(x, y) * f(x, y);
        }

      parameters.m_SignalGen.m_KspaceTimeSegments = mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS;
      const double automatic = RelativeDifference(Acquire(image, parameters), exact, nx, my, 2);
      parameters.m_SignalGen.m_KspaceTimeSegments = 2;
      const double few = RelativeDifference(Acquire(image, parameters), exact, nx, my, 2);
      MITK_INFO << "Time segments, " << c.name << ": distortion " << std::sqrt(distortion / norm) << ", error with the automatic count " << automatic << ", with 2 segments " << few;
      if (c.type == mitk::SignalGenerationParameters::SingleShotEpi && !c.ghostsAndAliasing)
        CPPUNIT_ASSERT_MESSAGE(c.name + ": the case distorts the image", std::sqrt(distortion / norm) > 0.3);
      CPPUNIT_ASSERT_MESSAGE(c.name + ": automatic count within 0.1 % of the exact sum", automatic < 1e-3);
      if (c.type == mitk::SignalGenerationParameters::SingleShotEpi)
        CPPUNIT_ASSERT_MESSAGE(c.name + ": 2 segments are far off", few > 1e-2);
    }
  }

  /** A 64x64 slice to acquire with KspaceImageFilter directly: EPI, two compartments in an ellipse. */
  struct PhaseCase
  {
    std::string name;
    double spacing;      ///< mm; with the line time it sets the phase spread
    float tLine;         ///< ms
    float eddyStrength;  ///< mT/m along (0.6, 0.8, 0) at b = 1000; 0 without eddy currents
    Pattern fieldMap;    ///< Hz; empty without distortions
    bool relaxation;     ///< T2 90 and 2000 ms; without, the compartments share one set of spectra
  };

  /** Complex k-space of the slice with the given m_KspaceTimeSegments; numSegments is the count used (0: exact sum). */
  static std::vector<std::complex<float>> AcquireSlice(const PhaseCase& c, int segments, unsigned int& numSegments)
  {
    const int n = 64;
    mitk::FiberfoxParameters parameters;
    auto& signalGen = parameters.m_SignalGen;
    signalGen.m_ImageRegion.SetSize(0, n);
    signalGen.m_ImageRegion.SetSize(1, n);
    signalGen.m_ImageRegion.SetSize(2, 1);
    signalGen.m_CroppedRegion = signalGen.m_ImageRegion;
    signalGen.m_ImageSpacing.Fill(c.spacing);
    signalGen.m_tLine = c.tLine;
    signalGen.m_tEcho = 2 * n * c.tLine;
    signalGen.m_DoSimulateRelaxation = c.relaxation;
    signalGen.m_SignalScale = 1;
    signalGen.m_KspaceTimeSegments = segments;
    signalGen.m_EddyStrength = c.eddyStrength;
    parameters.m_Misc.m_DoAddEddyCurrents = c.eddyStrength > 0;
    parameters.m_Misc.m_DoAddDistortions = static_cast<bool>(c.fieldMap);
    if (c.fieldMap)
    {
      auto fmap = FloatImageType::New();
      fmap->SetRegions(signalGen.m_ImageRegion);
      fmap->SetSpacing(signalGen.m_ImageSpacing);
      fmap->Allocate();
      for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
          fmap->SetPixel({{x, y, 0}}, static_cast<float>(c.fieldMap(x, y)));
      signalGen.m_FrequencyMap = fmap;
    }

    std::vector<itk::KspaceImageFilter<float>::InputImagePointerType> compartments;
    for (int i = 0; i < 2; ++i)
    {
      auto compartment = itk::KspaceImageFilter<float>::InputImageType::New();
      itk::ImageRegion<2> region;
      region.SetSize(0, n);
      region.SetSize(1, n);
      compartment->SetRegions(region);
      compartment->Allocate();
      for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
        {
          const double dx = (x - 0.5 * n) / (0.4 * n), dy = (y - 0.52 * n) / (0.45 * n);
          const double r2 = dx * dx + dy * dy;
          double value = 0;
          if (r2 < 1)
            value = i == 0 ? 300 + 200 * std::cos(0.4 * x) * std::sin(0.3 * y) : (r2 < 0.15 ? 800 : 50);
          compartment->SetPixel({{x, y}}, static_cast<float>(value));
        }
      compartments.push_back(compartment);
    }

    auto filter = KspaceSegmentsProbe::New();
    filter->SetCompartmentImages(compartments);
    filter->SetT2({90, 2000});
    filter->SetT1({900, 4000});
    filter->SetParameters(&parameters);
    filter->SetZ(0);
    filter->SetZidx(0);
    itk::Vector<float, 3> coil;
    coil.Fill(0);
    filter->SetCoilPosition(coil);
    filter->SetMotion(mitk::FiberfoxMotion());
    itk::Vector<double, 3> gradient;
    gradient[0] = 0.6;
    gradient[1] = 0.8;
    gradient[2] = 0;
    filter->SetDiffusionGradientDirection(gradient);
    filter->SetRandSeed(0);
    filter->Update();
    numSegments = filter->NumSegments();
    const auto* data = filter->GetOutput()->GetBufferPointer();
    return std::vector<std::complex<float>>(data, data + filter->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels());
  }

  static double RelativeDifference(const std::vector<std::complex<float>>& test, const std::vector<std::complex<float>>& reference)
  {
    double difference = 0, norm = 0;
    for (std::size_t i = 0; i < test.size(); ++i)
    {
      difference += std::norm(Complex(test[i]) - Complex(reference[i]));
      norm += std::norm(Complex(reference[i]));
    }
    return std::sqrt(difference / norm);
  }

  /**
   * The default, the automatic count of time segments, stays within 0.1 % (relative RMS of the complex
   * k-space) of the exact sum also where the phase is strong: eddy currents, whose phase changes about
   * five times faster at the start of the readout than at the echo (an image of 128 mm at 4 mm and a
   * line time of 2 ms has the phase spread of 256 mm at 2 mm and 1 ms), a field map together with eddy
   * currents (more than 512 bins of the field histogram), a field map whose values cluster (a linear
   * field, one value per column, plus a bump), and a field map alone with a phase spread above 128
   * cycles, also with a fixed count of 200 segments. Where the time segments are used, the case must
   * test them (the slice must not have taken the exact sum).
   */
  void TimeSegmentsHoldForStrongPhase()
  {
    const auto bump = [](double amplitude, double linear)
    {
      return Pattern([amplitude, linear](int x, int y)
      {
        const double dx = x - 0.62 * 64, dy = y - 0.3 * 64, sigma = 0.12 * 64;
        return amplitude * std::exp(-(dx * dx + dy * dy) / (2 * sigma * sigma)) - 0.2 * amplitude + linear * (x - 32) / 64.0;
      });
    };
    struct Expectation { PhaseCase phase; bool segmented; int fixedSegments; };
    const std::vector<Expectation> cases = {
      {{"eddy currents 0.02 mT/m", 4, 2, 0.02f, Pattern(), true}, true, 0},
      {{"eddy currents 0.05 mT/m", 4, 2, 0.05f, Pattern(), true}, true, 0},
      {{"field map 200 Hz, eddy currents 0.1 mT/m", 2, 1, 0.1f, bump(200, 0), true}, true, 0},
      {{"field map 400 Hz, eddy currents 0.05 mT/m", 4, 2, 0.05f, bump(400, 0), true}, false, 0},
      {{"clustered field map", 2, 1, 0, bump(500, 1500), false}, true, 0},
      {{"field map, phase spread 140 cycles", 2, 2, 0, bump(400, 1300), true}, false, 200} };
    for (const auto& c : cases)
    {
      unsigned int numSegments = 0;
      const auto exact = AcquireSlice(c.phase, 0, numSegments);
      const auto automatic = AcquireSlice(c.phase, mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS, numSegments);
      const double error = RelativeDifference(automatic, exact);
      MITK_INFO << "Time segments, " << c.phase.name << ": " << numSegments << " segments, error " << error;
      CPPUNIT_ASSERT_MESSAGE(c.phase.name + ": automatic count within 0.1 % of the exact sum", error < 1e-3);
      if (c.segmented)
        CPPUNIT_ASSERT_MESSAGE(c.phase.name + ": time segments used", numSegments > 0);
      if (c.fixedSegments > 0)
      {
        const auto fixed = AcquireSlice(c.phase, c.fixedSegments, numSegments);
        const double fixedError = RelativeDifference(fixed, exact);
        MITK_INFO << "Time segments, " << c.phase.name << ": " << numSegments << " fixed segments, error " << fixedError;
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(c.fixedSegments), numSegments);
        CPPUNIT_ASSERT_MESSAGE(c.phase.name + ": fixed count within 0.1 % of the exact sum", fixedError < 1e-3);
      }
    }
  }

  /** Complex Gaussian noise of variance sigma^2 in both image components: Rayleigh magnitude in an empty image. */
  void NoiseFullFourier()
  {
    const double sigma2 = 100;
    auto image = MakeImage(32, 32, 4, 2, [](int, int) { return 0.0; });
    auto parameters = MakeParameters(image);
    parameters.m_Misc.m_DoAddNoise = true;
    parameters.m_SignalGen.m_NoiseVariance = sigma2;
    auto filter = Acquire(image, parameters);

    double varReal, varImag, meanMagnitude;
    NoiseStatistics(filter, varReal, varImag, meanMagnitude);
    MITK_INFO << "Noise: variance " << varReal << " (real), " << varImag << " (imaginary), mean magnitude " << meanMagnitude;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma2, varReal, 0.06 * sigma2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma2, varImag, 0.06 * sigma2);
    const double rayleighMean = std::sqrt(sigma2 * itk::Math::pi / 2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(rayleighMean, meanMagnitude, 0.03 * rayleighMean);
  }

  /**
   * The noise variance is that of a fully sampled image: partial Fourier fills the missing lines with
   * conjugates, which moves noise power from the imaginary to the real part but keeps the total.
   */
  void NoisePartialFourier()
  {
    const double sigma2 = 100;
    auto image = MakeImage(32, 32, 4, 2, [](int, int) { return 0.0; });
    auto parameters = MakeParameters(image);
    parameters.m_Misc.m_DoAddNoise = true;
    parameters.m_SignalGen.m_NoiseVariance = sigma2;
    parameters.m_SignalGen.m_PartialFourier = 0.75;
    auto filter = Acquire(image, parameters);

    double varReal, varImag, meanMagnitude;
    NoiseStatistics(filter, varReal, varImag, meanMagnitude);
    MITK_INFO << "Noise with partial Fourier 0.75: variance " << varReal << " (real), " << varImag << " (imaginary)";
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2 * sigma2, varReal + varImag, 0.06 * 2 * sigma2);
    CPPUNIT_ASSERT_MESSAGE("conjugate filling: more noise in the real part", varReal > 1.3 * varImag);
  }

  /** Root of the mean squared magnitude of n coils with noise variance sigma^2: sigma/sqrt(n) chi_2n distributed. */
  void NoiseMultiCoil()
  {
    const double sigma2 = 100;
    const unsigned int coils = 4;
    auto image = MakeImage(32, 32, 2, 2, [](int, int) { return 0.0; });
    auto parameters = MakeParameters(image);
    parameters.m_Misc.m_DoAddNoise = true;
    parameters.m_SignalGen.m_NoiseVariance = sigma2;
    parameters.m_SignalGen.m_NumberOfCoils = coils;
    auto filter = Acquire(image, parameters);

    auto output = filter->GetOutput();
    const std::size_t n = output->GetLargestPossibleRegion().GetNumberOfPixels() * output->GetNumberOfComponentsPerPixel();
    double sum = 0;
    for (std::size_t i = 0; i < n; ++i)
      sum += output->GetBufferPointer()[i];
    const double k = 2.0 * coils;
    const double expected = std::sqrt(sigma2 / coils) * std::sqrt(2.0) * std::tgamma((k + 1) / 2) / std::tgamma(k / 2);
    MITK_INFO << coils << " coils: mean magnitude " << sum / n << " (expected " << expected << ")";
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(coils), filter->GetOutputImagesReal().size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, sum / n, 0.03 * expected);
  }

  /** With a noise model (Rician, chi-square) the model adds the noise of the run after the k-space
      step, so the k-space step adds no complex Gaussian noise of its own, whatever its variance. */
  void NoKspaceNoiseWithNoiseModel()
  {
    auto image = MakeImage(32, 32, 2, 2, [](int, int) { return 0.0; });
    auto parameters = MakeParameters(image);
    parameters.m_Misc.m_DoAddNoise = true;
    parameters.m_SignalGen.m_NoiseVariance = 100;
    parameters.m_NoiseModel = std::make_shared< mitk::RicianNoiseModel<> >();
    parameters.m_NoiseModel->SetNoiseVariance(100);
    auto filter = Acquire(image, parameters);

    double varReal, varImag, meanMagnitude;
    NoiseStatistics(filter, varReal, varImag, meanMagnitude);
    CPPUNIT_ASSERT_EQUAL(0.0, meanMagnitude);
  }

  /**
   * Every readout scheme visits each k-space sample once, starts at the line read first (the highest
   * line in normal, line 0 in reversed phase direction), and has its echo in the middle of the
   * readout (single shot EPI) or of each line (spin echo).
   */
  void ReadoutSchemes()
  {
    for (auto type : { mitk::SignalGenerationParameters::SingleShotEpi, mitk::SignalGenerationParameters::ConventionalSpinEcho, mitk::SignalGenerationParameters::FastSpinEcho })
      for (bool reverse : { false, true })
        for (auto size : { std::array<int, 2>{16, 16}, std::array<int, 2>{15, 13} })
        {
          mitk::FiberfoxParameters parameters;
          parameters.m_SignalGen.m_AcquisitionType = type;
          parameters.m_SignalGen.m_ReversePhase = reverse;
          parameters.m_SignalGen.m_CroppedRegion.SetSize(0, size[0]);
          parameters.m_SignalGen.m_CroppedRegion.SetSize(1, size[1]);
          parameters.m_SignalGen.m_tLine = 0.8f;
          parameters.m_SignalGen.m_tEcho = 100;
          auto scheme = itk::KspaceImageFilter<float>::CreateReadoutScheme(&parameters);
          const int nx = size[0], ny = size[1], n = nx * ny;
          const std::string label = std::to_string(type) + (reverse ? " reversed " : " normal ") + std::to_string(nx) + "x" + std::to_string(ny);

          std::set<std::pair<int, int>> visited;
          for (int tick = 0; tick < n; ++tick)
          {
            const auto index = scheme->GetActualKspaceIndex(tick);
            CPPUNIT_ASSERT_MESSAGE(label + ": index inside k-space", index[0] >= 0 && index[0] < nx && index[1] >= 0 && index[1] < ny);
            visited.insert({static_cast<int>(index[0]), static_cast<int>(index[1])});
            // fast spin echo: TE is the echo of the central line of the echo train
            if (type != mitk::SignalGenerationParameters::FastSpinEcho)
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(label + ": time from RF", parameters.m_SignalGen.m_tEcho + scheme->GetTimeFromMaxEcho(tick), scheme->GetTimeFromRf(tick), 1e-4);
          }
          CPPUNIT_ASSERT_EQUAL_MESSAGE(label + ": every sample once", static_cast<std::size_t>(n), visited.size());
          CPPUNIT_ASSERT_EQUAL_MESSAGE(label + ": first line", reverse ? 0 : ny - 1, static_cast<int>(scheme->GetActualKspaceIndex(0)[1]));

          const double dt = parameters.m_SignalGen.m_tLine / nx;
          if (type == mitk::SignalGenerationParameters::SingleShotEpi)
          {
            // one echo in the middle of the whole readout, ny lines of tLine
            CPPUNIT_ASSERT_DOUBLES_EQUAL(-n * dt / 2 + dt / 2, scheme->GetTimeFromMaxEcho(0), 1e-4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(n * dt / 2 - dt / 2, scheme->GetTimeFromMaxEcho(n - 1), 1e-4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(ny * parameters.m_SignalGen.m_tLine, scheme->GetMinimumEchoTime(), 1e-4);
            // zig-zag: the second line runs backwards
            CPPUNIT_ASSERT_EQUAL(scheme->GetActualKspaceIndex(nx - 1)[0], scheme->GetActualKspaceIndex(nx)[0]);
          }
          else
          {
            // one echo in the middle of every line
            for (int line : { 0, ny / 2, ny - 1 })
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(-nx * dt / 2 + dt / 2, scheme->GetTimeFromMaxEcho(line * nx), 1e-4);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(nx * dt / 2 - dt / 2, scheme->GetTimeFromMaxEcho(line * nx + nx - 1), 1e-4);
              CPPUNIT_ASSERT_EQUAL(0, static_cast<int>(scheme->GetActualKspaceIndex(line * nx)[0]));
            }
            const double minimum = type == mitk::SignalGenerationParameters::FastSpinEcho ? parameters.m_SignalGen.m_EchoTrainLength * parameters.m_SignalGen.m_tLine : parameters.m_SignalGen.m_tLine;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(minimum, scheme->GetMinimumEchoTime(), 1e-4);
          }
        }
  }

  /** Field map = sum over the axes of (index * gradient + offset) plus Gaussians h exp(-d^2/(2 var)) at world positions. */
  void FieldmapGenerator()
  {
    using GeneratorType = itk::FieldmapGeneratorFilter<FloatImageType>;
    auto generator = GeneratorType::New();
    itk::ImageRegion<3> region;
    region.SetSize(0, 10);
    region.SetSize(1, 8);
    region.SetSize(2, 6);
    generator->SetImageRegion(region);
    itk::Vector<double> spacing;
    spacing[0] = 2; spacing[1] = 1.5; spacing[2] = 3;
    generator->SetSpacing(spacing);
    mitk::Point3D origin;
    origin[0] = -5; origin[1] = 7; origin[2] = 1;
    generator->SetOrigin(origin);
    GeneratorType::MatrixType direction;
    direction.SetIdentity();
    direction[0][0] = direction[1][1] = std::cos(0.3);
    direction[0][1] = -std::sin(0.3);
    direction[1][0] = std::sin(0.3);
    generator->SetDirectionMatrix(direction);
    vnl_vector_fixed<double, 3> gradient(0.5, -1.0, 2.0), offset(1.0, 2.0, -0.5);
    generator->SetGradient(gradient);
    generator->SetOffset(offset);
    mitk::Point3D p1, p2;
    p1[0] = 3; p1[1] = 12; p1[2] = 8;
    p2[0] = -1; p2[1] = 9; p2[2] = 4;
    generator->SetWorldPositions({ p1, p2 });
    generator->SetHeights({ 50, -20 });
    generator->SetVariances({ 9, 4 });
    generator->Update();
    auto fmap = generator->GetOutput();

    CPPUNIT_ASSERT(fmap->GetLargestPossibleRegion() == region);
    for (const auto& index : { itk::Index<3>{{0, 0, 0}}, itk::Index<3>{{3, 4, 2}}, itk::Index<3>{{9, 7, 5}}, itk::Index<3>{{4, 2, 1}} })
    {
      double expected = 0;
      for (unsigned int i = 0; i < 3; ++i)
        expected += index[i] * gradient[i] + offset[i];
      itk::Point<double, 3> world;
      for (unsigned int i = 0; i < 3; ++i)
      {
        world[i] = origin[i];
        for (unsigned int j = 0; j < 3; ++j)
          world[i] += direction[i][j] * spacing[j] * index[j];
      }
      expected += 50 * std::exp(-world.SquaredEuclideanDistanceTo(p1) / (2 * 9)) - 20 * std::exp(-world.SquaredEuclideanDistanceTo(p2) / (2 * 4));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, fmap->GetPixel(index), 1e-4 * std::max(1.0, std::fabs(expected)));
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxKspace)
