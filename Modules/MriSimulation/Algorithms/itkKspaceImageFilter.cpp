/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkKspaceImageFilter_txx
#define __itkKspaceImageFilter_txx
//#endif

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <cmath>

#include "itkKspaceImageFilter.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <mitkSingleShotEpi.h>
#include <mitkConventionalSpinEcho.h>
#include <mitkFastSpinEcho.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <itkImageFileWriter.h>
#include <limits>
#include <boost/lexical_cast.hpp>
#include <array>
#include <map>
#include <numeric>
#include <functional>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>

namespace itk {

  template< class ScalarType >
  KspaceImageFilter< ScalarType >::KspaceImageFilter()
    : m_Z(0)
    , m_RandSeed(-1)
    , m_SpikesPerSlice(0)
    , m_IsBaseline(true)
    , m_StoreTimings(false)
    , m_NumSegments(0)
    , m_TimeSegmentsCapped(false)
  {
    m_DiffusionGradientDirection.Fill(0.0);
    m_CoilPosition.Fill(0.0);

    
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >
  ::BeforeThreadedGenerateData()
  {
    m_Spike = std::complex<ScalarType>(0,0);
    m_SpikeTick = std::numeric_limits<long>::max();
    m_SpikeLog = "";

    kxMax = m_Parameters->m_SignalGen.m_CroppedRegion.GetSize(0);
    kyMax = m_Parameters->m_SignalGen.m_CroppedRegion.GetSize(1);
    xMax = m_CompartmentImages.at(0)->GetLargestPossibleRegion().GetSize(0); // scanner coverage in x-direction
    yMax = m_CompartmentImages.at(0)->GetLargestPossibleRegion().GetSize(1); // scanner coverage in y-direction
    yMaxFov = yMax;
    if (m_Parameters->m_Misc.m_DoAddAliasing)
    {
      // actual FOV in y-direction (in x-direction FOV=xMax)
      yMaxFov = static_cast<int>(yMaxFov * m_Parameters->m_SignalGen.m_CroppingFactor);
    }
    yMaxFov_half = (yMaxFov-1)/2;
    numPix = kxMax*kyMax;

    float ringing_factor = static_cast<float>(m_Parameters->m_SignalGen.m_ZeroRinging)/100.0;
    ringing_lines_x = static_cast<int>(ceil(kxMax/2 * ringing_factor));
    ringing_lines_y = static_cast<int>(ceil(kyMax/2 * ringing_factor));

    // Adjust noise variance since it is the intended variance in physical space and not in k-space.
    // It is the per-component image variance of a fully sampled acquisition: with partial Fourier the
    // conjugate filling already keeps the total noise power (moving it from the imaginary to the real
    // part), so the variance is not scaled with the acquired fraction.
    float noiseVar = m_Parameters->m_SignalGen.m_NoiseVariance/(kyMax*kxMax);

    m_RandGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    if (m_RandSeed>=0)  // always generate the same random numbers?
      m_RandGen->SetSeed(m_RandSeed);
    else
      m_RandGen->SetSeed();

    typename OutputImageType::Pointer outputImage = OutputImageType::New();
    itk::ImageRegion<2> region;
    region.SetSize(0, m_Parameters->m_SignalGen.m_CroppedRegion.GetSize(0));
    region.SetSize(1, m_Parameters->m_SignalGen.m_CroppedRegion.GetSize(1));
    outputImage->SetLargestPossibleRegion( region );
    outputImage->SetBufferedRegion( region );
    outputImage->SetRequestedRegion( region );
    outputImage->Allocate();
    std::complex<ScalarType> zero = std::complex<ScalarType>(0, 0);
    outputImage->FillBuffer(zero);
    // A noise model (Rician, noncentral chi) adds the noise of the run itself, after the k-space step.
    if (m_Parameters->m_SignalGen.m_NoiseVariance>0 && m_Parameters->m_Misc.m_DoAddNoise && m_Parameters->m_NoiseModel==nullptr)
    {
      ImageRegionIterator< OutputImageType > oit(outputImage, outputImage->GetLargestPossibleRegion());
      while( !oit.IsAtEnd() )
      {
        oit.Set(std::complex<ScalarType>(m_RandGen->GetNormalVariate(0, noiseVar), m_RandGen->GetNormalVariate(0, noiseVar)));
        ++oit;
      }
    }

    m_KSpaceImage = InputImageType::New();
    m_KSpaceImage->SetLargestPossibleRegion( region );
    m_KSpaceImage->SetBufferedRegion( region );
    m_KSpaceImage->SetRequestedRegion( region );
    m_KSpaceImage->Allocate();
    m_KSpaceImage->FillBuffer(0.0);

    if (m_StoreTimings)
    {
        m_TickImage = InputImageType::New();
        m_TickImage->SetLargestPossibleRegion( region );
        m_TickImage->SetBufferedRegion( region );
        m_TickImage->SetRequestedRegion( region );
        m_TickImage->Allocate();
        m_TickImage->FillBuffer(-1.0);

        m_RfImage = InputImageType::New();
        m_RfImage->SetLargestPossibleRegion( region );
        m_RfImage->SetBufferedRegion( region );
        m_RfImage->SetRequestedRegion( region );
        m_RfImage->Allocate();
        m_RfImage->FillBuffer(-1.0);
    }
    else
    {
        m_TickImage = nullptr;
        m_RfImage = nullptr;
    }

    // Gyromagnetic ratio in Hz/T: the eddy phase phi enters exp(i 2 pi (k x + phi)) in cycles, like the
    // field map term (Hz times s), so the field G*r in T times gamma-bar gives Hz; the angular
    // 2 pi * 42.576 MHz/T would make the eddy field 2 pi times stronger than the documented mT/m.
    m_Gamma = 42576000;
    if ( m_Parameters->m_SignalGen.m_EddyStrength>0 && m_DiffusionGradientDirection.GetNorm()>0.001)
    {
      // a separate member, so that m_DiffusionGradientDirection stays unscaled when Update() runs again
      m_EddyFieldGradient = m_DiffusionGradientDirection *  m_Parameters->m_SignalGen.m_EddyStrength/1000 *  m_Gamma;
      m_IsBaseline = false;
    }
    else {
      m_IsBaseline = true;
    }

    this->SetNthOutput(0, outputImage);

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        m_Transform[i][j] = m_Parameters->m_SignalGen.m_ImageDirection[i][j] * m_Parameters->m_SignalGen.m_ImageSpacing[j]/1000;

    float a = m_Parameters->m_SignalGen.m_ImageRegion.GetSize(0)*m_Parameters->m_SignalGen.m_ImageSpacing[0];
    float b = m_Parameters->m_SignalGen.m_ImageRegion.GetSize(1)*m_Parameters->m_SignalGen.m_ImageSpacing[1];
    float diagonal = sqrt(a*a+b*b)/1000;   // image diagonal in m

    switch (m_Parameters->m_SignalGen.m_CoilSensitivityProfile)
    {
      case mitk::SignalGenerationParameters::COIL_CONSTANT:
      {
        m_CoilSensitivityFactor = 1;                    // same signal everywhere
        break;
      }
      case mitk::SignalGenerationParameters::COIL_LINEAR:
      {
        m_CoilSensitivityFactor = -1/diagonal;          // about 50% of the signal in the image center remaining
        break;
      }
      case mitk::SignalGenerationParameters::COIL_EXPONENTIAL:
      {
        m_CoilSensitivityFactor = -log(0.1)/diagonal;   // about 32% of the signal in the image center remaining
        break;
      }
    }

    m_ReadoutScheme = CreateReadoutScheme(m_Parameters);

    m_MovedFmap = nullptr;
    if (m_Parameters->m_Misc.m_DoAddDistortions && m_Parameters->m_SignalGen.m_FrequencyMap.IsNotNull() && m_Parameters->m_SignalGen.m_DoAddMotion)
    {
      // the frequency map moves with the head: the moved map at x is the map at the unmoved position of x
      itk::LinearInterpolateImageFunction< itk::Image< float, 3 >, float >::Pointer   fmapInterpolator;
      fmapInterpolator = itk::LinearInterpolateImageFunction< itk::Image< float, 3 >, float >::New();
      fmapInterpolator->SetInputImage(m_Parameters->m_SignalGen.m_FrequencyMap);

      m_MovedFmap = itk::Image< ScalarType, 2 >::New();
      m_MovedFmap->SetLargestPossibleRegion( m_CompartmentImages.at(0)->GetLargestPossibleRegion() );
      m_MovedFmap->SetBufferedRegion( m_CompartmentImages.at(0)->GetLargestPossibleRegion() );
      m_MovedFmap->SetRequestedRegion( m_CompartmentImages.at(0)->GetLargestPossibleRegion() );
      m_MovedFmap->Allocate();
      m_MovedFmap->FillBuffer(0);

      ImageRegionIterator< InputImageType > it(m_MovedFmap, m_MovedFmap->GetLargestPossibleRegion() );
      while( !it.IsAtEnd() )
      {
        itk::Image<float, 3>::IndexType index;
        index[0] = it.GetIndex()[0];
        index[1] = it.GetIndex()[1];
        index[2] = m_Zidx;

        itk::Point<float, 3> point3D;
        m_Parameters->m_SignalGen.m_FrequencyMap->TransformIndexToPhysicalPoint(index, point3D);
        it.Set(mitk::imv::GetImageValue<float>(m_Motion.Inverse(point3D), true, fmapInterpolator));
        ++it;
      }
    }

    // calculate T1 relaxation (independent of actual readout)
    m_T1Relax.clear();
    if ( m_Parameters->m_SignalGen.m_DoSimulateRelaxation)
    {
      // the inversion-recovery signal below holds for TI < TR only
      m_Parameters->m_SignalGen.CheckInversionTime();
      for (unsigned int i=0; i<m_CompartmentImages.size(); i++)
      {
        // account for T1 relaxation (how much magnetization is available since last excitation?)
        float relaxation = (1.0-std::exp(-m_Parameters->m_SignalGen.m_tRep/m_T1[i]));

        // account for inversion pulse and TI
        const float tInv = m_Parameters->m_SignalGen.m_tInv;
        if (tInv > 0)
        {
          // inversion recovery with a 90 degree excitation every TR: 1 - 2 exp(-TI/T1) + exp(-TR/T1)
          relaxation = 1.0 - 2.0*std::exp(-tInv/m_T1[i]) + std::exp(-m_Parameters->m_SignalGen.m_tRep/m_T1[i]);
        }

        m_T1Relax.push_back(relaxation);
      }
    }

    m_CompartmentSpectra.clear();
    m_NumSegments = 0;
    m_TimeSegmentsCapped = false;
    const bool eddyCurrents = m_Parameters->m_Misc.m_DoAddEddyCurrents && m_Parameters->m_SignalGen.m_EddyStrength>0 && !m_IsBaseline;
    const bool distortions = m_Parameters->m_Misc.m_DoAddDistortions && (m_MovedFmap.IsNotNull() || m_Parameters->m_SignalGen.m_FrequencyMap.IsNotNull());
    if (!eddyCurrents && !distortions)
      ComputeCompartmentSpectra();
    else if (m_Parameters->m_SignalGen.m_KspaceTimeSegments != 0)
      ComputeTimeSegmentation(eddyCurrents, distortions);
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >::ComputeDftFactors(std::vector< std::complex<double> >& factorsX, std::vector< std::complex<double> >& factorsY, std::vector< double >& weight) const
  {
    const std::size_t nx = static_cast<std::size_t>(xMax);
    const std::size_t ny = static_cast<std::size_t>(yMax);
    const std::size_t nkx = static_cast<std::size_t>(kxMax);
    const std::size_t nky = static_cast<std::size_t>(kyMax);
    // centred indices as in the sum over pixels in DynamicThreadedGenerateData: (0 -- N) --> (-N/2 -- N/2)
    const double xShift = static_cast<double>(nx/2);
    const double yShift = static_cast<double>(ny/2);
    const double kxShift = static_cast<double>(nkx/2);
    const double kyShift = static_cast<double>(nky/2);
    const bool ghosts = m_Parameters->m_Misc.m_DoAddGhosts;
    const std::size_t parities = ghosts ? 2 : 1;

    // factors along x per line parity (the ghost offset shifts kx of odd and even lines apart)
    factorsX.assign(parities*nkx*nx, std::complex<double>(0, 0));
    for (std::size_t p=0; p<parities; ++p)
      for (std::size_t kx=0; kx<nkx; ++kx)
      {
        double k = static_cast<double>(kx) - kxShift;
        if (ghosts)
          k += (p == 1) ? -m_Parameters->m_SignalGen.m_KspaceLineOffset : m_Parameters->m_SignalGen.m_KspaceLineOffset;
        k /= xMax;
        for (std::size_t x=0; x<nx; ++x)
          factorsX[(p*nkx + kx)*nx + x] = std::polar(1.0, itk::Math::twopi * k * (static_cast<double>(x) - xShift));
      }

    // factors along y; signal from outside the FOV is wrapped back (aliasing)
    factorsY.assign(nky*ny, std::complex<double>(0, 0));
    for (std::size_t y=0; y<ny; ++y)
    {
      double position = static_cast<double>(y) - yShift;
      if (m_Parameters->m_Misc.m_DoAddAliasing)
      {
        if (position < -yMaxFov_half)
          position += yMaxFov;
        else if (position > yMaxFov_half)
          position -= yMaxFov;
      }
      for (std::size_t ky=0; ky<nky; ++ky)
        factorsY[ky*ny + y] = std::polar(1.0, itk::Math::twopi * (static_cast<double>(ky) - kyShift) / yMaxFov * position);
    }

    // signal scale and coil sensitivity depend on the position only
    weight.assign(nx*ny, m_Parameters->m_SignalGen.m_SignalScale);
    if (m_Parameters->m_SignalGen.m_CoilSensitivityProfile!=mitk::SignalGenerationParameters::COIL_CONSTANT)
      for (std::size_t y=0; y<ny; ++y)
        for (std::size_t x=0; x<nx; ++x)
        {
          VectorType pos;
          pos[0] = static_cast<float>(x) - static_cast<float>(xShift);
          pos[1] = static_cast<float>(y) - static_cast<float>(yShift);
          pos[2] = m_Z;
          pos = m_Transform*pos;
          weight[y*nx + x] *= CoilSensitivity(pos);
        }
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >::ComputeCompartmentSpectra()
  {
    const std::size_t nx = static_cast<std::size_t>(xMax);
    const std::size_t ny = static_cast<std::size_t>(yMax);
    const std::size_t nkx = static_cast<std::size_t>(kxMax);
    const std::size_t nky = static_cast<std::size_t>(kyMax);
    const bool ghosts = m_Parameters->m_Misc.m_DoAddGhosts;
    const std::size_t parities = ghosts ? 2 : 1;
    std::vector< std::complex<double> > factorsX, factorsY;
    std::vector< double > weight;
    ComputeDftFactors(factorsX, factorsY, weight);

    std::vector< std::complex<double> > rows(parities*ny*nkx);
    for (const auto& compartment : m_CompartmentImages)
    {
      itkAssertOrThrowMacro(compartment->GetBufferedRegion() == compartment->GetLargestPossibleRegion(), "The compartment slices must be buffered whole.");
      const ScalarType* signal = compartment->GetBufferPointer();

      // separable: first along x for every row and parity, then along y
      for (std::size_t p=0; p<parities; ++p)
        for (std::size_t y=0; y<ny; ++y)
          for (std::size_t kx=0; kx<nkx; ++kx)
          {
            std::complex<double> sum(0, 0);
            const std::complex<double>* factors = &factorsX[(p*nkx + kx)*nx];
            for (std::size_t x=0; x<nx; ++x)
              sum += (signal[y*nx + x] * weight[y*nx + x]) * factors[x];
            rows[(p*ny + y)*nkx + kx] = sum;
          }

      std::vector< std::complex<double> > spectrum(nkx*nky);
      for (std::size_t ky=0; ky<nky; ++ky)
      {
        const std::size_t p = ghosts ? ky%2 : 0;
        for (std::size_t kx=0; kx<nkx; ++kx)
        {
          std::complex<double> sum(0, 0);
          for (std::size_t y=0; y<ny; ++y)
            sum += factorsY[ky*ny + y] * rows[(p*ny + y)*nkx + kx];
          spectrum[ky*nkx + kx] = sum;
        }
      }
      m_CompartmentSpectra.push_back(std::move(spectrum));
    }
  }

  /*
   * Time segmentation with least-squares interpolators, as in Sutton, Noll and Fessler, "Fast,
   * iterative image reconstruction for MRI in the presence of field inhomogeneities", IEEE TMI
   * 22(2):178-188, 2003 (section II-C), here in the forward direction. The phase of pixel r at the
   * readout time of tick j is phi(r, j) = u(r) . v(j) with u = (field map in Hz, eddy-current field
   * G . pos(r) in Hz) and v = (t in s, eddy factor exp(-tRead/tau) t in s), as in the exact sum. With L
   * segment times t_l,
   *   exp(i 2 pi phi(r, j)) ~ sum_l a_l(j) exp(i 2 pi phi(r, t_l)),
   * so s(k_j) ~ sum_l a_l(j) S_l(k_j), where S_l is the spectrum of the signal times the phase factor
   * at t_l: L separable DFTs per compartment instead of a sum over all pixels per sample. a(j) is the
   * least-squares fit over the pixels' values of u, binned to a histogram weighted by signal energy:
   * minimise sum_b w_b |exp(i 2 pi u_b . v(j)) - sum_l a_l exp(i 2 pi u_b . v_l)|^2.
   *
   * The segment times are spaced evenly in accumulated phase, the largest phase change between two
   * pixels with signal summed along the readout: the eddy factor changes about five times faster at
   * the start of an EPI readout than at the echo, and segments evenly in time missed the 0.1 % bound
   * there (issue #47, review). The error falls steeply once L exceeds the accumulated phase S in
   * cycles (the time-bandwidth product of the fit), hence the automatic count 1.25 S + 8; it grows by
   * a quarter until the energy-weighted RMS interpolation residual of the pixels at sampled readout
   * times is below 1e-4. Where the time segments would take longer than the exact sum on one core
   * (estimated before each fit, with the fits done), or more memory than MAX_TIME_SEGMENT_BYTES, the
   * slice takes the exact sum (GetTimeSegmentsCapped).
   */
  template< class ScalarType >
  void KspaceImageFilter< ScalarType >::ComputeTimeSegmentation(bool eddyCurrents, bool distortions)
  {
    const std::size_t nx = static_cast<std::size_t>(xMax);
    const std::size_t ny = static_cast<std::size_t>(yMax);
    const std::size_t nkx = static_cast<std::size_t>(kxMax);
    const std::size_t nky = static_cast<std::size_t>(kyMax);
    const std::size_t numPixels = nx*ny;
    const std::size_t numTicks = nkx*nky;
    const bool ghosts = m_Parameters->m_Misc.m_DoAddGhosts;
    const std::size_t parities = ghosts ? 2 : 1;

    // v(tick) with the same float arithmetic as the exact sum in DynamicThreadedGenerateData
    std::vector< float > echoTime(numTicks);
    std::vector< std::array<double, 2> > v(numTicks);
    for (std::size_t tick=0; tick<numTicks; ++tick)
    {
      const float t = m_ReadoutScheme->GetTimeFromMaxEcho(static_cast<int>(tick));
      float eddyDecay = 0;
      if (eddyCurrents)
      {
        const float tRead = m_ReadoutScheme->GetTimeFromLastDiffusionGradient(static_cast<int>(tick));
        eddyDecay = std::exp(-tRead/m_Parameters->m_SignalGen.m_Tau ) * t/1000;
      }
      echoTime[tick] = t;
      v[tick] = { distortions ? static_cast<double>(t/1000) : 0.0, static_cast<double>(eddyDecay) };
    }

    // u(r) as in the exact sum
    const float xShift = static_cast<float>(nx/2);
    const float yShift = static_cast<float>(ny/2);
    std::vector< std::array<double, 2> > u(numPixels, {0.0, 0.0});
    for (std::size_t y=0; y<ny; ++y)
      for (std::size_t x=0; x<nx; ++x)
      {
        typename InputImageType::IndexType index;
        index[0] = static_cast<itk::IndexValueType>(x);
        index[1] = static_cast<itk::IndexValueType>(y);
        auto& pixel = u[y*nx + x];
        if (distortions)
        {
          if (m_MovedFmap.IsNotNull())
            pixel[0] = m_MovedFmap->GetPixel(index);
          else
          {
            itk::Image<float, 3>::IndexType index3D; index3D[0] = index[0]; index3D[1] = index[1]; index3D[2] = m_Zidx;
            pixel[0] = m_Parameters->m_SignalGen.m_FrequencyMap->GetPixel(index3D);
          }
        }
        if (eddyCurrents)
        {
          VectorType pos;
          pos[0] = static_cast<float>(x) - xShift; pos[1] = static_cast<float>(y) - yShift; pos[2] = m_Z;
          pos = m_Transform*pos;
          pixel[1] = m_EddyFieldGradient[0]*pos[0]+m_EddyFieldGradient[1]*pos[1]+m_EddyFieldGradient[2]*pos[2];
        }
      }

    std::vector< std::complex<double> > factorsX, factorsY;
    std::vector< double > weight;
    ComputeDftFactors(factorsX, factorsY, weight);

    // signal per spectrum group: every compartment with relaxation (its factor depends on the tick), their sum without
    std::vector< std::vector< double > > groups;
    for (const auto& compartment : m_CompartmentImages)
    {
      itkAssertOrThrowMacro(compartment->GetBufferedRegion() == compartment->GetLargestPossibleRegion(), "The compartment slices must be buffered whole.");
      const ScalarType* signal = compartment->GetBufferPointer();
      if (groups.empty() || m_Parameters->m_SignalGen.m_DoSimulateRelaxation)
        groups.emplace_back(numPixels, 0.0);
      for (std::size_t r=0; r<numPixels; ++r)
        groups.back()[r] += signal[r];
    }

    std::array<double, 2> span = {0.0, 0.0};
    for (unsigned int d=0; d<2; ++d)
    {
      double lo = v[0][d], hi = v[0][d];
      for (const auto& vt : v) { lo = std::min(lo, vt[d]); hi = std::max(hi, vt[d]); }
      span[d] = hi - lo;
    }
    std::vector< double > energy(numPixels, 0.0);
    for (std::size_t r=0; r<numPixels; ++r)
    {
      double magnitude = 0;
      for (const auto& group : groups)
        magnitude += std::fabs(group[r]);
      energy[r] = magnitude*weight[r]*magnitude*weight[r];
    }

    // range of u over the pixels with signal
    std::array<double, 2> lowU = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 2> highU = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (std::size_t r=0; r<numPixels; ++r)
      if (energy[r] > 0)
        for (unsigned int d=0; d<2; ++d)
        {
          lowU[d] = std::min(lowU[d], u[r][d]);
          highU[d] = std::max(highU[d], u[r][d]);
        }
    std::array<double, 2> axisSpread = {0.0, 0.0};
    for (unsigned int d=0; d<2; ++d)
      if (highU[d] > lowU[d])
        axisSpread[d] = (highU[d]-lowU[d])*span[d];

    // Ticks in the order of their readout time (v is a function of it in every readout scheme) and the
    // phase accumulated up to each: the largest phase change (cycles) between two pixels with signal.
    std::vector< std::size_t > order(numTicks);
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::stable_sort(order.begin(), order.end(), [&echoTime](std::size_t a, std::size_t b) { return echoTime[a] < echoTime[b]; });
    std::vector< double > accumulated(numTicks, 0.0);
    for (std::size_t i=1; i<numTicks; ++i)
    {
      double step = 0;
      for (unsigned int d=0; d<2; ++d)
        if (axisSpread[d] > 0)
          step += (highU[d]-lowU[d])*std::fabs(v[order[i]][d] - v[order[i-1]][d]);
      accumulated[i] = accumulated[i-1] + step;
    }
    const double phaseSpread = accumulated.back();
    const int requested = m_Parameters->m_SignalGen.m_KspaceTimeSegments;
    const std::size_t numGroups = groups.size();
    unsigned int numSegments = requested > 0 ? static_cast<unsigned int>(requested) :
                                               phaseSpread > 0 ? static_cast<unsigned int>(std::ceil(1.25*phaseSpread)) + 8 : 1;

    // Estimated times in ns on one core at -O3 (measured, issue #47; only their ratios matter): the
    // exact sum takes about 35 ns per pixel and sample, a complex multiply-add of the vectorized fit and
    // transforms about 0.4 ns, a complex exponential about 15 ns, the SVD about 7 ns per bin and squared
    // segment count. The choice between time segments and the exact sum uses the times on one core, so
    // that it, and with it the result, does not depend on the work units of the job.
    const double exactTime = 35.0*numTicks*numPixels;
    // In parallel when the slice job has more than one work unit (with few volumes the setup took longer
    // than the exact sum in all work units); each index writes its own result, so the result does not
    // depend on the threads. Running the small calls serially made larger slices slower.
    const auto parallelFor = [this](std::size_t count, const std::function<void(std::size_t)>& function)
    {
      if (this->GetNumberOfWorkUnits() > 1)
        this->GetMultiThreader()->ParallelizeArray(0, count, [&function](SizeValueType i) { function(i); }, nullptr);
      else
        for (std::size_t i=0; i<count; ++i)
          function(i);
    };
    // per segment: phase factors of the pixels, transforms along x and y, folding
    const double segmentTime = 15.0*numPixels + numGroups*(0.47*nkx*ny*(parities*nx + nky) + 3.0*numTicks);
    // Where already the transforms of the first count take longer than the exact sum (see below), the
    // slice takes it before the histogram and the fit are set up.
    if (requested <= 0 && numSegments*segmentTime > 0.8*exactTime)
    {
      m_TimeSegmentsCapped = true;
      return;
    }

    // Histogram of u weighted by signal energy. The interpolators fitted at the bins hold for the
    // pixels between them if the bins sample the phase finely: 4 bins per cycle of phase spread along
    // each axis, at least 16 in all for a small spread. Fewer bins break the fit (2 per cycle gave
    // errors above 100 % with a field map and eddy currents), so their number is not limited; it is at
    // most the number of pixels.
    std::map< std::array<long long, 2>, std::array<double, 3> > histogram;
    {
      std::array<double, 2> binWidth = {0.0, 0.0};
      const double totalSpread = axisSpread[0] + axisSpread[1];
      for (unsigned int d=0; d<2; ++d)
        if (axisSpread[d] > 0)
          binWidth[d] = (highU[d]-lowU[d])/std::ceil(std::max(4*axisSpread[d], 16*axisSpread[d]/totalSpread));
      for (std::size_t r=0; r<numPixels; ++r)
      {
        if (energy[r] <= 0)
          continue;
        std::array<long long, 2> key = {0, 0};
        for (unsigned int d=0; d<2; ++d)
          if (binWidth[d] > 0)
            key[d] = static_cast<long long>(std::floor((u[r][d]-lowU[d])/binWidth[d]));
        auto& bin = histogram[key];
        bin[0] += energy[r];
        bin[1] += energy[r]*u[r][0];
        bin[2] += energy[r]*u[r][1];
      }
      // A field whose values cluster (e.g. one value per image column of a linear field) leaves empty
      // bins between the clusters, where the fit is free and grows large; the pixels beside them then
      // got errors up to 11 %. Along a single axis the empty bins in the range are therefore fitted too,
      // at their centre with 1 % of the mean bin weight. With two axes the pixel residual check below
      // guards against this.
      if ((axisSpread[0] > 0) != (axisSpread[1] > 0) && !histogram.empty())
      {
        const unsigned int d = axisSpread[0] > 0 ? 0 : 1;
        double meanEnergy = 0;
        for (const auto& [key, bin] : histogram)
          meanEnergy += bin[0];
        meanEnergy /= histogram.size();
        const long long first = histogram.begin()->first[d];
        const long long last = histogram.rbegin()->first[d];
        for (long long i=first+1; i<last; ++i)
        {
          std::array<long long, 2> key = {0, 0};
          key[d] = i;
          if (histogram.count(key) > 0)
            continue;
          std::array<double, 3> bin = {0.01*meanEnergy, 0, 0};
          bin[1+d] = bin[0]*(lowU[d] + (i+0.5)*binWidth[d]);
          histogram.emplace(key, bin);
        }
      }
    }
    // the distinct v in readout order; ticks with the same v share their interpolators (spin echo lines)
    std::vector< std::size_t > distinct;
    std::vector< std::size_t > weightIndex(numTicks, 0);
    {
      std::map< std::array<double, 2>, std::size_t > seen;
      for (const std::size_t tick : order)
      {
        const auto [found, inserted] = seen.emplace(v[tick], distinct.size());
        if (inserted)
          distinct.push_back(tick);
        weightIndex[tick] = found->second;
      }
    }

    const std::size_t numBins = histogram.size();
    std::vector< std::array<double, 2> > binU;
    std::vector< double > binWeight;
    for (const auto& [key, bin] : histogram)
    {
      binWeight.push_back(std::sqrt(bin[0]));
      binU.push_back({bin[1]/bin[0], bin[2]/bin[0]});
    }
    // The residual check takes at most 1024 of the pixels with signal, evenly spread, and up to 256
    // readout times, 16 per segment: in a test application of 93 small slice jobs it otherwise took
    // longer than the transforms. Its estimate changed by 1 % with 1024 pixels of 9268.
    std::vector< std::size_t > checkPixels;
    double checkEnergy = 0;
    {
      const auto numSignalPixels = static_cast<std::size_t>(std::count_if(energy.begin(), energy.end(), [](double e) { return e > 0; }));
      const std::size_t step = std::max<std::size_t>(1, numSignalPixels/1024);
      std::size_t count = 0;
      for (std::size_t r=0; r<numPixels; ++r)
        if (energy[r] > 0 && count++ % step == 0)
        {
          checkPixels.push_back(r);
          checkEnergy += energy[r];
        }
    }
    const auto residualSamples = [&distinct](unsigned int numSegments) { return std::min<std::size_t>({256, 16*std::size_t(numSegments), distinct.size()}); };

    // memory of the interpolators, the matrices of the fit (basis, SVD, pseudo-inverse), the check's
    // phase factors, the folded spectra and the buffers of the transforms
    const auto memory = [&](double L)
    {
      return sizeof(std::complex<double>)*(L*(distinct.size() + 4.0*numBins + checkPixels.size()) + numGroups*numTicks + parities*(nx + ny)*nkx + 2.0*numPixels);
    };
    // time of the interpolators of one block of readout times, of the SVD, of the residual check and of
    // the interpolators and transforms after it
    const auto blockTime = [&](double L) { return SEGMENT_WEIGHT_BLOCK*numBins*(0.36*L + 10.4); };
    const auto blocks = [](double count) { return std::ceil(count/SEGMENT_WEIGHT_BLOCK); };
    const auto fitTime = [&](double L) { return 7.0*numBins*L*std::min<double>(numBins, L); };
    const auto checkTime = [&](unsigned int numSegments)
    {
      const double L = numSegments;
      const double samples = residualSamples(numSegments);
      return checkPixels.size()*15.0*L + blocks(samples)*blockTime(L) + samples*checkPixels.size()*(1.36*L + 25);
    };
    const auto finalTime = [&](double L) { return blocks(distinct.size())*blockTime(L) + L*segmentTime; };
    // estimated time of the fits and checks done so far
    double spentTime = 0;

    std::vector< std::array<double, 2> > segmentV;
    // the pseudo-inverse of the fit (L x bins), real and imaginary parts apart for SolveSegmentWeights
    std::vector< double > pinvRe, pinvIm;
    // interpolators at the given v, a[l*count + i] for count v, in blocks in parallel
    const auto solveWeights = [&](const std::vector< std::array<double, 2> >& vs, std::vector< std::complex<double> >& a)
    {
      const std::size_t count = vs.size();
      a.assign(numSegments*count, std::complex<double>(0, 0));
      if (numBins == 0)
        return;
      parallelFor((count + SEGMENT_WEIGHT_BLOCK - 1)/SEGMENT_WEIGHT_BLOCK, [&](std::size_t block)
      {
        const std::size_t first = block*SEGMENT_WEIGHT_BLOCK;
        SolveSegmentWeights(pinvRe, pinvIm, binU, binWeight, &vs[first], std::min(SEGMENT_WEIGHT_BLOCK, count - first), &a[first], count);
      });
    };
    for (;;)
    {
      // The exact sum (in DynamicThreadedGenerateData) where the time segments would need more memory
      // than allowed (any count), or where they would take longer, with the fits and checks done so far
      // (automatic count; 0.8 leaves the close cases to the exact sum, which has no approximation error).
      // Counting the time spent bounds the time a slice can lose to a refinement that ends in the exact sum.
      const double check = requested > 0 ? 0.0 : checkTime(numSegments);
      if (memory(numSegments) > MAX_TIME_SEGMENT_BYTES ||
          (requested <= 0 && spentTime + fitTime(numSegments) + check + finalTime(numSegments) > 0.8*exactTime))
      {
        m_TimeSegmentsCapped = true;
        return;
      }

      // segment times: the ticks whose accumulated phase is closest to L values evenly spaced over the readout
      segmentV.clear();
      for (unsigned int l=0; l<numSegments; ++l)
      {
        const double target = numSegments==1 ? 0.5*phaseSpread : phaseSpread*l/(numSegments-1.0);
        auto next = std::lower_bound(accumulated.begin(), accumulated.end(), target);
        if (next == accumulated.end() || (next != accumulated.begin() && target - *(next-1) <= *next - target))
          --next;
        segmentV.push_back(v[order[static_cast<std::size_t>(next - accumulated.begin())]]);
      }

      if (numBins == 0)
        break;

      vnl_matrix< std::complex<double> > basis(static_cast<unsigned int>(numBins), numSegments);
      for (std::size_t b=0; b<numBins; ++b)
        for (unsigned int l=0; l<numSegments; ++l)
          basis(static_cast<unsigned int>(b), l) = binWeight[b]*std::polar(1.0, itk::Math::twopi*(binU[b][0]*segmentV[l][0] + binU[b][1]*segmentV[l][1]));
      // nearly equal segment phases (a weak field) make the basis ill-conditioned; drop those directions
      vnl_svd< std::complex<double> > svd(basis, -1e-8);
      const vnl_matrix< std::complex<double> > pseudoInverse = svd.pinverse();
      pinvRe.resize(numSegments*numBins);
      pinvIm.resize(numSegments*numBins);
      for (unsigned int l=0; l<numSegments; ++l)
        for (std::size_t b=0; b<numBins; ++b)
        {
          pinvRe[l*numBins + b] = pseudoInverse(l, static_cast<unsigned int>(b)).real();
          pinvIm[l*numBins + b] = pseudoInverse(l, static_cast<unsigned int>(b)).imag();
        }

      if (requested > 0)
        break;
      // energy-weighted residual of the pixels' phase factors at readout times evenly spread over the
      // distinct ones, both ends included; its RMS estimates the relative RMS error of the k-space
      std::vector< std::complex<double> > segmentPhase(checkPixels.size()*numSegments);
      parallelFor(checkPixels.size(), [&](std::size_t i)
      {
        const auto& ur = u[checkPixels[i]];
        for (unsigned int l=0; l<numSegments; ++l)
          segmentPhase[i*numSegments + l] = std::polar(1.0, itk::Math::twopi*(ur[0]*segmentV[l][0] + ur[1]*segmentV[l][1]));
      });
      const std::size_t numResidualSamples = residualSamples(numSegments);
      std::vector< std::array<double, 2> > sampleV(numResidualSamples);
      for (std::size_t s=0; s<numResidualSamples; ++s)
        sampleV[s] = v[distinct[numResidualSamples==1 ? 0 : s*(distinct.size()-1)/(numResidualSamples-1)]];
      std::vector< std::complex<double> > sampleWeights;
      solveWeights(sampleV, sampleWeights);
      std::vector< double > residual(numResidualSamples, 0.0);
      parallelFor(numResidualSamples, [&](std::size_t s)
      {
        const auto& vt = sampleV[s];
        double sum = 0;
        for (std::size_t i=0; i<checkPixels.size(); ++i)
        {
          const auto& ur = u[checkPixels[i]];
          const std::complex<double>* phase = &segmentPhase[i*numSegments];
          std::complex<double> fitted(0, 0);
          for (unsigned int l=0; l<numSegments; ++l)
            fitted += sampleWeights[l*numResidualSamples + s]*phase[l];
          sum += energy[checkPixels[i]]*std::norm(std::polar(1.0, itk::Math::twopi*(ur[0]*vt[0] + ur[1]*vt[1])) - fitted);
        }
        residual[s] = sum;
      });
      if (std::accumulate(residual.begin(), residual.end(), 0.0) <= 1e-8*checkEnergy*numResidualSamples)
        break;
      spentTime += fitTime(numSegments) + check;
      numSegments = static_cast<unsigned int>(std::ceil(1.25*numSegments));
    }

    // least-squares interpolators a(v) of every distinct v, a[l*numDistinct + i]
    const std::size_t numDistinct = distinct.size();
    std::vector< std::array<double, 2> > distinctV(numDistinct);
    for (std::size_t i=0; i<numDistinct; ++i)
      distinctV[i] = v[distinct[i]];
    std::vector< std::complex<double> > weights;
    solveWeights(distinctV, weights);

    // the interpolators of every k-space sample: the tick as in DynamicThreadedGenerateData
    std::vector< std::size_t > weightOfSample(numTicks, 0);
    for (std::size_t tick=0; tick<numTicks; ++tick)
    {
      const auto kIdx = m_ReadoutScheme->GetActualKspaceIndex(static_cast<int>(tick));
      std::size_t timeTick = tick;
      if (nky%2 == 0 && !m_Parameters->m_SignalGen.m_ReversePhase)
        timeTick = (tick + nkx)%numTicks;
      weightOfSample[static_cast<std::size_t>(kIdx[1])*nkx + static_cast<std::size_t>(kIdx[0])] = weightIndex[timeTick];
    }

    // The spectra of every group times the phase factor at every segment time, folded with the
    // interpolators as they are computed: m_CompartmentSpectra[g][k] = sum_l a_l(k) S_gl(k), summed in
    // segment order, so a sample needs no sum over segments and only one spectrum per group is kept
    // (with L spectra per group, slices of 256x256 with relaxation took the exact sum above about 60
    // segments for memory). The x factors are transposed (kx innermost) and real and imaginary parts
    // kept apart, so that the transforms run over kx with independent sums, which the compiler
    // vectorizes; the sums over x, y and l keep their order, so the result is that of the complex
    // arithmetic of the separable transform.
    std::vector< double > factorsXRe(parities*nx*nkx), factorsXIm(parities*nx*nkx);
    for (std::size_t p=0; p<parities; ++p)
      for (std::size_t kx=0; kx<nkx; ++kx)
        for (std::size_t x=0; x<nx; ++x)
        {
          factorsXRe[(p*nx + x)*nkx + kx] = factorsX[(p*nkx + kx)*nx + x].real();
          factorsXIm[(p*nx + x)*nkx + kx] = factorsX[(p*nkx + kx)*nx + x].imag();
        }
    m_CompartmentSpectra.assign(numGroups, std::vector< std::complex<double> >(numTicks, std::complex<double>(0, 0)));
    std::vector< double > phaseRe(numPixels), phaseIm(numPixels), rowsRe(parities*ny*nkx), rowsIm(parities*ny*nkx);
    for (unsigned int l=0; l<numSegments; ++l)
    {
      parallelFor(ny, [&](std::size_t y)
      {
        for (std::size_t r=y*nx; r<(y+1)*nx; ++r)
        {
          const std::complex<double> phase = std::polar(1.0, itk::Math::twopi*(u[r][0]*segmentV[l][0] + u[r][1]*segmentV[l][1]));
          phaseRe[r] = phase.real();
          phaseIm[r] = phase.imag();
        }
      });
      for (std::size_t g=0; g<numGroups; ++g)
      {
        // along x for every row and parity
        parallelFor(parities*ny, [&](std::size_t row)
        {
          const std::size_t p = row/ny;
          const std::size_t y = row%ny;
          double* sumRe = &rowsRe[row*nkx];
          double* sumIm = &rowsIm[row*nkx];
          std::fill(sumRe, sumRe + nkx, 0.0);
          std::fill(sumIm, sumIm + nkx, 0.0);
          for (std::size_t x=0; x<nx; ++x)
          {
            const std::size_t r = y*nx + x;
            const double scale = groups[g][r]*weight[r];
            const double inRe = scale*phaseRe[r];
            const double inIm = scale*phaseIm[r];
            const double* fRe = &factorsXRe[(p*nx + x)*nkx];
            const double* fIm = &factorsXIm[(p*nx + x)*nkx];
            for (std::size_t kx=0; kx<nkx; ++kx)
            {
              sumRe[kx] += inRe*fRe[kx] - inIm*fIm[kx];
              sumIm[kx] += inRe*fIm[kx] + inIm*fRe[kx];
            }
          }
        });
        // along y, then folded into the group's spectrum
        parallelFor(nky, [&](std::size_t ky)
        {
          const std::size_t p = ghosts ? ky%2 : 0;
          std::vector< double > sumRe(nkx, 0.0), sumIm(nkx, 0.0);
          for (std::size_t y=0; y<ny; ++y)
          {
            const double fRe = factorsY[ky*ny + y].real();
            const double fIm = factorsY[ky*ny + y].imag();
            const double* rRe = &rowsRe[(p*ny + y)*nkx];
            const double* rIm = &rowsIm[(p*ny + y)*nkx];
            for (std::size_t kx=0; kx<nkx; ++kx)
            {
              sumRe[kx] += fRe*rRe[kx] - fIm*rIm[kx];
              sumIm[kx] += fRe*rIm[kx] + fIm*rRe[kx];
            }
          }
          std::complex<double>* folded = &m_CompartmentSpectra[g][ky*nkx];
          for (std::size_t kx=0; kx<nkx; ++kx)
          {
            const std::complex<double> a = weights[l*numDistinct + weightOfSample[ky*nkx + kx]];
            folded[kx] += std::complex<double>(a.real()*sumRe[kx] - a.imag()*sumIm[kx], a.real()*sumIm[kx] + a.imag()*sumRe[kx]);
          }
        });
      }
    }
    m_NumSegments = numSegments;
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >::SolveSegmentWeights(const std::vector< double >& pinvRe, const std::vector< double >& pinvIm,
                                                           const std::vector< std::array<double, 2> >& binU, const std::vector< double >& binWeight,
                                                           const std::array<double, 2>* v, std::size_t count, std::complex<double>* a, std::size_t stride)
  {
    // The sums of a block of v are independent, so the innermost loop runs over the block and the
    // compiler vectorizes it (more than twice as fast as one v at a time). A whole block is computed
    // also for fewer v, so every v gets the same arithmetic, whatever the blocks.
    constexpr std::size_t block = SEGMENT_WEIGHT_BLOCK;
    const std::size_t numBins = binU.size();
    const std::size_t numSegments = pinvRe.size()/numBins;
    std::vector< double > targetRe(numBins*block), targetIm(numBins*block);
    for (std::size_t b=0; b<numBins; ++b)
      for (std::size_t j=0; j<block; ++j)
      {
        const auto& vj = v[std::min(j, count-1)];
        const std::complex<double> phase = std::polar(1.0, itk::Math::twopi*(binU[b][0]*vj[0] + binU[b][1]*vj[1]));
        targetRe[b*block + j] = binWeight[b]*phase.real();
        targetIm[b*block + j] = binWeight[b]*phase.imag();
      }
    std::array< double, block > sumRe, sumIm;
    for (std::size_t l=0; l<numSegments; ++l)
    {
      sumRe.fill(0.0);
      sumIm.fill(0.0);
      const double* rowRe = &pinvRe[l*numBins];
      const double* rowIm = &pinvIm[l*numBins];
      for (std::size_t b=0; b<numBins; ++b)
      {
        const double pRe = rowRe[b];
        const double pIm = rowIm[b];
        const double* tRe = &targetRe[b*block];
        const double* tIm = &targetIm[b*block];
        for (std::size_t j=0; j<block; ++j)
        {
          sumRe[j] += pRe*tRe[j] - pIm*tIm[j];
          sumIm[j] += pRe*tIm[j] + pIm*tRe[j];
        }
      }
      for (std::size_t j=0; j<count; ++j)
        a[l*stride + j] = std::complex<double>(sumRe[j], sumIm[j]);
    }
  }

  template< class ScalarType >
  std::unique_ptr<mitk::AcquisitionType> KspaceImageFilter< ScalarType >::CreateReadoutScheme(const mitk::FiberfoxParameters* parameters)
  {
    switch (parameters->m_SignalGen.m_AcquisitionType)
    {
      case mitk::SignalGenerationParameters::ConventionalSpinEcho:
        return std::make_unique<mitk::ConventionalSpinEcho>(parameters);
      case mitk::SignalGenerationParameters::FastSpinEcho:
        return std::make_unique<mitk::FastSpinEcho>(parameters);
      case mitk::SignalGenerationParameters::SingleShotEpi:
      default:
        return std::make_unique<mitk::SingleShotEpi>(parameters);
    }
  }

  template< class ScalarType >
  float KspaceImageFilter< ScalarType >::CoilSensitivity(VectorType& pos) const
  {
    // The coil ring moves with the excited slice. A local copy: all work units call this at once,
    // and for oblique images pos[2] differs between pixels of one slice.
    VectorType coilPosition = m_CoilPosition;
    coilPosition[2] = pos[2];

    switch (m_Parameters->m_SignalGen.m_CoilSensitivityProfile)
    {
      case mitk::SignalGenerationParameters::COIL_CONSTANT:
      return 1;
      case mitk::SignalGenerationParameters::COIL_LINEAR:
      {
        VectorType diff = pos-coilPosition;
        float sens = diff.GetNorm()*m_CoilSensitivityFactor + 1;
        if (sens<0)
          sens = 0;
        return sens;
      }
      case mitk::SignalGenerationParameters::COIL_EXPONENTIAL:
      {
        VectorType diff = pos-coilPosition;
        float dist = static_cast<float>(diff.GetNorm());
        return std::exp(-dist*m_CoilSensitivityFactor);
      }
      default:
      return 1;
    }
  }

  template< class ScalarType >
  bool KspaceImageFilter< ScalarType >::IsSkippedLine(int ky, int numLines, double partialFourier, bool reversePhase)
  {
    const int half = numLines/2;  // index of the k-space centre
    const int acquired = std::clamp(static_cast<int>(std::lround(partialFourier*numLines)), half+1, numLines);
    if (reversePhase)
      return ky >= acquired;
    const bool even = numLines%2 == 0;
    if (even && ky == 0)
      return false;
    return ky < numLines - acquired + (even ? 1 : 0);
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >
  ::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread )
  {
    typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);

    // precalculate shifts for DFT
    float x_shift = 0;
    float y_shift = 0;
    if (static_cast<int>(xMax)%2==1)
        x_shift = (xMax-1)/2;
    else
        x_shift = xMax/2;
    if (static_cast<int>(yMax)%2==1)
        y_shift = (yMax-1)/2;
    else
        y_shift = yMax/2;

    float kx_shift = 0;
    float ky_shift = 0;
    if (static_cast<int>(kxMax)%2==1)
        kx_shift = (kxMax-1)/2;
    else
        kx_shift = kxMax/2;
    if (static_cast<int>(kyMax)%2==1)
        ky_shift = (kyMax-1)/2;
    else
        ky_shift = kyMax/2;

    // Largest noise-free sample of this work unit (the spike amplitude). With one work unit the first
    // of equal maxima in readout order wins; the merge below keeps that rule for any number of units.
    std::complex<ScalarType> spike(0,0);
    long spikeTick = -1;

    std::complex<ScalarType> zero = std::complex<ScalarType>(0, 0);
    // the spectra of the compartments, or with time segments of the groups, folded (ComputeTimeSegmentation)
    const bool precomputed = !m_CompartmentSpectra.empty();
    // per work unit, not per sample
    std::vector< float > relaxFactor(m_CompartmentImages.size(), 1.0f);
    std::vector< const ScalarType* > compartmentSignals;
    for (const auto& compartment : m_CompartmentImages)
      compartmentSignals.push_back(compartment->GetBufferPointer());
    const long inputWidth = static_cast<long>(xMax);
    const long numInputPixels = static_cast<long>(m_CompartmentImages[0]->GetLargestPossibleRegion().GetNumberOfPixels());

    while( !oit.IsAtEnd() )
    {
      int tick = oit.GetIndex()[1] * kxMax + oit.GetIndex()[0];
      const long readoutTick = tick;

      // get current k-space index (depends on the chosen k-space readout scheme)
      itk::Index< 2 > kIdx = m_ReadoutScheme->GetActualKspaceIndex(tick);

      // we have to adjust the ticks to obtain correct times since the DFT is not completely symmetric in the even  number of lines case
      if (static_cast<int>(kyMax)%2 == 0 && !m_Parameters->m_SignalGen.m_ReversePhase)
      {
        tick += kxMax;
        tick %= static_cast<int>(numPix);
      }

      // partial fourier: the lines read last are not acquired
      if (IsSkippedLine(kIdx[1], static_cast<int>(kyMax), m_Parameters->m_SignalGen.m_PartialFourier, m_Parameters->m_SignalGen.m_ReversePhase))
      {
        outputImage->SetPixel(kIdx, zero);
        ++oit;
        continue;
      }
      if (m_StoreTimings)
        m_TickImage->SetPixel(kIdx, tick);

      // gibbs ringing by setting high frequencies to zero (alternative to using smaller k-space than input image space)
      if (m_Parameters->m_SignalGen.m_DoAddGibbsRinging && m_Parameters->m_SignalGen.m_ZeroRinging>0)
      {
        if (kIdx[0] < ringing_lines_x || kIdx[1] < ringing_lines_y ||
            kIdx[0] >= kxMax - ringing_lines_x || kIdx[1] >= kyMax - ringing_lines_y)
        {
          outputImage->SetPixel(kIdx, zero);
          ++oit;
          continue;
        }
      }

      // time from maximum echo
      float t = m_ReadoutScheme->GetTimeFromMaxEcho(tick);

      // calculate eddy current decay factor
      float eddyDecay = 0;
      if ( m_Parameters->m_Misc.m_DoAddEddyCurrents && m_Parameters->m_SignalGen.m_EddyStrength>0 && !m_IsBaseline)
      {
        // time passed since k-space readout started
        float tRead = m_ReadoutScheme->GetTimeFromLastDiffusionGradient(tick);
        eddyDecay = std::exp(-tRead/m_Parameters->m_SignalGen.m_Tau ) * t/1000; // time in seconds here
      }

      // calcualte signal relaxation factors
      if ( m_Parameters->m_SignalGen.m_DoSimulateRelaxation)
      {
        // time passes since application of the RF pulse
        float tRf = m_ReadoutScheme->GetTimeFromRf(tick);
        if (m_StoreTimings)
            m_RfImage->SetPixel(kIdx, tRf);

        for (unsigned int i=0; i<m_CompartmentImages.size(); i++)
        {
          // account for T2 relaxation (how much transverse magnetization is left since applicatiohn of RF pulse?)
          relaxFactor[i] = m_T1Relax[i] * std::exp(-tRf/m_T2[i] -fabs(t)/ m_Parameters->m_SignalGen.m_tInhom);
        }
      }

      std::complex<ScalarType> s(0,0);
      if (precomputed)
      {
        const std::size_t k = static_cast<std::size_t>(kIdx[1])*static_cast<std::size_t>(kxMax) + static_cast<std::size_t>(kIdx[0]);
        std::complex<double> sum(0, 0);
        for (std::size_t i=0; i<m_CompartmentSpectra.size(); ++i)
          sum += static_cast<double>(relaxFactor[i]) * m_CompartmentSpectra[i][k];
        s = std::complex<ScalarType>(sum / static_cast<double>(numPix));
      }
      else
      {
        // shift k for DFT: (0 -- N) --> (-N/2 -- N/2)
        float kx = kIdx[0] - kx_shift;
        float ky = kIdx[1] - ky_shift;

        // add ghosting by adding gradient delay induced offset
        if (m_Parameters->m_Misc.m_DoAddGhosts)
        {
          if (kIdx[1]%2 == 1)
            kx -= m_Parameters->m_SignalGen.m_KspaceLineOffset;
          else
            kx += m_Parameters->m_SignalGen.m_KspaceLineOffset;
        }

        // pull stuff out of inner loop
        t /= 1000; // time in seconds
        kx /= xMax;
        ky /= yMaxFov;

        // calculate signal s at k-space position (kx, ky)
        for (long pixel=0; pixel<numInputPixels; ++pixel)
        {
          typename InputImageType::IndexType input_idx;
          input_idx[0] = pixel % inputWidth;
          input_idx[1] = pixel / inputWidth;

          // shift x,y for DFT: (0 -- N) --> (-N/2 -- N/2)
          float x = input_idx[0] - x_shift;
          float y = input_idx[1] - y_shift;

          // sum compartment signals and simulate relaxation
          ScalarType f_real = 0;
          for (unsigned int i=0; i<m_CompartmentImages.size(); i++)
            if ( m_Parameters->m_SignalGen.m_DoSimulateRelaxation)
              f_real += compartmentSignals[i][pixel] * relaxFactor[i];
            else
              f_real += compartmentSignals[i][pixel];

          // vector from image center to current position (in meter)
          // only necessary for eddy currents and non-constant coil sensitivity
          VectorType pos;
          if ((m_Parameters->m_Misc.m_DoAddEddyCurrents && m_Parameters->m_SignalGen.m_EddyStrength>0 && !m_IsBaseline) ||
              m_Parameters->m_SignalGen.m_CoilSensitivityProfile!=mitk::SignalGenerationParameters::COIL_CONSTANT)
          {
            pos[0] = x; pos[1] = y; pos[2] = m_Z;
            pos = m_Transform*pos;
          }

          if (m_Parameters->m_SignalGen.m_CoilSensitivityProfile!=mitk::SignalGenerationParameters::COIL_CONSTANT)
            f_real *= CoilSensitivity(pos);

          // simulate eddy currents and other distortions
          float phi = 0;   // phase shift
          if (  m_Parameters->m_Misc.m_DoAddEddyCurrents && m_Parameters->m_SignalGen.m_EddyStrength>0 && !m_IsBaseline)
          {
            // duration (tRead) already included in "eddyDecay"
            phi += (m_EddyFieldGradient[0]*pos[0]+m_EddyFieldGradient[1]*pos[1]+m_EddyFieldGradient[2]*pos[2]) * eddyDecay;
          }

          // simulate distortions
          if (m_Parameters->m_Misc.m_DoAddDistortions)
          {
            if (m_MovedFmap.IsNotNull())    // if we have headmotion, use moved map
              phi += m_MovedFmap->GetPixel(input_idx) * t;
            else if (m_Parameters->m_SignalGen.m_FrequencyMap.IsNotNull())
            {
              itk::Image<float, 3>::IndexType index; index[0] = input_idx[0]; index[1] = input_idx[1]; index[2] = m_Zidx;
              phi += m_Parameters->m_SignalGen.m_FrequencyMap->GetPixel(index) * t;
            }
          }

          // if signal comes from outside FOV, mirror it back (wrap-around artifact - aliasing
          if (m_Parameters->m_Misc.m_DoAddAliasing)
          {
            if (y<-yMaxFov_half)
              y += yMaxFov;
            else if (y>yMaxFov_half)
              y -= yMaxFov;
          }

          // actual DFT term
          std::complex<ScalarType> f(f_real * m_Parameters->m_SignalGen.m_SignalScale, 0);
          s += f * std::exp( std::complex<ScalarType>(0, itk::Math::twopi * (kx*x + ky*y + phi )) );
        }
        s /= numPix;
      }

      if (m_SpikesPerSlice>0 && sqrt(s.imag()*s.imag()+s.real()*s.real()) > sqrt(spike.imag()*spike.imag()+spike.real()*spike.real()) )
      {
        spike = s;
        spikeTick = readoutTick;
      }

      s += outputImage->GetPixel(kIdx); // add precalculated noise
      outputImage->SetPixel(kIdx, s);
      m_KSpaceImage->SetPixel(kIdx, sqrt(s.imag()*s.imag()+s.real()*s.real()) );

      ++oit;
    }

    if (spikeTick>=0)
    {
      std::lock_guard<std::mutex> lock(m_SpikeMutex);
      const auto magnitude = sqrt(spike.imag()*spike.imag()+spike.real()*spike.real());
      const auto maxMagnitude = sqrt(m_Spike.imag()*m_Spike.imag()+m_Spike.real()*m_Spike.real());
      if (magnitude > maxMagnitude || (magnitude == maxMagnitude && spikeTick < m_SpikeTick))
      {
        m_Spike = spike;
        m_SpikeTick = spikeTick;
      }
    }
  }

  template< class ScalarType >
  void KspaceImageFilter< ScalarType >
  ::AfterThreadedGenerateData()
  {
    typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    int kxMax = outputImage->GetLargestPossibleRegion().GetSize(0);  // k-space size in x-direction
    int kyMax = outputImage->GetLargestPossibleRegion().GetSize(1);  // k-space size in y-direction

    ImageRegionIterator< OutputImageType > oit(outputImage, outputImage->GetLargestPossibleRegion());
    while( !oit.IsAtEnd() ) // use hermitian k-space symmetry to fill empty k-space parts resulting from partial fourier acquisition
    {
      int tick = oit.GetIndex()[1] * kxMax + oit.GetIndex()[0];
      auto kIdx = m_ReadoutScheme->GetActualKspaceIndex(tick);

      if (IsSkippedLine(kIdx[1], kyMax, m_Parameters->m_SignalGen.m_PartialFourier, m_Parameters->m_SignalGen.m_ReversePhase))
      {
        // calculate symmetric index
        auto sym = m_ReadoutScheme->GetSymmetricIndex(kIdx);

        // use complex conjugate of symmetric index value at current index
        std::complex<ScalarType> s = outputImage->GetPixel(sym);
        s = std::complex<ScalarType>(s.real(), -s.imag());
        outputImage->SetPixel(kIdx, s);

        m_KSpaceImage->SetPixel(kIdx, sqrt(s.imag()*s.imag()+s.real()*s.real()) );
      }
      ++oit;
    }

    m_Spike *=  m_Parameters->m_SignalGen.m_SpikeAmplitude;
    itk::Index< 2 > spikeIdx;
    for (unsigned int i=0; i<m_SpikesPerSlice; i++)
    {
      spikeIdx[0] = m_RandGen->GetIntegerVariate()%kxMax;
      spikeIdx[1] = m_RandGen->GetIntegerVariate()%kyMax;
      outputImage->SetPixel(spikeIdx, m_Spike);
      m_SpikeLog += "[" + boost::lexical_cast<std::string>(spikeIdx[0]) + "," + boost::lexical_cast<std::string>(spikeIdx[1]) + "," + boost::lexical_cast<std::string>(m_Zidx) + "] Magnitude: " + boost::lexical_cast<std::string>(m_Spike.real()) + "+" + boost::lexical_cast<std::string>(m_Spike.imag()) + "i\n";
    }
    m_ReadoutScheme.reset();
  }
}
#endif
