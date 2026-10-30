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
    m_SegmentWeights.clear();
    m_NumSegments = 0;
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
   * segment times t_l (ticks evenly spaced in t),
   *   exp(i 2 pi phi(r, j)) ~ sum_l a_l(j) exp(i 2 pi phi(r, t_l)),
   * so s(k_j) ~ sum_l a_l(j) S_l(k_j), where S_l is the spectrum of the signal times the phase factor
   * at t_l: L separable DFTs per compartment instead of a sum over all pixels per sample. a(j) is the
   * least-squares fit over the pixels' values of u, binned to a histogram weighted by signal energy:
   * minimise sum_b w_b |exp(i 2 pi u_b . v(j)) - sum_l a_l exp(i 2 pi u_b . v_l)|^2. The error falls
   * steeply once L exceeds the phase spread range(u) . range(v) in cycles (the time-bandwidth product
   * of the fit), hence the automatic count below.
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

    // Phase spread: the largest phase difference (cycles) between two pixels with signal over the
    // readout. The automatic count is the spread plus 8; the k-space test and the measurements of
    // issue #47 found a relative RMS error below 0.1 % of the image with it.
    std::array<double, 2> lowU = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 2> highU = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (std::size_t r=0; r<numPixels; ++r)
      if (energy[r] > 0)
        for (unsigned int d=0; d<2; ++d)
        {
          lowU[d] = std::min(lowU[d], u[r][d]);
          highU[d] = std::max(highU[d], u[r][d]);
        }
    double phaseSpread = 0;
    for (unsigned int d=0; d<2; ++d)
      if (highU[d] >= lowU[d])
        phaseSpread += (highU[d]-lowU[d])*span[d];
    const unsigned int numSegments = m_Parameters->m_SignalGen.m_KspaceTimeSegments > 0 ?
          static_cast<unsigned int>(m_Parameters->m_SignalGen.m_KspaceTimeSegments) : static_cast<unsigned int>(std::ceil(phaseSpread)) + 8;

    // segment times: ticks whose readout time is closest to L times evenly spaced over the readout
    const auto [minTime, maxTime] = std::minmax_element(echoTime.begin(), echoTime.end());
    std::vector< std::array<double, 2> > segmentV;
    for (unsigned int l=0; l<numSegments; ++l)
    {
      const double target = numSegments==1 ? 0.5*(*minTime + *maxTime) : *minTime + (*maxTime - *minTime)*l/(numSegments-1.0);
      std::size_t best = 0;
      for (std::size_t tick=1; tick<numTicks; ++tick)
        if (std::fabs(echoTime[tick] - target) < std::fabs(echoTime[best] - target))
          best = tick;
      segmentV.push_back(v[best]);
    }

    // Histogram of u weighted by signal energy. The interpolators fitted at the bins hold for the
    // pixels between them if the bins sample the phase finely: at least 8 bins per cycle of phase
    // spread and 16 over the range of u along each axis (coarser if that gives more than maxBins bins).
    const std::size_t maxBins = 512;
    std::map< std::array<long long, 2>, std::array<double, 3> > histogram;
    for (double coarsening = 1; ; coarsening *= 2)
    {
      std::array<double, 2> binWidth = {0.0, 0.0};
      bool single = true;
      for (unsigned int d=0; d<2; ++d)
        if (highU[d] > lowU[d] && span[d] > 0)
        {
          const double numBins = std::ceil(std::max(8*(highU[d]-lowU[d])*span[d], 16.0)/coarsening);
          binWidth[d] = (highU[d]-lowU[d])/numBins;
          single = single && numBins <= 1;
        }
      histogram.clear();
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
      if (histogram.size() <= maxBins || single)
        break;
    }

    // least-squares interpolators a(tick), the same for ticks with the same v (spin echo lines)
    m_SegmentWeights.assign(numTicks*numSegments, std::complex<double>(0, 0));
    if (!histogram.empty())
    {
      const std::size_t numBins = histogram.size();
      std::vector< std::array<double, 2> > binU;
      vnl_vector< double > binWeight(static_cast<unsigned int>(numBins));
      for (const auto& [key, bin] : histogram)
      {
        binWeight[static_cast<unsigned int>(binU.size())] = std::sqrt(bin[0]);
        binU.push_back({bin[1]/bin[0], bin[2]/bin[0]});
      }
      vnl_matrix< std::complex<double> > basis(static_cast<unsigned int>(numBins), numSegments);
      for (std::size_t b=0; b<numBins; ++b)
        for (unsigned int l=0; l<numSegments; ++l)
          basis(static_cast<unsigned int>(b), l) = binWeight[static_cast<unsigned int>(b)]*std::polar(1.0, itk::Math::twopi*(binU[b][0]*segmentV[l][0] + binU[b][1]*segmentV[l][1]));
      // nearly equal segment phases (a weak field) make the basis ill-conditioned; drop those directions
      vnl_svd< std::complex<double> > svd(basis, -1e-8);
      const vnl_matrix< std::complex<double> > pseudoInverse = svd.pinverse();

      std::map< std::array<double, 2>, std::size_t > solved;
      vnl_vector< std::complex<double> > target(static_cast<unsigned int>(numBins));
      for (std::size_t tick=0; tick<numTicks; ++tick)
      {
        const auto found = solved.find(v[tick]);
        if (found != solved.end())
        {
          std::copy_n(&m_SegmentWeights[found->second*numSegments], numSegments, &m_SegmentWeights[tick*numSegments]);
          continue;
        }
        for (std::size_t b=0; b<numBins; ++b)
          target[static_cast<unsigned int>(b)] = binWeight[static_cast<unsigned int>(b)]*std::polar(1.0, itk::Math::twopi*(binU[b][0]*v[tick][0] + binU[b][1]*v[tick][1]));
        const vnl_vector< std::complex<double> > a = pseudoInverse*target;
        std::copy_n(a.data_block(), numSegments, &m_SegmentWeights[tick*numSegments]);
        solved.emplace(v[tick], tick);
      }
    }

    // spectra of every group times the phase factor at every segment time
    std::vector< std::complex<double> > input(numPixels);
    std::vector< std::complex<double> > rows(parities*ny*nkx);
    for (const auto& group : groups)
      for (unsigned int l=0; l<numSegments; ++l)
      {
        for (std::size_t r=0; r<numPixels; ++r)
          input[r] = (group[r]*weight[r]) * std::polar(1.0, itk::Math::twopi*(u[r][0]*segmentV[l][0] + u[r][1]*segmentV[l][1]));

        for (std::size_t p=0; p<parities; ++p)
          for (std::size_t y=0; y<ny; ++y)
            for (std::size_t kx=0; kx<nkx; ++kx)
            {
              std::complex<double> sum(0, 0);
              const std::complex<double>* factors = &factorsX[(p*nkx + kx)*nx];
              const std::complex<double>* line = &input[y*nx];
              for (std::size_t x=0; x<nx; ++x)
                sum += line[x] * factors[x];
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
    m_NumSegments = numSegments;
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
    const bool precomputed = m_NumSegments==0 && !m_CompartmentSpectra.empty();
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
      else if (m_NumSegments>0)
      {
        // time-segmented approximation of the phase, see ComputeTimeSegmentation
        const std::size_t k = static_cast<std::size_t>(kIdx[1])*static_cast<std::size_t>(kxMax) + static_cast<std::size_t>(kIdx[0]);
        const std::complex<double>* a = &m_SegmentWeights[static_cast<std::size_t>(tick)*m_NumSegments];
        const std::size_t numGroups = m_CompartmentSpectra.size()/m_NumSegments;
        std::complex<double> sum(0, 0);
        for (std::size_t g=0; g<numGroups; ++g)
        {
          std::complex<double> group(0, 0);
          for (unsigned int l=0; l<m_NumSegments; ++l)
            group += a[l] * m_CompartmentSpectra[g*m_NumSegments + l][k];
          // without relaxation one group holds all compartments and relaxFactor[0] stays 1
          sum += static_cast<double>(relaxFactor[g]) * group;
        }
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
