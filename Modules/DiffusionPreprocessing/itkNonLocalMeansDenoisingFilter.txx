/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkNonLocalMeansDenoisingFilter_txx
#define itkNonLocalMeansDenoisingFilter_txx

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <mitkLog.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace itk {


template< class TPixelType >
NonLocalMeansDenoisingFilter< TPixelType >
::NonLocalMeansDenoisingFilter()
  : m_SearchRadius(4),
    m_ComparisonRadius(1),
    m_UseJointInformation(false),
    m_UseRicianAdaption(false),
    m_Variance(1),
    m_Mask(nullptr)
{
  this->SetNumberOfRequiredInputs( 1 );

  
}

template< class TPixelType >
void
NonLocalMeansDenoisingFilter< TPixelType >
::BeforeThreadedGenerateData()
{

  MITK_INFO << "SearchRadius: " << m_SearchRadius;
  MITK_INFO << "ComparisonRadius: " << m_ComparisonRadius;
  MITK_INFO << "Noisevariance: " << m_Variance;
  MITK_INFO << "Use Rician Adaption: " << std::boolalpha << m_UseRicianAdaption;
  MITK_INFO << "Use Joint Information: " << std::boolalpha << m_UseJointInformation;


  // The weights are exp(-distance / variance).
  if (!(m_Variance > 0))
    itkExceptionMacro("The noise variance must be positive.");

  typename InputImageType::Pointer inputImagePointer = static_cast< InputImageType * >( this->ProcessObject::GetInput(0) );
  if (m_Mask.IsNull())
  {
    // Without a mask every voxel is denoised. The full mask is kept for this run only, so that a
    // later run on another input does not inherit it.
    m_UsedMask = MaskImageType::New();
    m_UsedMask->SetRegions(inputImagePointer->GetLargestPossibleRegion());
    m_UsedMask->Allocate();
    m_UsedMask->FillBuffer(1);
  }
  else
  {
    // Mask and image are iterated voxel by voxel.
    if (m_Mask->GetLargestPossibleRegion() != inputImagePointer->GetLargestPossibleRegion())
      itkExceptionMacro("The mask does not have the size of the image.");
    m_UsedMask = m_Mask;

    // Calculation of the smallest masked region

    typename OutputImageType::Pointer outputImage =
            static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
    ImageRegionIterator< OutputImageType > oit(outputImage, inputImagePointer->GetLargestPossibleRegion());
    oit.GoToBegin();
    ImageRegionIterator< MaskImageType > mit(m_Mask, m_Mask->GetLargestPossibleRegion());
    mit.GoToBegin();
    typename MaskImageType::IndexType minIndex = m_Mask->GetLargestPossibleRegion().GetUpperIndex();
    typename MaskImageType::IndexType maxIndex = m_Mask->GetLargestPossibleRegion().GetIndex();
    bool foundMaskedVoxel = false;
    typename OutputImageType::PixelType outpix;
    outpix.SetSize(inputImagePointer->GetVectorLength());
    outpix.Fill(0);
    while (!mit.IsAtEnd())
    {

      if (mit.Get())
      {
        // calculation of the start & end index of the smallest masked region
        foundMaskedVoxel = true;
        minIndex[0] = minIndex[0] < mit.GetIndex()[0] ? minIndex[0] : mit.GetIndex()[0];
        minIndex[1] = minIndex[1] < mit.GetIndex()[1] ? minIndex[1] : mit.GetIndex()[1];
        minIndex[2] = minIndex[2] < mit.GetIndex()[2] ? minIndex[2] : mit.GetIndex()[2];

        maxIndex[0] = maxIndex[0] > mit.GetIndex()[0] ? maxIndex[0] : mit.GetIndex()[0];
        maxIndex[1] = maxIndex[1] > mit.GetIndex()[1] ? maxIndex[1] : mit.GetIndex()[1];
        maxIndex[2] = maxIndex[2] > mit.GetIndex()[2] ? maxIndex[2] : mit.GetIndex()[2];
      }
      else
      {
        oit.Set(outpix);
      }
      ++mit;
      ++oit;
    }

    // calculation of the masked region; an empty mask leaves the output all zero (set above), and
    // the threads find no masked voxel in the full region
    if (foundMaskedVoxel)
    {
      typename OutputImageType::SizeType size;
      size[0] = maxIndex[0] - minIndex[0] + 1;
      size[1] = maxIndex[1] - minIndex[1] + 1;
      size[2] = maxIndex[2] - minIndex[2] + 1;

      typename OutputImageType::RegionType region (minIndex, size);

      outputImage->SetRequestedRegion(region);
    }
  }

  m_CurrentVoxelCount = 0;
}

template< class TPixelType >
void
NonLocalMeansDenoisingFilter< TPixelType >
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread )
{
  OutputImageType* outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  const InputImageType* inputImage = static_cast< const InputImageType * >(this->ProcessObject::GetInput(0));

  ImageRegionIteratorWithIndex< OutputImageType > oit(outputImage, outputRegionForThread);
  ImageRegionConstIterator< MaskImageType > mit(m_UsedMask, outputRegionForThread);

  const unsigned int vectorLength = inputImage->GetVectorLength();
  const typename InputImageType::IndexType lower = inputImage->GetLargestPossibleRegion().GetIndex();
  const typename InputImageType::IndexType upper = inputImage->GetLargestPossibleRegion().GetUpperIndex();

  // Patches are read from the pixel buffer directly (offsets as in Image::ComputeOffset), and the
  // search window and patches are clipped to the image once instead of testing every element:
  // GetPixel and IsInside per patch element and channel took most of the run time. Every sum runs
  // over the search and patch positions in x-y-z loop order, the order the NLM references in
  // MITK-Data were computed in.
  const TPixelType* buffer = inputImage->GetBufferPointer();
  const typename InputImageType::IndexType bufferIndex = inputImage->GetBufferedRegion().GetIndex();
  const OffsetValueType* offsetTable = inputImage->GetOffsetTable();
  auto pixelAt = [&](IndexValueType x, IndexValueType y, IndexValueType z)
  {
    return buffer + ((x - bufferIndex[0]) + (y - bufferIndex[1]) * offsetTable[1] + (z - bufferIndex[2]) * offsetTable[2]) * vectorLength;
  };
  const OffsetValueType zStep = offsetTable[2] * vectorLength;

  // Work buffers, reused for every voxel. Without joint information, every channel has its own
  // weights, computed from that channel alone.
  const unsigned int weightsPerCandidate = m_UseJointInformation ? 1 : vectorLength;
  std::vector<const TPixelType*> candidates;
  std::vector<double> weights;
  std::vector<double> distances(weightsPerCandidate);
  std::vector<double> weightSums(weightsPerCandidate);
  std::vector<double> sums(vectorLength);
  typename OutputImageType::PixelType outpix(vectorLength);

  for (; !oit.IsAtEnd(); ++oit, ++mit)
  {
    if (mit.Get() != 0 && !this->GetAbortGenerateData())
    {
      const typename OutputImageType::IndexType index = oit.GetIndex();
      candidates.clear();
      weights.clear();
      std::fill(weightSums.begin(), weightSums.end(), 0.0);

      for (IndexValueType x = std::max(index[0] - m_SearchRadius, lower[0]); x <= std::min(index[0] + m_SearchRadius, upper[0]); ++x)
      {
        for (IndexValueType y = std::max(index[1] - m_SearchRadius, lower[1]); y <= std::min(index[1] + m_SearchRadius, upper[1]); ++y)
        {
          for (IndexValueType z = std::max(index[2] - m_SearchRadius, lower[2]); z <= std::min(index[2] + m_SearchRadius, upper[2]); ++z)
          {
            // Compare the patches around index and (x, y, z) at the offsets d for which both
            // index + d and (x, y, z) + d lie in the image.
            const IndexValueType candidate[3] = {x, y, z};
            IndexValueType dMin[3];
            IndexValueType dMax[3];
            for (unsigned int k = 0; k < 3; ++k)
            {
              dMin[k] = std::max({static_cast<IndexValueType>(-m_ComparisonRadius), lower[k] - index[k], lower[k] - candidate[k]});
              dMax[k] = std::min({static_cast<IndexValueType>(m_ComparisonRadius), upper[k] - index[k], upper[k] - candidate[k]});
            }

            std::fill(distances.begin(), distances.end(), 0.0);
            for (IndexValueType dx = dMin[0]; dx <= dMax[0]; ++dx)
            {
              for (IndexValueType dy = dMin[1]; dy <= dMax[1]; ++dy)
              {
                const TPixelType* valuesI = pixelAt(index[0] + dx, index[1] + dy, index[2] + dMin[2]);
                const TPixelType* valuesJ = pixelAt(x + dx, y + dy, z + dMin[2]);
                for (IndexValueType dz = dMin[2]; dz <= dMax[2]; ++dz, valuesI += zStep, valuesJ += zStep)
                {
                  // In double: the square of an int difference overflows for |diff| > 46340, and
                  // subtracting in TPixelType would wrap.
                  if (m_UseJointInformation)
                  {
                    double squaredNorm = 0;
                    for (unsigned int c = 0; c < vectorLength; ++c)
                    {
                      const double diff = static_cast<double>(valuesI[c]) - static_cast<double>(valuesJ[c]);
                      squaredNorm += diff*diff;
                    }
                    distances[0] += squaredNorm;
                  }
                  else
                  {
                    for (unsigned int c = 0; c < vectorLength; ++c)
                    {
                      const double diff = static_cast<double>(valuesI[c]) - static_cast<double>(valuesJ[c]);
                      distances[c] += diff*diff;
                    }
                  }
                }
              }
            }

            // weight all neighborhoods
            double size = static_cast<double>((dMax[0] - dMin[0] + 1) * (dMax[1] - dMin[1] + 1) * (dMax[2] - dMin[2] + 1));
            candidates.push_back(pixelAt(x, y, z));
            if (m_UseJointInformation)
            {
              size *= vectorLength + 1;
              const double w = std::exp( - (distances[0] / size) / m_Variance);
              weights.push_back(w);
              weightSums[0] += w;
            }
            else
            {
              for (unsigned int c = 0; c < vectorLength; ++c)
              {
                const double w = std::exp( - distances[c] / size / m_Variance);
                weights.push_back(w);
                weightSums[c] += w;
              }
            }
          }
        }
      }

      std::fill(sums.begin(), sums.end(), 0.0);
      for (std::size_t n = 0; n < candidates.size(); ++n)
      {
        const TPixelType* pixelJ = candidates[n];
        for (unsigned int c = 0; c < vectorLength; ++c)
        {
          const double value = m_UseRicianAdaption ? static_cast<double>(pixelJ[c] * pixelJ[c]) : static_cast<double>(pixelJ[c]);
          if (m_UseJointInformation)
            sums[c] += (weights[n] / weightSums[0]) * value;
          else
            sums[c] += (weights[n * vectorLength + c] / weightSums[c]) * value;
        }
      }

      for (unsigned int c = 0; c < vectorLength; ++c)
      {
        double a = sums[c];
        if (m_UseRicianAdaption)
        {
          a -= 2 * m_Variance;
        }
        if (a < 0)
        {
          a = 0;
        }
        if (m_UseRicianAdaption)
        {
          outpix[c] = static_cast<TPixelType>(std::floor(std::sqrt(a) + 0.5));
        }
        else
        {
          outpix[c] = static_cast<TPixelType>(std::floor(a + 0.5));
        }
      }
    }
    else
    {
      outpix.Fill(0);
    }

    oit.Set(outpix);
    ++m_CurrentVoxelCount;
  }
}

template< class TPixelType >
void
NonLocalMeansDenoisingFilter< TPixelType >
::AfterThreadedGenerateData()
{
  m_UsedMask = nullptr;
}

template< class TPixelType >
void NonLocalMeansDenoisingFilter< TPixelType >::SetInputImage(const InputImageType* image)
{
  this->SetNthInput(0, const_cast< InputImageType* >(image));
}

template< class TPixelType >
void NonLocalMeansDenoisingFilter< TPixelType >::SetInputMask(MaskImageType* mask)
{
  m_Mask = mask;
}
}

#endif
