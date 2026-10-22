/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionConstIterator.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>

#ifdef _OPENMP
#include <omp.h>
#endif

class mitkDiffusionModellingAlgorithmsTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingAlgorithmsTestSuite);
  MITK_TEST(EigenvalueCorrection_ThreadCount_SameResult);
  CPPUNIT_TEST_SUITE_END();

  using TensorImageType = itk::Image<itk::DiffusionTensor3D<float>, 3>;

  mitk::Image::Pointer LoadDwi(const std::string &path)
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath(path), &functor);
  }

  TensorImageType::Pointer CorrectEigenvalues(mitk::Image::Pointer dwi)
  {
    auto itkDwi = itk::VectorImage<short, 3>::New();
    mitk::CastToItkImage(dwi, itkDwi);

    using FilterType = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>;
    auto filter = FilterType::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi);
    filter->SetB0Threshold(0);
    filter->Update();
    return filter->GetOutput();
  }

public:
  // The filter smooths the DWI in place, voxel after voxel, so its result must not depend on the
  // number of OpenMP threads. diffusionImage.hdwi has about 2000 voxels with negative eigenvalues,
  // so the correction loop runs several iterations.
  void EigenvalueCorrection_ThreadCount_SameResult()
  {
    auto dwi = LoadDwi("DiffusionImaging/diffusionImage.hdwi");

#ifdef _OPENMP
    const int maxThreads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    auto serial = CorrectEigenvalues(dwi);
#ifdef _OPENMP
    omp_set_num_threads(std::max(maxThreads, 4));
#endif
    auto parallel1 = CorrectEigenvalues(dwi);
    auto parallel2 = CorrectEigenvalues(dwi);
#ifdef _OPENMP
    omp_set_num_threads(maxThreads);
#endif

    unsigned int differing = 0;
    unsigned int nonZero = 0;
    itk::ImageRegionConstIterator<TensorImageType> it0(serial, serial->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> it1(parallel1, parallel1->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TensorImageType> it2(parallel2, parallel2->GetLargestPossibleRegion());
    for (; !it0.IsAtEnd(); ++it0, ++it1, ++it2)
    {
      for (unsigned int i = 0; i < 6; ++i)
      {
        if (it0.Get()[i] != it1.Get()[i] || it0.Get()[i] != it2.Get()[i])
          ++differing;
        if (it0.Get()[i] != 0)
          ++nonZero;
      }
    }
    CPPUNIT_ASSERT_MESSAGE("Corrected tensors are not all zero", nonZero > 0);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Tensor components that differ between 1 and N threads", 0u, differing);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingAlgorithms)
