/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkTestingMacros.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <itkTractsToDWIImageFilter.h>
#include <mitkFiberfoxParameters.h>
#include <mitkStickModel.h>
#include <mitkTensorModel.h>
#include <mitkBallModel.h>
#include <mitkDotModel.h>
#include <mitkAstroStickModel.h>
#include <mitkImage.h>
#include <itkTestingComparisonImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <mitkRicianNoiseModel.h>
#include <mitkIOUtil.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkProperties.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <itkVectorImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include "mitkDiffusionTestData.h"
#include "mitkTestFixture.h"

class mitkFiberfoxSignalGenerationBrainSliceTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkFiberfoxSignalGenerationBrainSliceTestSuite);
  MITK_TEST(Test1);
  MITK_TEST(Test2);
  MITK_TEST(Test3);
  MITK_TEST(Test4);
  MITK_TEST(Test5);
  MITK_TEST(Test4Exact);
  CPPUNIT_TEST_SUITE_END();

  typedef itk::VectorImage< short, 3>   ItkDwiType;

private:

public:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/
  mitk::FiberBundle::Pointer m_FiberBundle;
  mitk::Image::Pointer m_Template;

  void setUp() override
  {
    m_FiberBundle = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("Fiberfox/SignalGen_BrainSlice.fib"));
    m_Template = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/SignalGen_BrainSliceTemplate.nii.gz"));
  }

  void adjust_to_template(mitk::FiberfoxParameters& parameters)
  {
    typedef itk::Image< short, 3 >    ItkImageType;
    ItkImageType::Pointer itkImagePointer = ItkImageType::New();
    mitk::CastToItkImage(m_Template, itkImagePointer);
    parameters.m_SignalGen.m_ImageRegion = itkImagePointer->GetLargestPossibleRegion();
    parameters.m_SignalGen.m_ImageSpacing = itkImagePointer->GetSpacing();
    parameters.m_SignalGen.m_ImageOrigin = itkImagePointer->GetOrigin();
    parameters.m_SignalGen.m_ImageDirection = itkImagePointer->GetDirection();
  }

  void tearDown() override
  {

  }

  /** True if no value differs by more than maxDifference. */
  bool CompareDwi(itk::VectorImage< short, 3 >* dwi1, itk::VectorImage< short, 3 >* dwi2, short maxDifference)
  {
    bool out = true;
    typedef itk::VectorImage< short, 3 > DwiImageType;
    // the iterators below walk both images in step
    if (dwi1->GetLargestPossibleRegion().GetSize() != dwi2->GetLargestPossibleRegion().GetSize()
        || dwi1->GetVectorLength() != dwi2->GetVectorLength())
    {
      MITK_INFO << "Image size " << dwi1->GetLargestPossibleRegion().GetSize() << " x " << dwi1->GetVectorLength()
                << " differs from the reference size " << dwi2->GetLargestPossibleRegion().GetSize() << " x " << dwi2->GetVectorLength();
      return false;
    }
    try{
      itk::ImageRegionIterator< DwiImageType > it1(dwi1, dwi1->GetLargestPossibleRegion());
      itk::ImageRegionIterator< DwiImageType > it2(dwi2, dwi2->GetLargestPossibleRegion());
      int count = 0;
      while(!it1.IsAtEnd())
      {
        for (unsigned int i=0; i<dwi1->GetVectorLength(); ++i)
        {
          short d = abs(it1.Get()[i]-it2.Get()[i]);

          if (d>maxDifference)
          {
            if (count<10)
            {
              MITK_INFO << "Index: " << it1.GetIndex() << " Gradient: " << i << " - " << it1.Get()[i] << " should be " << it2.Get()[i];
            }
            out = false;
            count++;
          }
        }
        ++it1;
        ++it2;
      }
      if (count>=10)
        MITK_INFO << "Skipping errors.";
      MITK_INFO << "Errors detected: " << count;
    }
    catch(...)
    {
      return false;
    }
    return out;
  }

  /** maxDifference 1 allows rounding of the int16 output to differ; 0 asks for the reference bit for bit. */
  void StartSimulation(mitk::FiberfoxParameters parameters, mitk::Image::Pointer refImage, std::string out, short maxDifference = 1)
  {
    itk::TractsToDWIImageFilter< short >::Pointer tractsToDwiFilter = itk::TractsToDWIImageFilter< short >::New();
    tractsToDwiFilter->SetUseConstantRandSeed(true);
    tractsToDwiFilter->SetParameters(parameters);
    tractsToDwiFilter->SetFiberBundle(m_FiberBundle);
    tractsToDwiFilter->Update();

    mitk::Image::Pointer testImage = mitk::GrabItkImageMemory( tractsToDwiFilter->GetOutput() );
    mitk::DiffusionPropertyHelper::SetGradientContainer(testImage, parameters.m_SignalGen.GetItkGradientContainer());
    mitk::DiffusionPropertyHelper::SetReferenceBValue(testImage, parameters.m_SignalGen.GetBvalue());

    mitk::DiffusionPropertyHelper::InitializeImage( testImage );

    if (refImage.IsNotNull())
    {
      CPPUNIT_ASSERT_MESSAGE("Reference should be a diffusion-weighted image", mitk::DiffusionPropertyHelper::GetGradientContainer(refImage).IsNotNull());
      {
        ItkDwiType::Pointer itkTestImagePointer = ItkDwiType::New();
        mitk::CastToItkImage(testImage, itkTestImagePointer);
        ItkDwiType::Pointer itkRefImagePointer = ItkDwiType::New();
        mitk::CastToItkImage(refImage, itkRefImagePointer);

        bool cond = CompareDwi(itkTestImagePointer, itkRefImagePointer, maxDifference);
        if (!cond)
        {
          MITK_INFO << "Saving test image to " << mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + out;
          mitk::IOUtil::Save(testImage, "DWI_NIFTI", mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + out);
          mitk::IOUtil::Save(refImage, "DWI_NIFTI", mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + "REF_" + out);
        }
        CPPUNIT_ASSERT_MESSAGE("Simulated images should be equal", cond);
      }
    }
  }

  void Test1()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_1.ffp"), true);
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());

    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_1_OUT_v2.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_1_OUT.nii.gz");
  }

  void Test2()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_2.ffp"), true);
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_2_OUT_v3.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_2_OUT.nii.gz");
  }

  void Test3()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_3.ffp"), true);
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_3_OUT_v2.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_3_OUT.nii.gz");
  }

  void Test4()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_4.ffp"), true);
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_4_OUT_v4.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_4_OUT.nii.gz");
  }

  void Test5()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_5_v2.ffp"), true);
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_5_OUT_v4.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_5_OUT.nii.gz");
  }

  /**
   * The reference of BrainSlice_4 (eddy currents) was simulated with the exact sum of the k-space
   * phase, which is not the default (Test4 simulates with the time segments). The exact sum
   * reproduces it within 1, as before #47: 145 values differ by 1 from a reference older than the
   * current code, while the exact path itself gives the same image bit for bit as before #47.
   * mitkFiberfoxSignalGenerationTest checks the exact sum bit for bit (param1, param6, param8).
   */
  void Test4Exact()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_4.ffp"), true);
    parameters.m_SignalGen.m_KspaceTimeSegments = 0;
    adjust_to_template(parameters);
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer refImage = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_4_OUT_v4.nii.gz"), &functor);
    StartSimulation(parameters, refImage, "BrainSlice_4_OUT_exact.nii.gz");
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxSignalGenerationBrainSlice)
