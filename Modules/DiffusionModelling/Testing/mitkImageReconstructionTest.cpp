/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkImageCast.h>
#include <mitkTensorImage.h>
#include <mitkOdfImage.h>
#include <mitkIOUtil.h>
#include <itkDiffusionTensor3D.h>
#include <mitkTestingMacros.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkDiffusionTensor3D.h>
#include <itkDiffusionQballReconstructionImageFilter.h>
#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <mitkImage.h>
#include <mitkDiffusionPropertyHelper.h>


#include <mitkTestFixture.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

class mitkImageReconstructionTestSuite : public mitk::TestFixture
{

  CPPUNIT_TEST_SUITE(mitkImageReconstructionTestSuite);
  MITK_TEST(DtiTest);
  MITK_TEST(NumericQball);
  MITK_TEST(StandardQball);
  MITK_TEST(CsaQball);
  MITK_TEST(ADC);
  MITK_TEST(RawSignal);
  CPPUNIT_TEST_SUITE_END();

  typedef itk::Image<float, 3> ItkFloatImgType;

  private:

  /** Members used inside the different (sub-)tests. All members are initialized via setUp().*/
      mitk::Image::Pointer dwi;
      float b_value;
      itk::VectorImage<short,3>::Pointer itkVectorImagePointer;
      mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer gradients;

  public:

  void setUp() override
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"), &functor);
    b_value = mitk::DiffusionPropertyHelper::GetReferenceBValue( dwi );
    gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(dwi);

    itkVectorImagePointer = itk::VectorImage<short,3>::New();
    mitk::CastToItkImage(dwi, itkVectorImagePointer);
  }

  void tearDown() override
  {

  }

  void DtiTest()
  {
    auto tensorImage = mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_dti.dti"));

    typedef itk::DiffusionTensor3DReconstructionImageFilter< short, short, float > TensorReconstructionImageFilterType;
    TensorReconstructionImageFilterType::Pointer filter = TensorReconstructionImageFilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( dynamic_cast<mitk::GradientDirectionsProperty *>(dwi->GetProperty(mitk::DiffusionPropertyHelper::GetGradientContainerPropertyName().c_str()).GetPointer())->GetGradientDirectionsContainerCopy(), itkVectorImagePointer );
    filter->Update();
    mitk::TensorImage::Pointer testImage = mitk::TensorImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *tensorImage, 0.0001, true), "tensor reconstruction test.");
  }

  void NumericQball()
  {
    MITK_INFO << "Numerical Q-ball reconstruction";
    mitk::OdfImage::Pointer odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QN0.qbi"));
    typedef itk::DiffusionQballReconstructionImageFilter<short, short, float, ODF_SAMPLING_SIZE> QballReconstructionImageFilterType;
    QballReconstructionImageFilterType::Pointer filter = QballReconstructionImageFilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( gradients, itkVectorImagePointer );
    filter->SetNormalizationMethod(QballReconstructionImageFilterType::QBR_STANDARD);
    filter->Update();
    mitk::OdfImage::Pointer testImage = mitk::OdfImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *odfImage, 0.0001, true), "Numerical Q-ball reconstruction test.");
  }

  void StandardQball()
  {
    MITK_INFO << "Standard Q-ball reconstruction";
    mitk::OdfImage::Pointer odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA0.qbi"));
    typedef itk::AnalyticalDiffusionQballReconstructionImageFilter<short,short,float,4,ODF_SAMPLING_SIZE> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( gradients, itkVectorImagePointer );
    filter->SetLambda(0.006);
    filter->SetNormalizationMethod(FilterType::QBAR_STANDARD);
    filter->Update();
    mitk::OdfImage::Pointer testImage = mitk::OdfImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *odfImage, 0.0001, true), "Standard Q-ball reconstruction test.");
  }

  void CsaQball()
  {
    MITK_INFO << "CSA Q-ball reconstruction";
    mitk::OdfImage::Pointer odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA6.qbi"));
    typedef itk::AnalyticalDiffusionQballReconstructionImageFilter<short,short,float,4,ODF_SAMPLING_SIZE> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( gradients, itkVectorImagePointer );
    filter->SetLambda(0.006);
    filter->SetNormalizationMethod(FilterType::QBAR_SOLID_ANGLE);
    filter->Update();
    mitk::OdfImage::Pointer testImage = mitk::OdfImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *odfImage, 0.0001, true), "CSA Q-ball reconstruction test.");
  }

  void ADC()
  {
    MITK_INFO << "ADC profile reconstruction";
    mitk::OdfImage::Pointer odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA4.qbi"));
    typedef itk::AnalyticalDiffusionQballReconstructionImageFilter<short,short,float,4,ODF_SAMPLING_SIZE> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( gradients, itkVectorImagePointer );
    filter->SetLambda(0.006);
    filter->SetNormalizationMethod(FilterType::QBAR_ADC_ONLY);
    filter->Update();
    mitk::OdfImage::Pointer testImage = mitk::OdfImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *odfImage, 0.0001, true), "ADC profile reconstruction test.");
  }

  void RawSignal()
  {
    MITK_INFO << "Raw signal modeling";
    mitk::OdfImage::Pointer odfImage = mitk::IOUtil::Load<mitk::OdfImage>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_QA5.qbi"));
    typedef itk::AnalyticalDiffusionQballReconstructionImageFilter<short,short,float,4,ODF_SAMPLING_SIZE> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetBValue( b_value );
    filter->SetGradientImage( gradients, itkVectorImagePointer );
    filter->SetLambda(0.006);
    filter->SetNormalizationMethod(FilterType::QBAR_RAW_SIGNAL);
    filter->Update();
    mitk::OdfImage::Pointer testImage = mitk::OdfImage::New();
    testImage->InitializeByItk( filter->GetOutput() );
    testImage->SetVolume( filter->GetOutput()->GetBufferPointer() );
    MITK_TEST_CONDITION_REQUIRED(mitk::Equal(*testImage, *odfImage, 0.1, true), "Raw signal modeling test.");
  }

};

MITK_TEST_SUITE_REGISTRATION(mitkImageReconstruction)
