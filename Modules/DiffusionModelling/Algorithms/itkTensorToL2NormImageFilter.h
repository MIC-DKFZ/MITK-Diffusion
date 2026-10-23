/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef _itk_TensorToL2NormImageFilter_h_
#define _itk_TensorToL2NormImageFilter_h_

#include "itkImageToImageFilter.h"

namespace itk
{

  template <class TInputImage, class TOutputImage>
    class TensorToL2NormImageFilter
    : public ImageToImageFilter<TInputImage, TOutputImage>
  {

  public:

    typedef TensorToL2NormImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkOverrideGetNameOfClassMacro(TensorToL2NormImageFilter);

    itkStaticConstMacro (ImageDimension, unsigned int,
                         TOutputImage::ImageDimension);

    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    typedef TInputImage                           InputImageType;
    typedef typename InputImageType::PixelType    InputPixelType;
    typedef TOutputImage                          OutputImageType;
    typedef typename OutputImageType::PixelType   OutputPixelType;
    typedef typename OutputImageType::RegionType  OutputImageRegionType;

  protected:
    TensorToL2NormImageFilter(){
      
    };
    ~TensorToL2NormImageFilter() override = default;

    void PrintSelf (std::ostream& os, Indent indent) const override
    {
      Superclass::PrintSelf (os, indent);
    }


    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override;


  private:
    ITK_DISALLOW_COPY_AND_MOVE(TensorToL2NormImageFilter);

    };


} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorToL2NormImageFilter.txx"
#endif


#endif
