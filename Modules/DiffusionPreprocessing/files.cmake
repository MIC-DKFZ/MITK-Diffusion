set(H_FILES

  itkB0ImageExtractionImageFilter.h
  itkB0ImageExtractionToSeparateImageFilter.h
  itkExtractDwiChannelFilter.h

  itkElectrostaticRepulsionDiffusionGradientReductionFilter.h
  itkMergeDiffusionImagesFilter.h
  itkResampleDwiImageFilter.h
  itkDwiGradientLengthCorrectionFilter.h
  itkRemoveDwiChannelFilter.h
  itkNonLocalMeansDenoisingFilter.h
  itkVectorImageToImageFilter.h
  mitkDiffusionImageDenoising.h
  mitkDiffusionImageEditing.h
)

set(CPP_FILES

  itkB0ImageExtractionImageFilter.txx
  itkB0ImageExtractionToSeparateImageFilter.txx
  itkExtractDwiChannelFilter.txx

  itkDwiGradientLengthCorrectionFilter.cpp
  mitkDiffusionImageDenoising.cpp
  mitkDiffusionImageEditing.cpp
)

set(RESOURCE_FILES

)
