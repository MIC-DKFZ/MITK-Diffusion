set(H_FILES

  itkB0ImageExtractionImageFilter.h
  itkB0ImageExtractionToSeparateImageFilter.h

  itkElectrostaticRepulsionDiffusionGradientReductionFilter.h
  itkMergeDiffusionImagesFilter.h
  itkResampleDwiImageFilter.h
  itkRemoveDwiChannelFilter.h
  itkNonLocalMeansDenoisingFilter.h
  mitkB0Indices.h
  mitkDiffusionImageDenoising.h
  mitkDiffusionImageEditing.h
)

set(CPP_FILES

  itkB0ImageExtractionImageFilter.txx
  itkB0ImageExtractionToSeparateImageFilter.txx

  mitkDiffusionImageDenoising.cpp
  mitkDiffusionImageEditing.cpp
)

set(RESOURCE_FILES

)
