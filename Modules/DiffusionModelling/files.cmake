set(CPP_FILES

  mitkDiffusionModellingHelperFunctions.cpp

  # DataStructures
  IO/mitkRgbCachedImage.cpp
  IO/mitkOdfImage.cpp
  IO/mitkShImage.cpp
  IO/mitkTensorImage.cpp
  IO/mitkPeakImage.cpp

  IO/mitkDiffusionServiceActivator.cpp
  IO/mitkDiffusionIOMimeTypes.cpp

  IO/mitkNrrdTensorImageReader.cpp
  IO/mitkNrrdTensorImageWriter.cpp
  IO/mitkTensorImageSerializer.cpp
  IO/mitkNrrdOdfImageReader.cpp
  IO/mitkNrrdOdfImageWriter.cpp
  IO/mitkOdfImageSerializer.cpp
  IO/mitkShImageReader.cpp
  IO/mitkShImageWriter.cpp
  IO/mitkShImageSerializer.cpp
  IO/mitkPeakImageReader.cpp
  IO/mitkPeakImageSerializer.cpp

  # Mapper etc.
  Rendering/mitkOdfNormalizationMethodProperty.cpp
  Rendering/mitkOdfScaleByProperty.cpp
  Rendering/mitkCompositeMapper.cpp
  Rendering/mitkOdfGlyphGenerator.cpp
  Rendering/mitkOdfVtkMapper2D.cpp
  Rendering/mitkPeakImageMapper2D.cpp
  Rendering/mitkPeakImageMapper3D.cpp

  # Algorithms
  Algorithms/mitkLbfgsbMutex.cpp
  Algorithms/Reconstruction/mitkDiffusionReconstruction.cpp
  Algorithms/Reconstruction/mitkPointShell.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.cpp
)

set(H_FILES

  mitkDiffusionModellingHelperFunctions.h

  # DataStructures
  IO/mitkRgbCachedImage.h
  IO/mitkOdfImage.h
  IO/mitkShImage.h
  IO/mitkTensorImage.h
  IO/mitkPeakImage.h

  IO/mitkDiffusionIOMimeTypes.h

  IO/mitkNrrdTensorImageReader.h
  IO/mitkNrrdTensorImageWriter.h
  IO/mitkTensorImageSerializer.h
  IO/mitkNrrdOdfImageReader.h
  IO/mitkNrrdOdfImageWriter.h
  IO/mitkOdfImageSerializer.h
  IO/mitkShImageReader.h
  IO/mitkShImageWriter.h
  IO/mitkShImageSerializer.h
  IO/mitkPeakImageReader.h
  IO/mitkPeakImageSerializer.h

  # Rendering
  Rendering/mitkOdfGlyphGenerator.h
  Rendering/mitkOdfVtkMapper2D.h
  Rendering/mitkCompositeMapper.h
  Rendering/mitkPeakImageMapper2D.h
  Rendering/mitkPeakImageMapper3D.h

  # Reconstruction
  Algorithms/Reconstruction/mitkDiffusionReconstruction.h
  Algorithms/Reconstruction/mitkSingleShellGradients.h
  Algorithms/Reconstruction/mitkShRegularization.h
  Algorithms/Reconstruction/itkDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkAnalyticalDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkDiffusionMultiShellQballReconstructionImageFilter.h
  Algorithms/Reconstruction/mitkPointShell.h
  Algorithms/Reconstruction/itkOrientationDistributionFunction.h

  # MultishellProcessing
  Algorithms/Reconstruction/MultishellProcessing/itkRadialMultishellToSingleshellImageFilter.h
  Algorithms/Reconstruction/MultishellProcessing/itkDWIVoxelFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.h

  # Algorithms
  Algorithms/itkDiffusionOdfGeneralizedFaImageFilter.h
  Algorithms/itkTensorDerivedMeasurementsFilter.h
  Algorithms/itkTensorImageToDiffusionImageFilter.h
  Algorithms/itkTensorToL2NormImageFilter.h
  Algorithms/itkTensorImageToOdfImageFilter.h
  Algorithms/itkTensorToRgbImageFilter.h
  Algorithms/itkOdfToRgbImageFilter.h
  Algorithms/itkShToRgbImageFilter.h
  Algorithms/mitkDirectionColor.h
  Algorithms/mitkImageGeometry4D.h
  Algorithms/mitkLbfgsbMutex.h
  Algorithms/itkDiffusionTensorPrincipalDirectionImageFilter.h
  Algorithms/itkTensorReconstructionWithEigenvalueCorrectionFilter.h
  Algorithms/itkShCoefficientImageImporter.h
  Algorithms/itkShCoefficientImageExporter.h
  Algorithms/itkOdfMaximaExtractionFilter.h
  Algorithms/itkAdcImageFilter.h
  Algorithms/itkFlipPeaksFilter.h
  Algorithms/itkShToOdfImageFilter.h
)

set( TOOL_FILES
)

