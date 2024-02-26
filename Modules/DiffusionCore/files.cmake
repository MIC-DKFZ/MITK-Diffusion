set(CPP_FILES

  mitkDiffusionFunctionCollection.cpp

  # DataStructures
  IO/mitkOdfImage.cpp
  IO/mitkShImage.cpp
  IO/mitkTensorImage.cpp
  IO/mitkPeakImage.cpp

  IO/mitkDiffusionServiceActivator.cpp
  IO/mitkDiffusionIOMimeTypes.cpp
  IO/mitkDiffusionCoreObjectFactory.cpp

  IO/mitkNrrdTensorImageReader.cpp
  IO/mitkNrrdTensorImageWriter.cpp
  IO/mitkTensorImageSerializer.cpp
  IO/mitkNrrdOdfImageReader.cpp
  IO/mitkNrrdOdfImageWriter.cpp
  IO/mitkOdfImageSerializer.cpp
  IO/mitkShImageReader.cpp
  IO/mitkShImageWriter.cpp
  IO/mitkShImageSerializer.cpp
  IO/mitkPeakImageReader.h
  IO/mitkPeakImageSerializer.h

  # Mapper etc.
  Rendering/vtkMaskedProgrammableGlyphFilter.cpp
  Rendering/mitkVectorImageVtkGlyphMapper3D.cpp
  Rendering/vtkOdfSource.cxx
  Rendering/vtkThickPlane.cxx
  Rendering/mitkOdfNormalizationMethodProperty.cpp
  Rendering/mitkOdfScaleByProperty.cpp
  Rendering/mitkCompositeMapper.cpp
  Rendering/mitkPeakImageMapper2D.cpp
  Rendering/mitkPeakImageMapper3D.cpp

  # Algorithms
  
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkADCFitFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.cpp
)

set(H_FILES

  mitkDiffusionFunctionCollection.h

  # DataStructures
  IO/mitkOdfImage.h
  IO/mitkShImage.h
  IO/mitkTensorImage.h
  IO/mitkPeakImage.h

  IO/mitkDiffusionIOMimeTypes.h
  IO/mitkDiffusionCoreObjectFactory.h

  IO/mitkNrrdTensorImageReader.h
  IO/mitkNrrdTensorImageWriter.h
  IO/mitkTensorImageSerializer.h
  IO/mitkNrrdOdfImageReader.h
  IO/mitkNrrdOdfImageWriter.h
  IO/mitkOdfImageSerializer.h
  IO/mitkShImageReader.h
  IO/mitkShImageWriter.h
  IO/mitkShImageSerializer.h
  IO/mitkPeakImageReader.cpp
  IO/mitkPeakImageSerializer.cpp

  # Rendering
  Rendering/mitkOdfVtkMapper2D.h
  Rendering/mitkCompositeMapper.h
  Rendering/mitkPeakImageMapper2D.h
  Rendering/mitkPeakImageMapper3D.h

  # Reconstruction
  Algorithms/Reconstruction/itkDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/mitkTeemDiffusionTensor3DReconstructionImageFilter.h
  Algorithms/Reconstruction/itkAnalyticalDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkDiffusionMultiShellQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkPointShell.h
  Algorithms/Reconstruction/itkOrientationDistributionFunction.h

  # MultishellProcessing
  Algorithms/Reconstruction/MultishellProcessing/itkRadialMultishellToSingleshellImageFilter.h
  Algorithms/Reconstruction/MultishellProcessing/itkDWIVoxelFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkADCFitFunctor.h

  # Algorithms
  Algorithms/itkDiffusionOdfGeneralizedFaImageFilter.h
  Algorithms/itkTensorDerivedMeasurementsFilter.h
  Algorithms/itkTensorImageToDiffusionImageFilter.h
  Algorithms/itkTensorToL2NormImageFilter.h
  Algorithms/itkGaussianInterpolateImageFunction.h
  
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

