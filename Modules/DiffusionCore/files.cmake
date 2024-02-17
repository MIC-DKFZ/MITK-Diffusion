set(CPP_FILES

  mitkDiffusionFunctionCollection.cpp

  # DicomImport
  # DicomImport/mitkGroupDiffusionHeadersFilter.cpp
  DicomImport/mitkDicomDiffusionImageHeaderReader.cpp
  DicomImport/mitkGEDicomDiffusionImageHeaderReader.cpp
  DicomImport/mitkPhilipsDicomDiffusionImageHeaderReader.cpp
  DicomImport/mitkSiemensDicomDiffusionImageHeaderReader.cpp
  DicomImport/mitkSiemensMosaicDicomDiffusionImageHeaderReader.cpp

  DicomImport/mitkDiffusionDICOMFileReader.cpp
  DicomImport/mitkDiffusionHeaderDICOMFileReader.cpp
  DicomImport/mitkDiffusionHeaderSiemensDICOMFileReader.cpp
  DicomImport/mitkDiffusionHeaderSiemensDICOMFileHelper.cpp
  DicomImport/mitkDiffusionHeaderSiemensMosaicDICOMFileReader.cpp
  DicomImport/mitkDiffusionHeaderGEDICOMFileReader.cpp
  DicomImport/mitkDiffusionHeaderPhilipsDICOMFileReader.cpp

  # DataStructures
  IODataStructures/DiffusionWeightedImages/mitkDiffusionImageHeaderInformation.cpp
  IODataStructures/DiffusionWeightedImages/mitkDiffusionImageCorrectionFilter.cpp
  IODataStructures/DiffusionWeightedImages/mitkDiffusionImageCreationFilter.cpp
  IODataStructures/mitkOdfImage.cpp
  IODataStructures/mitkShImage.cpp
  IODataStructures/mitkTensorImage.cpp
  IODataStructures/mitkPeakImage.cpp
  IODataStructures/mitkFiberBundle.cpp
  IODataStructures/mitkTrackvis.cpp
  IODataStructures/mitkPlanarFigureComposite.cpp


  IODataStructures/mitkStatisticModel.cpp

  # Properties
  IODataStructures/Properties/mitkBValueMapProperty.cpp
  IODataStructures/Properties/mitkGradientDirectionsProperty.cpp
  IODataStructures/Properties/mitkMeasurementFrameProperty.cpp
  IODataStructures/Properties/mitkDiffusionPropertyHelper.cpp
  IODataStructures/Properties/mitkNodePredicateIsDWI.cpp

  # Serializer
  IODataStructures/Properties/mitkBValueMapPropertySerializer.cpp
  IODataStructures/Properties/mitkGradientDirectionsPropertySerializer.cpp
  IODataStructures/Properties/mitkMeasurementFramePropertySerializer.cpp

  # Mapper etc.
  Rendering/vtkMaskedProgrammableGlyphFilter.cpp
  Rendering/mitkVectorImageVtkGlyphMapper3D.cpp
  Rendering/vtkOdfSource.cxx
  Rendering/vtkThickPlane.cxx
  Rendering/mitkOdfNormalizationMethodProperty.cpp
  Rendering/mitkOdfScaleByProperty.cpp
  Rendering/mitkFiberBundleMapper2D.cpp
  Rendering/mitkFiberBundleMapper3D.cpp
  Rendering/mitkPeakImageMapper2D.cpp
  Rendering/mitkPeakImageMapper3D.cpp
  Rendering/mitkCompositeMapper.cpp

  # Algorithms
  Algorithms/itkDwiGradientLengthCorrectionFilter.cpp
  Algorithms/Registration/mitkRegistrationWrapper.cpp
  Algorithms/Registration/mitkPyramidImageRegistrationMethod.cpp
  # Algorithms/Registration/mitkRegistrationMethodITK4.cpp
  Algorithms/Registration/mitkDWIHeadMotionCorrectionFilter.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkADCFitFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.cpp
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.cpp
)

set(H_FILES
  mitkDiffusionFunctionCollection.h

  # Rendering
  Rendering/mitkOdfVtkMapper2D.h

  # Reconstruction
  Algorithms/Reconstruction/itkDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/mitkTeemDiffusionTensor3DReconstructionImageFilter.h
  Algorithms/Reconstruction/itkAnalyticalDiffusionQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkDiffusionMultiShellQballReconstructionImageFilter.h
  Algorithms/Reconstruction/itkPointShell.h
  Algorithms/Reconstruction/itkOrientationDistributionFunction.h
  Algorithms/Reconstruction/itkBallAndSticksImageFilter.h
  Algorithms/Reconstruction/itkMultiTensorImageFilter.h

  # Fitting functions
  Algorithms/Reconstruction/FittingFunctions/mitkAbstractFitter.h
  Algorithms/Reconstruction/FittingFunctions/mitkMultiTensorFitter.h
  Algorithms/Reconstruction/FittingFunctions/mitkBallStickFitter.h


  # MultishellProcessing
  Algorithms/Reconstruction/MultishellProcessing/itkRadialMultishellToSingleshellImageFilter.h
  Algorithms/Reconstruction/MultishellProcessing/itkDWIVoxelFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkADCAverageFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkKurtosisFitFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkBiExpFitFunctor.h
  Algorithms/Reconstruction/MultishellProcessing/itkADCFitFunctor.h

  IODataStructures/mitkFiberBundle.h
  IODataStructures/mitkTrackvis.h
  IODataStructures/mitkPlanarFigureComposite.h

  IODataStructures/Properties/mitkBValueMapProperty.h
  IODataStructures/Properties/mitkGradientDirectionsProperty.h
  IODataStructures/Properties/mitkMeasurementFrameProperty.h
  IODataStructures/Properties/mitkDiffusionPropertyHelper.h

  IODataStructures/DiffusionWeightedImages/mitkDiffusionImageTransformedCreationFilter.h

  # Algorithms
  Algorithms/itkDiffusionOdfGeneralizedFaImageFilter.h
  Algorithms/itkDiffusionOdfPrepareVisualizationImageFilter.h
  Algorithms/itkElectrostaticRepulsionDiffusionGradientReductionFilter.h
  Algorithms/itkTensorDerivedMeasurementsFilter.h
  Algorithms/itkBrainMaskExtractionImageFilter.h
  Algorithms/itkB0ImageExtractionImageFilter.h
  Algorithms/itkB0ImageExtractionToSeparateImageFilter.h
  Algorithms/itkTensorImageToDiffusionImageFilter.h
  Algorithms/itkTensorToL2NormImageFilter.h
  Algorithms/itkGaussianInterpolateImageFunction.h
  
  Algorithms/itkDiffusionTensorPrincipalDirectionImageFilter.h
  Algorithms/itkCartesianToPolarVectorImageFilter.h
  Algorithms/itkPolarToCartesianVectorImageFilter.h
  Algorithms/itkExtractChannelFromRgbaImageFilter.h
  Algorithms/itkTensorReconstructionWithEigenvalueCorrectionFilter.h
  Algorithms/itkMergeDiffusionImagesFilter.h
  Algorithms/itkShCoefficientImageImporter.h
  Algorithms/itkShCoefficientImageExporter.h
  Algorithms/itkOdfMaximaExtractionFilter.h
  Algorithms/itkResampleDwiImageFilter.h
  Algorithms/itkDwiGradientLengthCorrectionFilter.h
  Algorithms/itkAdcImageFilter.h
  Algorithms/itkDwiNormilzationFilter.h
  Algorithms/itkSplitDWImageFilter.h
  Algorithms/itkRemoveDwiChannelFilter.h
  Algorithms/itkExtractDwiChannelFilter.h
  Algorithms/itkFlipPeaksFilter.h
  Algorithms/itkShToOdfImageFilter.h
  Algorithms/itkFourDToVectorImageFilter.h
  Algorithms/itkVectorImageToFourDImageFilter.h
  Algorithms/Registration/mitkDWIHeadMotionCorrectionFilter.h
  Algorithms/itkNonLocalMeansDenoisingFilter.h
  Algorithms/itkVectorImageToImageFilter.h
  Algorithms/itkSplitVectorImageFilter.h
  Algorithms/itkRegularizedIVIMLocalVariationImageFilter.h
  Algorithms/itkRegularizedIVIMReconstructionFilter.h
  Algorithms/itkRegularizedIVIMReconstructionSingleIteration.h
)

set( TOOL_FILES
)

