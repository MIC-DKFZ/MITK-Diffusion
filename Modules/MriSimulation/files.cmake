set(CPP_FILES

  mitkFiberfoxParameters.cpp
)

set(H_FILES
  mitkFiberfoxParameters.h

  itkFitFibersToImageFilter.h
  itkFibersFromPlanarFiguresFilter.h
  itkTractsToDWIImageFilter.h
  itkKspaceImageFilter.h
  itkDftImageFilter.h
  itkFieldmapGeneratorFilter.h
  itkRandomPhantomFilter.h

  SignalModels/mitkDiffusionSignalModel.h
  SignalModels/mitkTensorModel.h
  SignalModels/mitkBallModel.h
  SignalModels/mitkDotModel.h
  SignalModels/mitkAstroStickModel.h
  SignalModels/mitkStickModel.h
  SignalModels/mitkRawShModel.h
  SignalModels/mitkDiffusionNoiseModel.h
  SignalModels/mitkRicianNoiseModel.h
  SignalModels/mitkChiSquareNoiseModel.h

  Sequences/mitkAcquisitionType.h
  Sequences/mitkSingleShotEpi.h
  Sequences/mitkConventionalSpinEcho.h
  Sequences/mitkFastSpinEcho.h
)
