set(CPP_FILES

  mitkFiberfoxParameters.cpp
)

set(H_FILES
  mitkFiberfoxParameters.h

  Algorithms/itkFitFibersToImageFilter.h
  Algorithms/itkFibersFromPlanarFiguresFilter.h
  Algorithms/itkTractsToDWIImageFilter.h
  Algorithms/itkKspaceImageFilter.h
  Algorithms/itkDftImageFilter.h
  Algorithms/itkFieldmapGeneratorFilter.h
  Algorithms/itkRandomPhantomFilter.h

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
