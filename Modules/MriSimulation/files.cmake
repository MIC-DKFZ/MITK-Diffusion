set(CPP_FILES

  mitkFiberfoxParameters.cpp

  Algorithms/itkFitFibersToImageFilter.cpp
  Algorithms/itkFibersFromPlanarFiguresFilter.cpp
  Algorithms/itkRandomPhantomFilter.cpp
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
  Algorithms/mitkFiberfoxMotion.h

  SignalModels/mitkDiffusionSignalModel.h
  SignalModels/mitkTensorModel.h
  SignalModels/mitkTensorModel.txx
  SignalModels/mitkBallModel.h
  SignalModels/mitkBallModel.txx
  SignalModels/mitkDotModel.h
  SignalModels/mitkDotModel.txx
  SignalModels/mitkAstroStickModel.h
  SignalModels/mitkAstroStickModel.txx
  SignalModels/mitkStickModel.h
  SignalModels/mitkStickModel.txx
  SignalModels/mitkRawShModel.h
  SignalModels/mitkRawShModel.txx
  SignalModels/mitkDiffusionNoiseModel.h
  SignalModels/mitkRicianNoiseModel.h
  SignalModels/mitkRicianNoiseModel.txx
  SignalModels/mitkNoncentralChiNoiseModel.h
  SignalModels/mitkNoncentralChiNoiseModel.txx

  Sequences/mitkAcquisitionType.h
  Sequences/mitkSingleShotEpi.h
  Sequences/mitkConventionalSpinEcho.h
  Sequences/mitkFastSpinEcho.h
)
