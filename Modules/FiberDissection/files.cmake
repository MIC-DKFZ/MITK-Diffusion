set(H_FILES

    #Interactor
    Interactor/mitkSphereInteractor.h
    Interactor/mitkSphereSegmentationInteractor.h
    Interactor/mitkStreamlineInteractorBrush.h

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.h
)

set(CPP_FILES

    #Interactor
    Interactor/mitkSphereInteractor.cpp
    Interactor/mitkSphereSegmentationInteractor.cpp
    Interactor/mitkStreamlineInteractorBrush.cpp

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.cpp
)

set(RESOURCE_FILES
  # "Interactions" prefix forced by mitk::StateMachine

  Interactions/SphereSegmentationInteractionsStates.xml
  Interactions/SphereSegmentationInteractionsConfig.xml
  Interactions/SphereInteractionsStates.xml
  Interactions/SphereInteractionsConfig.xml
  Interactions/StreamlineBrush3DStates.xml
  Interactions/StreamlineBrush3DConfig.xml
)
