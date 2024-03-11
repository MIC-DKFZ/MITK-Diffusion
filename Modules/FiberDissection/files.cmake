set(H_FILES

    #Interactor
    Interactor/mitkSphereInteractor.h
    Interactor/mitkSphereSegmentationInteractor.h
    Interactor/mitkStreamlineInteractor.h
    Interactor/mitkStreamlineInteractorBrush.h

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.h
    
    MachineLearning/mitkStatisticModel.h
)

set(CPP_FILES

    #Interactor
    Interactor/mitkSphereInteractor.cpp
    Interactor/mitkSphereSegmentationInteractor.cpp
    Interactor/mitkStreamlineInteractor.cpp
    Interactor/mitkStreamlineInteractorBrush.cpp

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.cpp
    
    MachineLearning/mitkStatisticModel.cpp
)

set(RESOURCE_FILES
  # "Interactions" prefix forced by mitk::StateMachine

  Interactions/SphereSegmentationInteractionsStates.xml
  Interactions/SphereSegmentationInteractionsConfig.xml
  Interactions/SphereInteractionsStates.xml
  Interactions/SphereInteractionsConfig.xml
  Interactions/Streamline3DStates.xml
  Interactions/Streamline3DConfig.xml
  Interactions/StreamlineBrush3DStates.xml
  Interactions/StreamlineBrush3DConfig.xml
)
