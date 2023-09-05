set(H_FILES

    #Interactor
    Interactor/mitkSphereInteractor.h
    Interactor/mitkStreamlineInteractor.h
    Interactor/mitkStreamlineInteractorBrush.h

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.h
)

set(CPP_FILES

    #Interactor
    Interactor/mitkSphereInteractor.cpp
    Interactor/mitkStreamlineInteractor.cpp
    Interactor/mitkStreamlineInteractorBrush.cpp

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.cpp
)

set(RESOURCE_FILES
  # "Interactions" prefix forced by mitk::StateMachine

  Interactions/SphereInteractionsStates.xml
  Interactions/SphereInteractionsConfig.xml
  Interactions/Streamline3DStates.xml
  Interactions/Streamline3DConfig.xml
  Interactions/StreamlineBrush3DStates.xml
  Interactions/StreamlineBrush3DConfig.xml
)
