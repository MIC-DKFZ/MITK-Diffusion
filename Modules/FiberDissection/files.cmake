set(H_FILES

    #Interactor
    Interactor/mitkStreamlineInteractor.h

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.h
)

set(CPP_FILES

    #Interactor
    Interactor/mitkStreamlineInteractor.cpp

    #MachineLearning
    MachineLearning/mitkStreamlineFeatureExtractor.cpp
)

set(RESOURCE_FILES
  # "Interactions" prefix forced by mitk::StateMachine
  Interactions/Streamline3DStates.xml
  Interactions/Streamline3DConfig.xml
)
