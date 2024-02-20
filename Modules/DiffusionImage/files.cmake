set(H_FILES

    IO/mitkDiffusionImageMimeTypes.h
    IO/mitkDiffusionImageObjectFactory.h

    IO/mitkDiffusionImageDicomReader.h
    IO/mitkDiffusionImageNrrdReader.h
    IO/mitkDiffusionImageNrrdWriter.h
    IO/mitkDiffusionImageNiftiReader.h
    IO/mitkDiffusionImageNiftiWriter.h

    Properties/mitkBValueMapProperty.h
    Properties/mitkGradientDirectionsProperty.h
    Properties/mitkMeasurementFrameProperty.h
    Properties/mitkDiffusionPropertyHelper.h
)

set(CPP_FILES

    IO/mitkDiffusionImageServiceActivator.cpp
    IO/mitkDiffusionImageMimeTypes.cpp
    IO/mitkDiffusionImageObjectFactory.cpp

    IO/mitkDiffusionImageDicomReader.cpp
    IO/mitkDiffusionImageNrrdReader.cpp
    IO/mitkDiffusionImageNrrdWriter.cpp
    IO/mitkDiffusionImageNiftiReader.cpp
    IO/mitkDiffusionImageNiftiWriter.cpp


    # Properties
    Properties/mitkBValueMapProperty.cpp
    Properties/mitkGradientDirectionsProperty.cpp
    Properties/mitkMeasurementFrameProperty.cpp
    Properties/mitkDiffusionPropertyHelper.cpp
    Properties/mitkNodePredicateIsDWI.cpp

    # Serializer
    Properties/mitkBValueMapPropertySerializer.cpp
    Properties/mitkGradientDirectionsPropertySerializer.cpp
    Properties/mitkMeasurementFramePropertySerializer.cpp
)

set(RESOURCE_FILES

)
