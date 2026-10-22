set(H_FILES

    mitkDiffusionImageHelperFunctions.h

    IO/mitkDiffusionImageMimeTypes.h

    IO/mitkDiffusionImageDicomReader.h
    IO/mitkDiffusionImageNrrdReader.h
    IO/mitkDiffusionImageNrrdWriter.h
    IO/mitkDiffusionImageNiftiReader.h
    IO/mitkDiffusionImageNiftiWriter.h

    # DicomImport
    IO/DicomImport/mitkDiffusionDICOMFileReader.h
    IO/DicomImport/mitkDiffusionHeaderDICOMFileReader.h
    IO/DicomImport/mitkDiffusionHeaderSiemensDICOMFileReader.h
    IO/DicomImport/mitkDiffusionHeaderSiemensDICOMFileHelper.h
    IO/DicomImport/mitkDiffusionHeaderSiemensMosaicDICOMFileReader.h
    IO/DicomImport/mitkDiffusionHeaderGEDICOMFileReader.h
    IO/DicomImport/mitkDiffusionHeaderPhilipsDICOMFileReader.h

    # Properties
    Properties/mitkBValueMapProperty.h
    Properties/mitkGradientDirectionsProperty.h
    Properties/mitkMeasurementFrameProperty.h
    Properties/mitkDiffusionPropertyHelper.h
    Properties/mitkNodePredicateIsDWI.h

)

set(CPP_FILES

    mitkDiffusionImageHelperFunctions.cpp

    IO/mitkDiffusionImageServiceActivator.cpp
    IO/mitkDiffusionImageMimeTypes.cpp

    IO/mitkDiffusionImageDicomReader.cpp
    IO/mitkDiffusionImageNrrdReader.cpp
    IO/mitkDiffusionImageNrrdWriter.cpp
    IO/mitkDiffusionImageNiftiReader.cpp
    IO/mitkDiffusionImageNiftiWriter.cpp

    # DicomImport
    IO/DicomImport/mitkDiffusionDICOMFileReader.cpp
    IO/DicomImport/mitkDiffusionHeaderDICOMFileReader.cpp
    IO/DicomImport/mitkDiffusionHeaderSiemensDICOMFileReader.cpp
    IO/DicomImport/mitkDiffusionHeaderSiemensDICOMFileHelper.cpp
    IO/DicomImport/mitkDiffusionHeaderSiemensMosaicDICOMFileReader.cpp
    IO/DicomImport/mitkDiffusionHeaderGEDICOMFileReader.cpp
    IO/DicomImport/mitkDiffusionHeaderPhilipsDICOMFileReader.cpp

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
