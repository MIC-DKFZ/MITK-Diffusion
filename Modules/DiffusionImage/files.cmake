set(H_FILES

    mitkDiffusionImageHeaderInformation.h
    mitkDiffusionImageHelperFunctions.h
    itkB0ImageExtractionImageFilter.h
    itkB0ImageExtractionToSeparateImageFilter.h
    itkExtractDwiChannelFilter.h

    IO/mitkDiffusionImageMimeTypes.h
    IO/mitkDiffusionImageObjectFactory.h

    IO/mitkDiffusionImageDicomReader.h
    IO/mitkDiffusionImageNrrdReader.h
    IO/mitkDiffusionImageNrrdWriter.h
    IO/mitkDiffusionImageNiftiReader.h
    IO/mitkDiffusionImageNiftiWriter.h

    # DicomImport
    IO/DicomImport/mitkDicomDiffusionImageHeaderReader.h
    IO/DicomImport/mitkGEDicomDiffusionImageHeaderReader.h
    IO/DicomImport/mitkPhilipsDicomDiffusionImageHeaderReader.h
    IO/DicomImport/mitkSiemensDicomDiffusionImageHeaderReader.h
    IO/DicomImport/mitkSiemensMosaicDicomDiffusionImageHeaderReader.h
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

)

set(CPP_FILES

    mitkDiffusionImageHeaderInformation.cpp
    mitkDiffusionImageHelperFunctions.cpp
    itkB0ImageExtractionImageFilter.txx
    itkB0ImageExtractionToSeparateImageFilter.txx
    itkExtractDwiChannelFilter.txx

    IO/mitkDiffusionImageServiceActivator.cpp
    IO/mitkDiffusionImageMimeTypes.cpp
    IO/mitkDiffusionImageObjectFactory.cpp

    IO/mitkDiffusionImageDicomReader.cpp
    IO/mitkDiffusionImageNrrdReader.cpp
    IO/mitkDiffusionImageNrrdWriter.cpp
    IO/mitkDiffusionImageNiftiReader.cpp
    IO/mitkDiffusionImageNiftiWriter.cpp

    # DicomImport
    IO/DicomImport/mitkDicomDiffusionImageHeaderReader.cpp
    IO/DicomImport/mitkGEDicomDiffusionImageHeaderReader.cpp
    IO/DicomImport/mitkPhilipsDicomDiffusionImageHeaderReader.cpp
    IO/DicomImport/mitkSiemensDicomDiffusionImageHeaderReader.cpp
    IO/DicomImport/mitkSiemensMosaicDicomDiffusionImageHeaderReader.cpp
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
