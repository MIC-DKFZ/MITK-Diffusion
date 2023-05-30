set(CPP_FILES

  mitkDiffusionServiceActivator.cpp
  mitkDiffusionIOMimeTypes.cpp

  ObjectFactories/mitkDiffusionCoreObjectFactory.cpp
  ObjectFactories/mitkConnectomicsObjectFactory.cpp

  ReaderWriter/mitkConnectomicsNetworkReader.cpp
  ReaderWriter/mitkConnectomicsNetworkWriter.cpp
  ReaderWriter/mitkConnectomicsNetworkCSVWriter.cpp
  ReaderWriter/mitkConnectomicsNetworkMatrixWriter.cpp
  ReaderWriter/mitkConnectomicsNetworkSerializer.cpp
  ReaderWriter/mitkConnectomicsNetworkDefinitions.cpp

  ReaderWriter/mitkFiberBundleDicomReader.cpp
  ReaderWriter/mitkFiberBundleDicomWriter.cpp
  ReaderWriter/mitkFiberBundleTckReader.cpp
  ReaderWriter/mitkFiberBundleTrackVisReader.cpp
  ReaderWriter/mitkFiberBundleTrackVisWriter.cpp
  ReaderWriter/mitkFiberBundleVtkReader.cpp
  ReaderWriter/mitkFiberBundleVtkWriter.cpp
  ReaderWriter/mitkFiberBundleSerializer.cpp

  ReaderWriter/mitkPlanarFigureCompositeWriter.cpp
  ReaderWriter/mitkPlanarFigureCompositeReader.cpp
  ReaderWriter/mitkPlanarFigureCompositeSerializer.cpp

  ReaderWriter/mitkPeakImageReader.cpp
  ReaderWriter/mitkPeakImageSerializer.cpp

  ReaderWriter/mitkNrrdTensorImageReader.cpp
  ReaderWriter/mitkNrrdTensorImageWriter.cpp
  ReaderWriter/mitkTensorImageSerializer.cpp

  ReaderWriter/mitkNrrdOdfImageReader.cpp
  ReaderWriter/mitkNrrdOdfImageWriter.cpp
  ReaderWriter/mitkOdfImageSerializer.cpp

  ReaderWriter/mitkShImageReader.cpp
  ReaderWriter/mitkShImageWriter.cpp
  ReaderWriter/mitkShImageSerializer.cpp

  ReaderWriter/mitkDiffusionImageDicomReader.cpp
  ReaderWriter/mitkDiffusionImageNrrdReader.cpp
  ReaderWriter/mitkDiffusionImageNrrdWriter.cpp
  ReaderWriter/mitkDiffusionImageNiftiReader.cpp
  ReaderWriter/mitkDiffusionImageNiftiWriter.cpp
)

