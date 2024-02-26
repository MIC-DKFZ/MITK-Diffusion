set(H_FILES

  mitkPlanarFigureComposite.h
  IO/mitkPlanarFigureCompositeWriter.h
  IO/mitkPlanarFigureCompositeReader.h
  IO/mitkPlanarFigureCompositeSerializer.h


  mitkFiberBundle.h

  IO/mitkFiberBundleObjectFactory.h
  IO/mitkFiberBundleMimeTypes.h

  IO/mitkTrackvis.h
  IO/mitkFiberBundleDicomReader.h
  IO/mitkFiberBundleDicomWriter.h
  IO/mitkFiberBundleTckReader.h
  IO/mitkFiberBundleTrackVisReader.h
  IO/mitkFiberBundleTrackVisWriter.h
  IO/mitkFiberBundleVtkReader.h
  IO/mitkFiberBundleVtkWriter.h
  IO/mitkFiberBundleSerializer.h

  Rendering/mitkFiberBundleMapper2D.h
  Rendering/mitkFiberBundleMapper3D.h
)

set(CPP_FILES

  mitkPlanarFigureComposite.cpp
  IO/mitkPlanarFigureCompositeWriter.cpp
  IO/mitkPlanarFigureCompositeReader.cpp
  IO/mitkPlanarFigureCompositeSerializer.cpp

  mitkFiberBundle.cpp

  IO/mitkFiberBundleServiceActivator.cpp
  IO/mitkFiberBundleMimeTypes.cpp
  IO/mitkFiberBundleObjectFactory.cpp

  IO/mitkTrackvis.cpp
  IO/mitkFiberBundleDicomReader.cpp
  IO/mitkFiberBundleDicomWriter.cpp
  IO/mitkFiberBundleTckReader.cpp
  IO/mitkFiberBundleTrackVisReader.cpp
  IO/mitkFiberBundleTrackVisWriter.cpp
  IO/mitkFiberBundleVtkReader.cpp
  IO/mitkFiberBundleVtkWriter.cpp
  IO/mitkFiberBundleSerializer.cpp

  Rendering/mitkFiberBundleMapper2D.cpp
  Rendering/mitkFiberBundleMapper3D.cpp
)

set(RESOURCE_FILES

)
