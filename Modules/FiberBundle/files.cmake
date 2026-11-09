set(H_FILES

  mitkPlanarFigureComposite.h
  IO/mitkPlanarFigureCompositeWriter.h
  IO/mitkPlanarFigureCompositeReader.h
  IO/mitkPlanarFigureCompositeSerializer.h


  mitkFiberBundle.h

  IO/mitkFiberBundleMimeTypes.h
  IO/mitkFiberBundleIOHelper.h
  IO/mitkFiberBundleDicomReferenceTags.h

  IO/mitkFiberBundleDicomReader.h
  IO/mitkFiberBundleDicomWriter.h
  IO/mitkFiberBundleTckReader.h
  IO/mitkFiberBundleTrxReader.h
  IO/mitkFiberBundleTrxWriter.h
  IO/mitkTrxFile.h
  IO/mitkFiberBundleTrackVisReader.h
  IO/mitkFiberBundleTrackVisWriter.h
  IO/mitkFiberBundleVtkReader.h
  IO/mitkFiberBundleVtkWriter.h
  IO/mitkFiberBundleSerializer.h

  Rendering/mitkFiberBundleMapper2D.h
  Rendering/mitkFiberBundleMapper3D.h
  Rendering/mitkFiberTubeShader.h
  Rendering/mitkFiberSlabIndex.h
  Rendering/mitkFiberSlabTiles.h
)

set(CPP_FILES

  mitkPlanarFigureComposite.cpp
  IO/mitkPlanarFigureCompositeWriter.cpp
  IO/mitkPlanarFigureCompositeReader.cpp
  IO/mitkPlanarFigureCompositeSerializer.cpp

  mitkFiberBundle.cpp
  mitkFiberBundleColoring.cpp
  mitkFiberBundleResampling.cpp
  mitkFiberBundleRoiExtraction.cpp
  mitkFiberBundleSelection.cpp
  mitkFiberBundleTransforms.cpp
  mitkFiberPolyDataBuilder.cpp

  IO/mitkFiberBundleServiceActivator.cpp
  IO/mitkFiberBundleMimeTypes.cpp
  IO/mitkFiberBundleIOHelper.cpp
  IO/mitkFiberBundleDicomReferenceTags.cpp

  IO/mitkTrackVisFile.cpp
  IO/mitkFiberBundleDicomReader.cpp
  IO/mitkFiberBundleDicomWriter.cpp
  IO/mitkFiberBundleTckReader.cpp
  IO/mitkFiberBundleTrxReader.cpp
  IO/mitkFiberBundleTrxWriter.cpp
  IO/mitkTrxFile.cpp
  IO/mitkFiberBundleTrackVisReader.cpp
  IO/mitkFiberBundleTrackVisWriter.cpp
  IO/mitkFiberBundleVtkReader.cpp
  IO/mitkFiberBundleVtkWriter.cpp
  IO/mitkFiberBundleSerializer.cpp

  Rendering/mitkFiberBundleMapper2D.cpp
  Rendering/mitkFiberBundleMapper3D.cpp
  Rendering/mitkFiberTubeShader.cpp
  Rendering/mitkFiberSlabIndex.cpp
  Rendering/mitkFiberSlabTiles.cpp
)

set(RESOURCE_FILES

)
