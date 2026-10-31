/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkIOUtil.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkITKImageImport.h>
#include <mitkOdfGlyphGenerator.h>
#include <mitkOdfImage.h>
#include <mitkOdfNormalizationMethodProperty.h>
#include <mitkOdfScaleByProperty.h>
#include <mitkOdfVtkMapper2D.h>
#include <mitkProperties.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkShImage.h>
#include <mitkSlicedGeometry3D.h>
#include <mitkTensorImage.h>
#include <mitkDiffusionTestData.h>
#include <mitkOmpThreadGuard.h>
#include <mitkPointShell.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>

#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

/**
 * Rendering classes of DiffusionModelling that work without a render window: the glyph
 * generator and the display properties.
 */
class mitkDiffusionModellingRenderingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingRenderingTestSuite);
  MITK_TEST(OdfProperties_CloneKeepsType);
  MITK_TEST(OdfPropertySerializers_MissingValue_Null);
  MITK_TEST(GlyphGenerator_SliceOfStandardPlanes);
  MITK_TEST(GlyphGenerator_PlaneOutsideImage_NoSlice);
  MITK_TEST(GlyphGenerator_RotatedImage_OrientationFromWorldAxes);
  MITK_TEST(GlyphGenerator_PositionsAndWindow);
  MITK_TEST(GlyphGenerator_MaxGlyphs_SpreadsOverSlice);
  MITK_TEST(GlyphGenerator_RandomPlacement_StableSubset);
  MITK_TEST(GlyphGenerator_NormalizationScaleAndColours);
  MITK_TEST(GlyphGenerator_ShMatchesOdfOfBothConventions);
  MITK_TEST(GlyphGenerator_Order2ShIsNoTensor);
  MITK_TEST(GlyphGenerator_TensorOdfs_MatchOdfImage);
  MITK_TEST(GlyphGenerator_TensorEllipsoid_Semiaxes);
  MITK_TEST(GlyphGenerator_UnsupportedImage_Throws);
  MITK_TEST(GlyphGenerator_FullMesh_RecordedOutput);
  MITK_TEST(GlyphGenerator_ValueColours_FromLookupTable);
  MITK_TEST(GlyphGenerator_ThreadCount_SameOutput);
  MITK_TEST(GlyphGenerator_MeshLevel_ByGlyphWidth);
  MITK_TEST(GlyphGenerator_CoarseMesh_ClosedSurface);
  MITK_TEST(GlyphGenerator_CoarseMesh_TensorValuesAtItsDirections);
  MITK_TEST(GlyphGenerator_CoarseMesh_ShValuesAtItsDirections);
  MITK_TEST(GlyphGenerator_CoarseMesh_OdfValuesOfNearestDirection);
  MITK_TEST(GlyphGenerator_CoarseMesh_NoLargerThanFullGlyph);
  MITK_TEST(GlyphGenerator_EllipsoidInstances_MatchTheTensors);
  MITK_TEST(GlyphGenerator_EllipsoidInstances_OnlyForTheirSettings);
  MITK_TEST(GlyphGenerator_InstanceSphere_MeshOfTheLevel);
  MITK_TEST(OdfMapper_DefaultProperties_KeepExistingValues);
  MITK_TEST(OdfMapper_GenerationRequired);
  CPPUNIT_TEST_SUITE_END();

public:
  // Clone() returned nullptr; cloning a property list made the mapper fall back to other settings.
  void OdfProperties_CloneKeepsType()
  {
    auto normalization = mitk::OdfNormalizationMethodProperty::New(mitk::ODFN_NONE);
    auto normalizationClone = normalization->Clone();
    CPPUNIT_ASSERT(normalizationClone.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(int(mitk::ODFN_NONE), int(normalizationClone->GetNormalization()));
    mitk::BaseProperty::Pointer baseClone = static_cast<mitk::BaseProperty *>(normalization.GetPointer())->Clone();
    CPPUNIT_ASSERT(dynamic_cast<mitk::OdfNormalizationMethodProperty *>(baseClone.GetPointer()) != nullptr);

    auto scaleBy = mitk::OdfScaleByProperty::New(mitk::ODFSB_GFA);
    auto scaleByClone = scaleBy->Clone();
    CPPUNIT_ASSERT(scaleByClone.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(int(mitk::ODFSB_GFA), int(scaleByClone->GetScaleBy()));
  }

  // A scene element without "value" built a std::string from nullptr.
  void OdfPropertySerializers_MissingValue_Null()
  {
    tinyxml2::XMLDocument document;
    auto element = document.NewElement("property");
    CPPUNIT_ASSERT(mitk::OdfNormalizationMethodPropertySerializer::New()->Deserialize(element).IsNull());
    CPPUNIT_ASSERT(mitk::OdfScaleByPropertySerializer::New()->Deserialize(element).IsNull());

    element->SetAttribute("value", "Maximum");
    auto property = mitk::OdfNormalizationMethodPropertySerializer::New()->Deserialize(element);
    CPPUNIT_ASSERT(property.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(std::string("Maximum"), property->GetValueAsString());
  }

  using OdfVectorImageType = itk::Image<itk::Vector<float, ODF_SAMPLING_SIZE>, 3>;
  using Orientation = mitk::OdfGlyphGenerator::Orientation;

  static mitk::PlaneGeometry::Pointer Plane(double x, double y, double z, double nx, double ny, double nz)
  {
    mitk::Point3D origin;
    origin[0] = x;
    origin[1] = y;
    origin[2] = z;
    mitk::Vector3D normal;
    normal[0] = nx;
    normal[1] = ny;
    normal[2] = nz;
    auto plane = mitk::PlaneGeometry::New();
    plane->InitializePlane(origin, normal);
    return plane;
  }

  /** ODF image of the given size, spacing (1, 2, 3) and origin (10, 20, 30), every ODF 'value'. */
  static mitk::OdfImage::Pointer MakeOdfImage(unsigned int sx, unsigned int sy, unsigned int sz, float value)
  {
    auto itkImage = OdfVectorImageType::New();
    OdfVectorImageType::RegionType region;
    region.SetSize({{sx, sy, sz}});
    itkImage->SetRegions(region);
    const double spacing[3] = {1, 2, 3};
    const double origin[3] = {10, 20, 30};
    itkImage->SetSpacing(spacing);
    itkImage->SetOrigin(origin);
    itkImage->Allocate();
    itk::Vector<float, ODF_SAMPLING_SIZE> odf;
    odf.Fill(value);
    itkImage->FillBuffer(odf);
    auto image = mitk::OdfImage::New();
    image->InitializeByItk(itkImage.GetPointer());
    image->SetVolume(itkImage->GetBufferPointer());
    return image;
  }

  static mitk::OdfGlyphGenerator::Settings AllGlyphs()
  {
    mitk::OdfGlyphGenerator::Settings settings;
    settings.maxGlyphs = 1000000;
    settings.randomPlacement = false;
    return settings;
  }

  static mitk::Point3D GlyphCenter(vtkPolyData *glyphs, vtkIdType glyph)
  {
    const vtkIdType n = ODF_SAMPLING_SIZE;
    mitk::Point3D center;
    center.Fill(0);
    for (vtkIdType i = glyph * n; i < (glyph + 1) * n; ++i)
    {
      for (unsigned int c = 0; c < 3; ++c)
        center[c] += glyphs->GetPoint(i)[c] / n;
    }
    return center;
  }

  static double Distance(const double *a, const mitk::Point3D &b)
  {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
  }

  // Standard planes through voxel centres of tensorImage.dti (64x64x3, spacing 3, origin 0), as
  // the slice navigation creates them.
  void GlyphGenerator_SliceOfStandardPlanes()
  {
    auto image = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti"));
    const std::tuple<mitk::AnatomicalPlane, unsigned int, Orientation, int> expected[] = {
      {mitk::AnatomicalPlane::Axial, 2, Orientation::Axial, 1},
      {mitk::AnatomicalPlane::Sagittal, 0, Orientation::Sagittal, 20},
      {mitk::AnatomicalPlane::Coronal, 1, Orientation::Coronal, 63}};
    for (const auto &[anatomicalPlane, normalAxis, orientation, sliceNumber] : expected)
    {
      auto planes = mitk::SlicedGeometry3D::New();
      planes->InitializePlanes(image->GetGeometry(), anatomicalPlane);
      const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *planes->GetPlaneGeometry(sliceNumber));
      CPPUNIT_ASSERT(slice.has_value());
      CPPUNIT_ASSERT_EQUAL(normalAxis, slice->normalAxis);
      CPPUNIT_ASSERT(orientation == slice->orientation);
      CPPUNIT_ASSERT_EQUAL(3.0, slice->glyphSpacing);
      // Planes are numbered along their normal, which points against the image axis for some.
      const auto index = slice->index;
      CPPUNIT_ASSERT(index == sliceNumber || index == static_cast<int>(image->GetDimension(normalAxis)) - 1 - sliceNumber);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(static_cast<double>(index), slice->planeIndex, 1e-6);
    }
  }

  // Outside the image the former mapper showed the border slice or stale glyphs.
  void GlyphGenerator_PlaneOutsideImage_NoSlice()
  {
    auto image = MakeOdfImage(3, 2, 4, 1);
    // z of voxel layer k: 30 + 3k; the image covers k in [-0.5, 3.5].
    const std::pair<double, std::optional<int>> cases[] = {
      {30 - 3 * 0.6, std::nullopt}, {30 - 3 * 0.4, 0}, {30 + 3 * 1.4, 1}, {30 + 3 * 1.6, 2},
      {30 + 3 * 3.4, 3}, {30 + 3 * 3.6, std::nullopt}, {30 - 3 * 10.0, std::nullopt}, {30 + 3 * 10.0, std::nullopt}};
    for (const auto &[z, expectedIndex] : cases)
    {
      const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, z, 0, 0, 1));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(std::to_string(z), expectedIndex.has_value(), slice.has_value());
      if (expectedIndex)
      {
        CPPUNIT_ASSERT_EQUAL(static_cast<itk::IndexValueType>(*expectedIndex), slice->index);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((z - 30) / 3, slice->planeIndex, 1e-9);
      }
    }

    // Planes not parallel to two image axes show no glyphs.
    CPPUNIT_ASSERT(!mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(10, 20, 33, 0, 1, 1)));
  }

  // The orientation and the glyph size come from the image axis perpendicular to the plane, not
  // from the name of the render window.
  void GlyphGenerator_RotatedImage_OrientationFromWorldAxes()
  {
    auto image = MakeOdfImage(3, 2, 4, 1);
    // Image axis 0 along world z, axis 1 along world x, axis 2 along world y.
    mitk::Matrix3D matrix;
    matrix.Fill(0);
    matrix[2][0] = 1;
    matrix[0][1] = 2;
    matrix[1][2] = 3;
    auto transform = mitk::AffineTransform3D::New();
    transform->SetMatrix(matrix);
    transform->SetOffset(image->GetGeometry()->GetIndexToWorldTransform()->GetOffset());
    image->GetGeometry()->SetIndexToWorldTransform(transform);

    const auto axial = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 31, 0, 0, 1));
    CPPUNIT_ASSERT(axial.has_value());
    CPPUNIT_ASSERT_EQUAL(0u, axial->normalAxis);
    CPPUNIT_ASSERT(Orientation::Axial == axial->orientation);
    CPPUNIT_ASSERT_EQUAL(2.0, axial->glyphSpacing);

    const auto sagittal = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(12, 0, 0, 1, 0, 0));
    CPPUNIT_ASSERT(sagittal.has_value());
    CPPUNIT_ASSERT_EQUAL(1u, sagittal->normalAxis);
    CPPUNIT_ASSERT(Orientation::Sagittal == sagittal->orientation);
    CPPUNIT_ASSERT_EQUAL(1.0, sagittal->glyphSpacing);

    const auto coronal = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 23, 0, 0, 1, 0));
    CPPUNIT_ASSERT(coronal.has_value());
    CPPUNIT_ASSERT_EQUAL(2u, coronal->normalAxis);
    CPPUNIT_ASSERT(Orientation::Coronal == coronal->orientation);
  }

  // One glyph per voxel of the slice, centred on the plane; isotropic ODFs give spheres of radius
  // scaling * glyph spacing / 2.
  void GlyphGenerator_PositionsAndWindow()
  {
    auto image = MakeOdfImage(3, 2, 4, 1);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30 + 3 * 1.2, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    auto settings = AllGlyphs();
    settings.scaling = 2;
    auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(6 * ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
    CPPUNIT_ASSERT(glyphs->GetNumberOfPolys() > 0);
    CPPUNIT_ASSERT(glyphs->GetPointData()->GetArray("ODF_COLORS") != nullptr);
    CPPUNIT_ASSERT(glyphs->GetPointData()->GetNormals() != nullptr);
    for (vtkIdType g = 0; g < 6; ++g)
    {
      mitk::Point3D expected;
      expected[0] = 10 + (g % 3);
      expected[1] = 20 + 2 * (g / 3);
      expected[2] = 30 + 3 * 1.2;
      const auto center = GlyphCenter(glyphs, g);
      CPPUNIT_ASSERT(center.EuclideanDistanceTo(expected) < 1e-4);
      for (vtkIdType i = g * ODF_SAMPLING_SIZE; i < (g + 1) * ODF_SAMPLING_SIZE; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2 * 1.0 / 2, Distance(glyphs->GetPoint(i), expected), 1e-4);
    }

    // A window around the voxel (1, 1) keeps only its glyph.
    mitk::OdfGlyphGenerator::Window window;
    window.center[0] = 11;
    window.center[1] = 22;
    window.center[2] = 33.6;
    window.axis1[0] = 1;
    window.axis1[1] = 0;
    window.axis1[2] = 0;
    window.axis2[0] = 0;
    window.axis2[1] = 1;
    window.axis2[2] = 0;
    window.halfSize1 = 0.5;
    window.halfSize2 = 1;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, &window, settings);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
    CPPUNIT_ASSERT(GlyphCenter(glyphs, 0).EuclideanDistanceTo(window.center) < 1e-4);
  }

  // At most maxGlyphs glyphs; in striding mode they spread over the whole slice.
  void GlyphGenerator_MaxGlyphs_SpreadsOverSlice()
  {
    auto image = MakeOdfImage(10, 10, 1, 1);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30, 0, 0, 1));
    auto settings = AllGlyphs();
    settings.maxGlyphs = 30;
    auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    const vtkIdType numGlyphs = glyphs->GetNumberOfPoints() / ODF_SAMPLING_SIZE;
    CPPUNIT_ASSERT(numGlyphs > 0 && numGlyphs <= 30);
    CPPUNIT_ASSERT(GlyphCenter(glyphs, numGlyphs - 1)[1] >= 20 + 2 * 9 - 1e-4);

    settings.randomPlacement = true;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    CPPUNIT_ASSERT(glyphs->GetNumberOfPoints() / ODF_SAMPLING_SIZE <= 30);

    settings.maxGlyphs = 0;
    CPPUNIT_ASSERT_EQUAL(vtkIdType(0), mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings)->GetNumberOfPoints());
  }

  /** The voxels (x and y index) of the glyphs of an image of MakeOdfImage with constant ODFs. */
  static std::set<std::pair<int, int>> GlyphVoxels(vtkPolyData *glyphs)
  {
    std::set<std::pair<int, int>> voxels;
    for (vtkIdType g = 0; g < glyphs->GetNumberOfPoints() / ODF_SAMPLING_SIZE; ++g)
    {
      const auto center = GlyphCenter(glyphs, g);
      voxels.emplace(static_cast<int>(std::lround(center[0] - 10)), static_cast<int>(std::lround((center[1] - 20) / 2)));
    }
    return voxels;
  }

  /** Window over size x size voxels of an image of MakeOdfImage, from voxel (firstX, firstY). */
  static mitk::OdfGlyphGenerator::Window VoxelWindow(int firstX, int firstY, int size)
  {
    mitk::OdfGlyphGenerator::Window window;
    window.center[0] = 10 + firstX + (size - 1) / 2.0;
    window.center[1] = 20 + 2 * (firstY + (size - 1) / 2.0);
    window.center[2] = 30;
    window.axis1.Fill(0);
    window.axis1[0] = 1;
    window.axis2.Fill(0);
    window.axis2[1] = 1;
    window.halfSize1 = size / 2.0;
    window.halfSize2 = size;
    return window;
  }

  // Random placement keeps exactly maxGlyphs voxels, chosen by a fixed hash per voxel: the same
  // voxels at every call, a voxel kept with fewer glyphs is kept with more, a voxel kept in a view
  // is kept in a smaller view that contains it, and a pan keeps the glyphs that stay in view,
  // except near the selection threshold.
  void GlyphGenerator_RandomPlacement_StableSubset()
  {
    auto image = MakeOdfImage(40, 40, 1, 1);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    auto settings = AllGlyphs();
    settings.randomPlacement = true;
    const auto keptVoxels = [&](int maxGlyphs, const mitk::OdfGlyphGenerator::Window *window) {
      settings.maxGlyphs = maxGlyphs;
      return GlyphVoxels(mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, window, settings));
    };

    std::set<std::pair<int, int>> fewer;
    for (const int maxGlyphs : {1, 37, 400, 1599, 1600, 5000})
    {
      const auto kept = keptVoxels(maxGlyphs, nullptr);
      const std::string name = std::to_string(maxGlyphs) + " glyphs";
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, std::min<std::size_t>(maxGlyphs, 40 * 40), kept.size());
      CPPUNIT_ASSERT_MESSAGE(name, std::includes(kept.begin(), kept.end(), fewer.begin(), fewer.end()));
      CPPUNIT_ASSERT_MESSAGE(name, kept == keptVoxels(maxGlyphs, nullptr));
      fewer = kept;
    }

    // Zoomed in on 20 x 20 voxels.
    const auto window = VoxelWindow(10, 10, 20);
    const auto inWindow = keptVoxels(200, &window);
    CPPUNIT_ASSERT_EQUAL(std::size_t(200), inWindow.size());
    for (const auto &voxel : keptVoxels(200, nullptr))
    {
      if (voxel.first >= 10 && voxel.first < 30 && voxel.second >= 10 && voxel.second < 30)
        CPPUNIT_ASSERT(inWindow.count(voxel) == 1);
    }

    // Panned by 5 voxels. A new random subset of a quarter of the voxels would keep about a
    // quarter of the glyphs that stay in view.
    const auto panned = VoxelWindow(15, 10, 20);
    const auto before = keptVoxels(100, &window);
    const auto after = keptVoxels(100, &panned);
    CPPUNIT_ASSERT_EQUAL(std::size_t(100), after.size());
    int stayInView = 0;
    int kept = 0;
    for (const auto &voxel : before)
    {
      if (voxel.first >= 15)
      {
        ++stayInView;
        kept += static_cast<int>(after.count(voxel));
      }
    }
    CPPUNIT_ASSERT(stayInView > 50);
    CPPUNIT_ASSERT_MESSAGE(std::to_string(kept) + " of " + std::to_string(stayInView) + " glyphs kept", kept >= 0.9 * stayInView);
  }

  /** Radius of each point of the first glyph around its centre. */
  static std::vector<double> Radii(vtkPolyData *glyphs, const mitk::Point3D &center)
  {
    std::vector<double> radii;
    for (vtkIdType i = 0; i < ODF_SAMPLING_SIZE; ++i)
      radii.push_back(Distance(glyphs->GetPoint(i), center));
    return radii;
  }

  // The glyph arithmetic of the former vtkOdfSource: point = direction * value * scale *
  // glyph spacing / 2, colours from the normalized value (1 red, 0 blue) or the principal
  // direction.
  void GlyphGenerator_NormalizationScaleAndColours()
  {
    auto image = MakeOdfImage(1, 1, 1, 0);
    {
      mitk::ImagePixelWriteAccessor<itk::Vector<float, ODF_SAMPLING_SIZE>, 3> access(image);
      itk::Vector<float, ODF_SAMPLING_SIZE> odf;
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
        odf[i] = 1 + 3 * std::fabs(itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>::GetDirection(i)[0]);
      access.SetPixelByIndex({{0, 0, 0}}, odf);
    }
    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> odf;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = 1 + 3 * std::fabs(odf.GetDirection(i)[0]);
    const vtkIdType maxIndex = odf.GetPrincipalDiffusionDirectionIndex();
    vtkIdType minIndex = 0;
    for (unsigned int i = 1; i < ODF_SAMPLING_SIZE; ++i)
    {
      if (odf[i] < odf[minIndex])
        minIndex = i;
    }

    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30, 0, 0, 1));
    mitk::Point3D center;
    center[0] = 10;
    center[1] = 20;
    center[2] = 30;
    auto settings = AllGlyphs();
    settings.scaling = 4;

    settings.normalization = mitk::ODFN_MAX;
    auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    auto radii = Radii(glyphs, center);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4 * 1.0 / 2, radii[maxIndex], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(odf[minIndex] / odf[maxIndex] * 2, radii[minIndex], 1e-5);
    auto colors = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
    CPPUNIT_ASSERT_EQUAL(255, int(colors->GetTypedComponent(maxIndex, 0)));
    CPPUNIT_ASSERT_EQUAL(0, int(colors->GetTypedComponent(maxIndex, 2)));
    CPPUNIT_ASSERT_EQUAL(255, int(colors->GetTypedComponent(maxIndex, 3)));

    settings.normalization = mitk::ODFN_MINMAX;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    radii = Radii(glyphs, center);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, radii[maxIndex], 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, radii[minIndex], 1e-5);
    colors = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
    CPPUNIT_ASSERT_EQUAL(0, int(colors->GetTypedComponent(minIndex, 0)));
    CPPUNIT_ASSERT_EQUAL(255, int(colors->GetTypedComponent(minIndex, 2)));

    settings.normalization = mitk::ODFN_NONE;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    radii = Radii(glyphs, center);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(odf[maxIndex] * 2, radii[maxIndex], 1e-5);

    settings.normalization = mitk::ODFN_MAX;
    settings.scaleBy = mitk::ODFSB_GFA;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    radii = Radii(glyphs, center);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2 * odf.GetGeneralizedFractionalAnisotropy(), radii[maxIndex], 1e-5);

    settings.colourByDirection = true;
    glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    colors = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
    const auto direction = odf.GetDirection(maxIndex);
    for (vtkIdType i : {vtkIdType(0), vtkIdType(minIndex), vtkIdType(maxIndex)})
    {
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(std::fabs(direction[c]) * 255), int(colors->GetTypedComponent(i, c)));
    }
  }

  // Tensor glyphs in ODF mode are the ODFs that convert::GetOdfFromTensorImage computes, so the
  // generator reads the six tensor elements in the order of the tensor filters.
  void GlyphGenerator_TensorOdfs_MatchOdfImage()
  {
    auto tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti"));
    auto odfs = mitk::convert::GetOdfFromTensorImage(tensors);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*tensors, 0, *Plane(0, 0, tensors->GetGeometry()->GetCenter()[2], 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    auto tensorGlyphs = mitk::OdfGlyphGenerator::Generate(*tensors, 0, *slice, nullptr, AllGlyphs());
    auto odfGlyphs = mitk::OdfGlyphGenerator::Generate(*odfs, 0, *slice, nullptr, AllGlyphs());
    CPPUNIT_ASSERT(tensorGlyphs->GetNumberOfPoints() > 0);
    CPPUNIT_ASSERT_EQUAL(odfGlyphs->GetNumberOfPoints(), tensorGlyphs->GetNumberOfPoints());
    for (vtkIdType i = 0; i < tensorGlyphs->GetNumberOfPoints(); ++i)
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(odfGlyphs->GetPoint(i)[c], tensorGlyphs->GetPoint(i)[c], 1e-5);
  }

  // The ellipsoid of diag(3, 1, 1) has the semi-axes 3, 1, 1 along x, y, z: the glyph point in
  // direction u lies at 1 / sqrt(ux^2/9 + uy^2 + uz^2), normalized to the largest of the sampled
  // directions, times half the spacing from the voxel centre.
  void GlyphGenerator_TensorEllipsoid_Semiaxes()
  {
    using ItkTensorImageType = mitk::TensorImage::ItkTensorImageType;
    auto itkImage = ItkTensorImageType::New();
    ItkTensorImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    ItkTensorImageType::PixelType tensor;
    tensor.Fill(0);
    tensor(0, 0) = 3;
    tensor(1, 1) = 1;
    tensor(2, 2) = 1;
    itkImage->FillBuffer(tensor);
    auto image = mitk::TensorImage::New();
    image->InitializeByItk(itkImage.GetPointer());
    image->SetVolume(itkImage->GetBufferPointer());

    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    auto settings = AllGlyphs();
    settings.tensorEllipsoids = true;
    auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
    std::vector<double> radii;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      const auto u = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>::GetDirection(i);
      radii.push_back(1 / std::sqrt(u[0] * u[0] / 9 + u[1] * u[1] + u[2] * u[2]));
    }
    const double largest = *std::max_element(radii.begin(), radii.end());
    for (vtkIdType i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      const auto u = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>::GetDirection(i);
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(u[c] * radii[i] / largest * 0.5, glyphs->GetPoint(i)[c], 1e-5);
    }
  }

  static mitk::ShImage::Pointer LoadShImage()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"SH Image"}, std::vector<std::string>());
    return mitk::IOUtil::Load<mitk::ShImage>(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), &functor);
  }

  // The glyphs of an SH image follow its convention: they equal those of its ODF image.
  void GlyphGenerator_ShMatchesOdfOfBothConventions()
  {
    auto shImage = LoadShImage();
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*shImage, 0, *Plane(0, 0, shImage->GetGeometry()->GetOrigin()[2] + 2 * shImage->GetGeometry()->GetSpacing()[2], 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    vtkSmartPointer<vtkPolyData> previous;
    for (auto convention : {mitk::ShImage::MRTRIX, mitk::ShImage::FSL})
    {
      shImage->SetShConvention(convention);
      auto shGlyphs = mitk::OdfGlyphGenerator::Generate(*shImage, 0, *slice, nullptr, AllGlyphs());
      auto odfGlyphs = mitk::OdfGlyphGenerator::Generate(*mitk::convert::GetOdfFromShImage(shImage.GetPointer()), 0, *slice, nullptr, AllGlyphs());
      CPPUNIT_ASSERT(shGlyphs->GetNumberOfPoints() > 0);
      CPPUNIT_ASSERT_EQUAL(odfGlyphs->GetNumberOfPoints(), shGlyphs->GetNumberOfPoints());
      for (vtkIdType i = 0; i < shGlyphs->GetNumberOfPoints(); ++i)
      {
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(odfGlyphs->GetPoint(i)[c], shGlyphs->GetPoint(i)[c], 1e-4);
      }

      if (previous)
      {
        double maxDifference = 0;
        for (vtkIdType i = 0; i < shGlyphs->GetNumberOfPoints(); ++i)
          for (int c = 0; c < 3; ++c)
            maxDifference = std::max(maxDifference, std::fabs(previous->GetPoint(i)[c] - shGlyphs->GetPoint(i)[c]));
        CPPUNIT_ASSERT(maxDifference > 0.01);
      }
      previous = shGlyphs;
    }
  }

  // Six SH coefficients are order 2, not tensor elements: c0 alone is an isotropic ODF.
  void GlyphGenerator_Order2ShIsNoTensor()
  {
    using ShVectorImageType = itk::Image<itk::Vector<float, 6>, 3>;
    auto itkImage = ShVectorImageType::New();
    ShVectorImageType::RegionType region;
    region.SetSize({{1, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    itk::Vector<float, 6> coefficients;
    coefficients.Fill(0);
    coefficients[0] = 1;
    itkImage->FillBuffer(coefficients);
    auto shImage = mitk::ShImage::New();
    shImage->InitializeByItk(itkImage.GetPointer());
    shImage->SetVolume(itkImage->GetBufferPointer());

    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*shImage, 0, *Plane(0, 0, 0, 0, 0, 1));
    auto settings = AllGlyphs();
    settings.scaleBy = mitk::ODFSB_NONE;
    auto glyphs = mitk::OdfGlyphGenerator::Generate(*shImage, 0, *slice, nullptr, settings);
    mitk::Point3D center;
    center.Fill(0);
    for (const double radius : Radii(glyphs, center))
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, radius, 1e-5);
  }

  void GlyphGenerator_UnsupportedImage_Throws()
  {
    auto image = mitk::Image::New();
    image->Initialize(mitk::MakeScalarPixelType<float>(), 3, std::array<unsigned int, 3>{{2, 2, 2}}.data());
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    CPPUNIT_ASSERT_THROW(mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, AllGlyphs()), mitk::Exception);
  }

  using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;

  /** FNV-1a hash of raw bytes. */
  class Fingerprint
  {
  public:
    void Add(const void *data, std::size_t size)
    {
      const auto *bytes = static_cast<const unsigned char *>(data);
      for (std::size_t i = 0; i < size; ++i)
      {
        m_Hash ^= bytes[i];
        m_Hash *= 1099511628211ull;
      }
    }

    std::uint64_t Get() const { return m_Hash; }

  private:
    std::uint64_t m_Hash = 14695981039346656037ull;
  };

  /** Hashes of the float bits of the points, of the point ids of the triangles and of the RGBA colours. */
  struct GlyphFingerprints
  {
    std::uint64_t points = 0;
    std::uint64_t triangles = 0;
    std::uint64_t colours = 0;

    bool operator==(const GlyphFingerprints &) const = default;
  };

  static GlyphFingerprints Fingerprints(vtkPolyData *glyphs)
  {
    Fingerprint points;
    for (vtkIdType i = 0; i < glyphs->GetNumberOfPoints(); ++i)
    {
      double p[3];
      glyphs->GetPoint(i, p);
      const float coordinates[3] = {static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])};
      points.Add(coordinates, sizeof(coordinates));
    }

    Fingerprint triangles;
    auto ids = vtkSmartPointer<vtkIdList>::New();
    vtkCellArray *polys = glyphs->GetPolys();
    for (vtkIdType c = 0; c < polys->GetNumberOfCells(); ++c)
    {
      polys->GetCellAtId(c, ids);
      for (vtkIdType k = 0; k < ids->GetNumberOfIds(); ++k)
      {
        const std::int64_t id = ids->GetId(k);
        triangles.Add(&id, sizeof(id));
      }
    }

    Fingerprint colours;
    auto rgba = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
    CPPUNIT_ASSERT(rgba != nullptr);
    CPPUNIT_ASSERT_EQUAL(4, rgba->GetNumberOfComponents());
    CPPUNIT_ASSERT_EQUAL(glyphs->GetNumberOfPoints(), rgba->GetNumberOfTuples());
    for (vtkIdType i = 0; i < rgba->GetNumberOfTuples(); ++i)
    {
      for (int c = 0; c < 4; ++c)
      {
        const unsigned char value = rgba->GetTypedComponent(i, c);
        colours.Add(&value, 1);
      }
    }
    return {points.Get(), triangles.Get(), colours.Get()};
  }

  static std::uint64_t NormalsFingerprint(vtkPolyData *glyphs)
  {
    vtkDataArray *normals = glyphs->GetPointData()->GetNormals();
    CPPUNIT_ASSERT(normals != nullptr);
    Fingerprint fingerprint;
    for (vtkIdType i = 0; i < normals->GetNumberOfTuples(); ++i)
    {
      for (int c = 0; c < 3; ++c)
      {
        const auto value = static_cast<float>(normals->GetComponent(i, c));
        fingerprint.Add(&value, sizeof(value));
      }
    }
    return fingerprint.Get();
  }

  /** Every glyph is the ODF base mesh: its triangles are the base triangles on its own points. */
  static void CheckTopology(vtkPolyData *glyphs, const std::string &name)
  {
    vtkCellArray *baseTriangles = OdfType::GetBaseMesh()->GetPolys();
    const vtkIdType numBaseTriangles = baseTriangles->GetNumberOfCells();
    const vtkIdType numGlyphs = glyphs->GetNumberOfPoints() / ODF_SAMPLING_SIZE;
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name, numGlyphs * ODF_SAMPLING_SIZE, glyphs->GetNumberOfPoints());
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name, numGlyphs * numBaseTriangles, glyphs->GetNumberOfPolys());
    auto baseIds = vtkSmartPointer<vtkIdList>::New();
    auto ids = vtkSmartPointer<vtkIdList>::New();
    vtkIdType mismatches = 0;
    for (vtkIdType g = 0; g < numGlyphs; ++g)
    {
      for (vtkIdType t = 0; t < numBaseTriangles; ++t)
      {
        baseTriangles->GetCellAtId(t, baseIds);
        glyphs->GetPolys()->GetCellAtId(g * numBaseTriangles + t, ids);
        if (ids->GetNumberOfIds() != 3 || baseIds->GetNumberOfIds() != 3)
        {
          ++mismatches;
          continue;
        }
        for (vtkIdType k = 0; k < 3; ++k)
        {
          if (ids->GetId(k) != baseIds->GetId(k) + g * ODF_SAMPLING_SIZE)
            ++mismatches;
        }
      }
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name, vtkIdType(0), mismatches);
  }

  /**
   * Largest difference between the point normals of the glyphs and those vtkPolyDataNormals
   * computes on their points and triangles, set up as the generator once used it. Infinite where
   * only one of the two is NaN.
   */
  static double NormalsDeviation(vtkPolyData *glyphs)
  {
    vtkDataArray *normals = glyphs->GetPointData()->GetNormals();
    CPPUNIT_ASSERT(normals != nullptr);
    CPPUNIT_ASSERT_EQUAL(std::string("Normals"), std::string(normals->GetName()));
    CPPUNIT_ASSERT_EQUAL(3, normals->GetNumberOfComponents());

    // Points and triangles only: vtkPolyDataNormals keeps normals its input already has.
    auto mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->SetPoints(glyphs->GetPoints());
    mesh->SetPolys(glyphs->GetPolys());
    auto filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    filter->SetInputData(mesh);
    filter->SplittingOff();
    filter->ConsistencyOff();
    filter->AutoOrientNormalsOff();
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOff();
    filter->FlipNormalsOff();
    filter->NonManifoldTraversalOff();
    filter->Update();
    vtkDataArray *expected = filter->GetOutput()->GetPointData()->GetNormals();
    CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfTuples(), normals->GetNumberOfTuples());

    double deviation = 0;
    for (vtkIdType i = 0; i < normals->GetNumberOfTuples(); ++i)
    {
      for (int c = 0; c < 3; ++c)
      {
        const double actual = normals->GetComponent(i, c);
        const double reference = expected->GetComponent(i, c);
        if (std::isnan(actual) || std::isnan(reference))
        {
          if (!(std::isnan(actual) && std::isnan(reference)))
            deviation = std::numeric_limits<double>::infinity();
        }
        else
        {
          deviation = std::max(deviation, std::fabs(actual - reference));
        }
      }
    }
    return deviation;
  }

  struct GlyphCase
  {
    std::string name;
    mitk::Image::Pointer image;
    mitk::OdfGlyphGenerator::Slice slice;
    std::optional<mitk::OdfGlyphGenerator::Window> window;
    mitk::OdfGlyphGenerator::Settings settings;

    vtkSmartPointer<vtkPolyData> Generate() const
    {
      return mitk::OdfGlyphGenerator::Generate(*image, 0, slice, window ? &*window : nullptr, settings);
    }
  };

  /** ODF image of three voxels along x: a zero ODF, a constant ODF and one with negative values. */
  static mitk::OdfImage::Pointer MakeDegenerateOdfImage()
  {
    auto image = MakeOdfImage(3, 1, 1, 0);
    {
      mitk::ImagePixelWriteAccessor<itk::Vector<float, ODF_SAMPLING_SIZE>, 3> access(image);
      itk::Vector<float, ODF_SAMPLING_SIZE> odf;
      odf.Fill(2);
      access.SetPixelByIndex({{1, 0, 0}}, odf);
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
        odf[i] = static_cast<float>(1.5 * OdfType::GetDirection(i)[0] - 0.5);
      access.SetPixelByIndex({{2, 0, 0}}, odf);
    }
    return image;
  }

  /** Tensor image of three voxels along x: a null tensor, one with a negative eigenvalue and a rotated diag(3, 1, 1). */
  static mitk::TensorImage::Pointer MakeDegenerateTensorImage()
  {
    using ItkTensorImageType = mitk::TensorImage::ItkTensorImageType;
    auto itkImage = ItkTensorImageType::New();
    ItkTensorImageType::RegionType region;
    region.SetSize({{3, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    ItkTensorImageType::PixelType tensor;
    tensor.Fill(0);
    itkImage->FillBuffer(tensor);
    tensor(0, 0) = 1;
    tensor(1, 1) = 1;
    tensor(2, 2) = -1;
    itkImage->SetPixel({{1, 0, 0}}, tensor);
    // diag(3, 1, 1) rotated by 30 degrees about z.
    const double c = std::cos(itk::Math::pi / 6);
    const double s = std::sin(itk::Math::pi / 6);
    tensor.Fill(0);
    tensor(0, 0) = 3 * c * c + s * s;
    tensor(0, 1) = 2 * c * s;
    tensor(1, 1) = 3 * s * s + c * c;
    tensor(2, 2) = 1;
    itkImage->SetPixel({{2, 0, 0}}, tensor);
    auto image = mitk::TensorImage::New();
    image->InitializeByItk(itkImage.GetPointer());
    image->SetVolume(itkImage->GetBufferPointer());
    return image;
  }

  /**
   * One slice of each kind of glyph: tensors as ellipsoids and as tensor ODFs coloured by
   * direction, ODFs with value colours at each normalization and scale mode, SH in both
   * conventions, and glyphs that collapse to their centre (zero, constant, null, non-positive).
   * tensorImage.dti and qBallImage.qbi (64x64x3, spacing 3) through a window of 16x16 voxels of
   * layer 1, which holds both isotropic and directed voxels.
   */
  static std::vector<GlyphCase> GlyphCases()
  {
    std::vector<GlyphCase> cases;
    const auto add = [&cases](const std::string &name, mitk::Image *image, const mitk::PlaneGeometry &plane, const std::optional<mitk::OdfGlyphGenerator::Window> &window, const std::function<void(mitk::OdfGlyphGenerator::Settings &)> &set) {
      const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, plane);
      CPPUNIT_ASSERT_MESSAGE(name, slice.has_value());
      auto settings = AllGlyphs();
      set(settings);
      cases.push_back({name, image, *slice, window, settings});
    };

    mitk::OdfGlyphGenerator::Window window;
    window.center[0] = 94.5;
    window.center[1] = 94.5;
    window.center[2] = 3;
    window.axis1[0] = 1;
    window.axis1[1] = 0;
    window.axis1[2] = 0;
    window.axis2[0] = 0;
    window.axis2[1] = 1;
    window.axis2[2] = 0;
    window.halfSize1 = 24;
    window.halfSize2 = 24;
    const auto layer1 = Plane(0, 0, 3, 0, 0, 1);

    mitk::Image::Pointer tensors = mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti")).GetPointer();
    for (const bool ellipsoids : {true, false})
    {
      add(ellipsoids ? "tensor ellipsoids" : "tensor ODFs", tensors, *layer1, window, [ellipsoids](mitk::OdfGlyphGenerator::Settings &settings) {
        settings.tensorEllipsoids = ellipsoids;
        settings.colourByDirection = true;
        settings.scaling = 3;
        settings.scaleBy = mitk::ODFSB_GFA;
      });
    }

    mitk::Image::Pointer odfs = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi")).GetPointer();
    const std::pair<int, std::string> normalizations[] = {
      {mitk::ODFN_MINMAX, "min-max"}, {mitk::ODFN_MAX, "max"}, {mitk::ODFN_NONE, "none"}, {mitk::ODFN_GLOBAL_MAX, "global max"}};
    const std::pair<int, std::string> scaleModes[] = {{mitk::ODFSB_NONE, "unscaled"}, {mitk::ODFSB_GFA, "GFA"}, {mitk::ODFSB_PC, "PC"}};
    for (const auto &[normalization, normalizationName] : normalizations)
    {
      for (const auto &[scaleBy, scaleName] : scaleModes)
      {
        add("ODFs " + normalizationName + " " + scaleName, odfs, *layer1, window, [normalization, scaleBy](mitk::OdfGlyphGenerator::Settings &settings) {
          settings.scaling = 1.5;
          settings.normalization = normalization;
          settings.scaleBy = scaleBy;
          // The mapper's default angles (2, 1) give an empty range and a curvature of 0 everywhere.
          settings.indexParam1 = 0;
          settings.indexParam2 = 30;
        });
      }
    }

    for (const auto convention : {mitk::ShImage::MRTRIX, mitk::ShImage::FSL})
    {
      auto sh = LoadShImage();
      sh->SetShConvention(convention);
      mitk::Point3D layer1Index;
      layer1Index.Fill(0);
      layer1Index[2] = 1;
      mitk::Point3D origin;
      sh->GetGeometry()->IndexToWorld(layer1Index, origin);
      add(convention == mitk::ShImage::MRTRIX ? "SH MRtrix" : "SH FSL", sh.GetPointer(), *Plane(origin[0], origin[1], origin[2], 0, 0, 1), std::nullopt, [](mitk::OdfGlyphGenerator::Settings &settings) {
        settings.scaling = 1.5;
        settings.scaleBy = mitk::ODFSB_GFA;
      });
    }

    mitk::Image::Pointer degenerateOdfs = MakeDegenerateOdfImage().GetPointer();
    for (const auto &[normalization, normalizationName] : normalizations)
    {
      add("degenerate ODFs " + normalizationName, degenerateOdfs, *Plane(0, 0, 30, 0, 0, 1), std::nullopt, [normalization](mitk::OdfGlyphGenerator::Settings &settings) {
        settings.normalization = normalization;
      });
    }
    mitk::Image::Pointer degenerateTensors = MakeDegenerateTensorImage().GetPointer();
    for (const bool ellipsoids : {true, false})
    {
      add(ellipsoids ? "degenerate tensor ellipsoids" : "degenerate tensor ODFs", degenerateTensors, *Plane(0, 0, 0, 0, 0, 1), std::nullopt, [ellipsoids](mitk::OdfGlyphGenerator::Settings &settings) {
        settings.tensorEllipsoids = ellipsoids;
        settings.colourByDirection = true;
      });
    }
    return cases;
  }

  // The full-mesh output as the generator computed it with per-point VTK calls and
  // vtkPolyDataNormals (develop ea962a7): points bit for bit, the same triangles and colours,
  // normals within 1e-5 of vtkPolyDataNormals on the same points and triangles.
  void GlyphGenerator_FullMesh_RecordedOutput()
  {
    const std::map<std::string, GlyphFingerprints> recorded = {
      {"tensor ellipsoids", {17034852577485430665ull, 1096128973618458029ull, 490934932846469741ull}},
      {"tensor ODFs", {12851842336964751317ull, 1096128973618458029ull, 2635710187966396005ull}},
      {"ODFs min-max unscaled", {14189319946489299714ull, 1096128973618458029ull, 7742319593747850464ull}},
      {"ODFs min-max GFA", {17668013067538098052ull, 1096128973618458029ull, 7742319593747850464ull}},
      {"ODFs min-max PC", {16332114341006182612ull, 1096128973618458029ull, 7742319593747850464ull}},
      {"ODFs max unscaled", {9345829136554767028ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs max GFA", {4520074477892272862ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs max PC", {8738574457543285195ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs none unscaled", {9333783507626493847ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs none GFA", {6025879742403316006ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs none PC", {12122150632538098043ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs global max unscaled", {9345829136554767028ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs global max GFA", {4520074477892272862ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"ODFs global max PC", {8738574457543285195ull, 1096128973618458029ull, 9733053046497571580ull}},
      {"SH MRtrix", {12865566614912398697ull, 4467218837662528689ull, 1804916545293251655ull}},
      {"SH FSL", {10487421149490399071ull, 4467218837662528689ull, 4990161006980516369ull}},
      {"degenerate ODFs min-max", {5204507226439308963ull, 5386312073842232526ull, 17532495340253511400ull}},
      {"degenerate ODFs max", {17697105185732528215ull, 5386312073842232526ull, 16719975996596679099ull}},
      {"degenerate ODFs none", {8926904338823488585ull, 5386312073842232526ull, 16719975996596679099ull}},
      {"degenerate ODFs global max", {17697105185732528215ull, 5386312073842232526ull, 16719975996596679099ull}},
      {"degenerate tensor ellipsoids", {7152420504689445845ull, 5386312073842232526ull, 3723008392808557773ull}},
      {"degenerate tensor ODFs", {13356892293473713405ull, 5386312073842232526ull, 2899782707417845261ull}},
    };
    std::ostringstream mismatches;
    for (const auto &glyphCase : GlyphCases())
    {
      auto glyphs = glyphCase.Generate();
      CPPUNIT_ASSERT_MESSAGE(glyphCase.name, glyphs->GetNumberOfPoints() > 0);
      CheckTopology(glyphs, glyphCase.name);
      const double deviation = NormalsDeviation(glyphs);
      CPPUNIT_ASSERT_MESSAGE(glyphCase.name + " normals deviate by " + std::to_string(deviation), deviation <= 1e-5);
      const auto fingerprints = Fingerprints(glyphs);
      const auto expected = recorded.find(glyphCase.name);
      if (expected == recorded.end() || !(expected->second == fingerprints))
      {
        mismatches << "\n{\"" << glyphCase.name << "\", {" << fingerprints.points << "ull, " << fingerprints.triangles
                   << "ull, " << fingerprints.colours << "ull}},";
      }
    }
    CPPUNIT_ASSERT_MESSAGE("Glyphs that differ from the recorded output:" + mismatches.str(), mismatches.str().empty());
  }

  // Value colours are the colours of vtkLookupTable::GetColor for 1 - the normalized value,
  // converted with 255 * colour; values outside [0, 1] get the colour of the nearer end.
  void GlyphGenerator_ValueColours_FromLookupTable()
  {
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetRange(0, 1);
    lut->Build();
    const auto expectedColour = [&lut](float normalizedValue) {
      double rgb[3];
      lut->GetColor(1 - normalizedValue, rgb);
      std::array<int, 4> rgba = {0, 0, 0, 255};
      for (unsigned int c = 0; c < 3; ++c)
        rgba[c] = static_cast<unsigned char>(255.0 * rgb[c]);
      return rgba;
    };
    const auto checkColours = [&expectedColour](vtkPolyData *glyphs, vtkIdType glyph, const OdfType &colourOdf, const std::string &name) {
      auto colours = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
      CPPUNIT_ASSERT(colours != nullptr);
      for (vtkIdType j = 0; j < ODF_SAMPLING_SIZE; ++j)
      {
        const auto rgba = expectedColour(colourOdf[j]);
        for (int c = 0; c < 4; ++c)
          CPPUNIT_ASSERT_EQUAL_MESSAGE(name + " point " + std::to_string(j), rgba[c], int(colours->GetTypedComponent(glyph * ODF_SAMPLING_SIZE + j, c)));
      }
    };
    // Min-max colours the normalized values; the other normalizations colour the values divided by the maximum.
    const auto colourOdf = [](const OdfType &odf, int normalization) {
      return normalization == mitk::ODFN_MINMAX ? odf.MinMaxNormalize() : odf.MaxNormalize();
    };

    // Values from -2 to 1: below 0 after division by the maximum.
    auto image = MakeDegenerateOdfImage();
    OdfType odf;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = static_cast<float>(1.5 * OdfType::GetDirection(i)[0] - 0.5);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30, 0, 0, 1));
    for (const int normalization : {mitk::ODFN_MINMAX, mitk::ODFN_MAX, mitk::ODFN_NONE, mitk::ODFN_GLOBAL_MAX})
    {
      auto settings = AllGlyphs();
      settings.normalization = normalization;
      auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
      CPPUNIT_ASSERT_EQUAL(vtkIdType(3 * ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
      checkColours(glyphs, 2, colourOdf(odf, normalization), "synthetic, normalization " + std::to_string(normalization));
    }

    // The 16x16 voxels of qBallImage.qbi around the centre of layer 1, row by row.
    auto odfs = mitk::IOUtil::Load<mitk::OdfImage>(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi"));
    const auto cases = GlyphCases();
    const auto odfCase = std::find_if(cases.begin(), cases.end(), [](const GlyphCase &c) { return c.name == "ODFs min-max unscaled"; });
    CPPUNIT_ASSERT(odfCase != cases.end());
    mitk::ImagePixelReadAccessor<itk::Vector<float, ODF_SAMPLING_SIZE>, 3> access(odfs.GetPointer());
    for (const int normalization : {mitk::ODFN_MINMAX, mitk::ODFN_MAX})
    {
      auto glyphCase = *odfCase;
      glyphCase.settings.normalization = normalization;
      auto glyphs = glyphCase.Generate();
      CPPUNIT_ASSERT_EQUAL(vtkIdType(256 * ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
      for (vtkIdType g = 0; g < 256; ++g)
      {
        const auto values = access.GetPixelByIndex({{static_cast<itk::IndexValueType>(24 + g % 16), static_cast<itk::IndexValueType>(24 + g / 16), 1}});
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = values[i];
        checkColours(glyphs, g, colourOdf(odf, normalization), "qBallImage glyph " + std::to_string(g));
      }
    }
  }

  // Each glyph is computed on its own, so the output does not depend on the number of threads.
  void GlyphGenerator_ThreadCount_SameOutput()
  {
    auto cases = GlyphCases();
    // The whole layer, more glyphs than threads by far.
    cases.front().window.reset();
    // The cases at the coarse meshes too.
    for (std::size_t i = 0, numCases = cases.size(); i < numCases; ++i)
    {
      for (const auto &[width, n] : coarseLevels)
      {
        auto coarse = cases[i];
        coarse.name += ", " + std::to_string(n) + " points";
        coarse.settings.pixelsPerMm = PixelsPerMm(coarse.slice, coarse.settings.scaling, width);
        cases.push_back(coarse);
      }
    }
    mitk::OmpThreadGuard threads(1);
    for (const auto &glyphCase : cases)
    {
      threads.SetNumberOfThreads(1);
      auto reference = glyphCase.Generate();
      const auto referenceFingerprints = Fingerprints(reference);
      const auto referenceNormals = NormalsFingerprint(reference);
      for (const int numThreads : {3, 8})
      {
        threads.SetNumberOfThreads(numThreads);
        auto glyphs = glyphCase.Generate();
        const std::string name = glyphCase.name + ", " + std::to_string(numThreads) + " threads";
        CPPUNIT_ASSERT_MESSAGE(name, referenceFingerprints == Fingerprints(glyphs));
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name, referenceNormals, NormalsFingerprint(glyphs));
      }
    }
  }

  /** Pixels per mm that make the glyphs of the slice widthInPixels wide at the given scaling. */
  static double PixelsPerMm(const mitk::OdfGlyphGenerator::Slice &slice, double scaling, double widthInPixels)
  {
    return widthInPixels / (scaling * slice.glyphSpacing);
  }

  /** Direction j of the mesh with n points: the pointshell table of n directions. */
  static vnl_vector_fixed<double, 3> MeshDirection(unsigned int n, unsigned int j)
  {
    const auto &directions = mitk::pointshell::Directions(n);
    return vnl_vector_fixed<double, 3>(directions(0, j), directions(1, j), directions(2, j));
  }

  // Glyphs narrower than 8 and 20 pixels on screen (scaling x glyph spacing x pixels per mm)
  // are geodesic spheres of 42 and 162 directions; pixels per mm 0 keeps the full mesh.
  void GlyphGenerator_MeshLevel_ByGlyphWidth()
  {
    std::vector<std::pair<std::string, mitk::Image::Pointer>> images = {
      {"tensor", MakeDegenerateTensorImage().GetPointer()}, {"SH", LoadShImage().GetPointer()}, {"ODF", MakeOdfImage(2, 1, 1, 1).GetPointer()}};
    const std::tuple<double, double, vtkIdType> cases[] = {
      {1, 0, ODF_SAMPLING_SIZE}, {1, 0.01, 42}, {1, 7.9, 42}, {1, 8.1, 162}, {1, 19.9, 162}, {1, 20.1, ODF_SAMPLING_SIZE},
      {2.5, 7.9, 42}, {2.5, 8.1, 162}, {2.5, 19.9, 162}, {2.5, 20.1, ODF_SAMPLING_SIZE}};
    for (const auto &[name, image] : images)
    {
      mitk::Point3D layer1;
      image->GetGeometry()->IndexToWorld(mitk::Point3D(0.0), layer1);
      const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(layer1[0], layer1[1], layer1[2], 0, 0, 1));
      CPPUNIT_ASSERT_MESSAGE(name, slice.has_value());
      auto settings = AllGlyphs();
      settings.tensorEllipsoids = true;
      const vtkIdType numGlyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings)->GetNumberOfPoints() / ODF_SAMPLING_SIZE;
      CPPUNIT_ASSERT_MESSAGE(name, numGlyphs > 0);
      for (const auto &[scaling, width, pointsPerGlyph] : cases)
      {
        settings.scaling = static_cast<float>(scaling);
        settings.pixelsPerMm = PixelsPerMm(*slice, scaling, width);
        auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
        const std::string message = name + ", scaling " + std::to_string(scaling) + ", " + std::to_string(width) + " pixels";
        CPPUNIT_ASSERT_EQUAL_MESSAGE(message, numGlyphs * pointsPerGlyph, glyphs->GetNumberOfPoints());
        CPPUNIT_ASSERT_EQUAL_MESSAGE(message, numGlyphs * (2 * pointsPerGlyph - 4), glyphs->GetNumberOfPolys());
      }
    }

  }

  /** Glyph widths in pixels that give the coarse meshes, and their numbers of points. */
  static constexpr std::pair<double, vtkIdType> coarseLevels[] = {{5, 42}, {15, 162}};

  // A coarse glyph is a closed surface on the 42 or 162 directions of mitk::pointshell, in their
  // order, with triangles counterclockwise seen from outside and the normals of vtkPolyDataNormals.
  void GlyphGenerator_CoarseMesh_ClosedSurface()
  {
    auto image = MakeDegenerateTensorImage();
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));
    for (const auto &[width, n] : coarseLevels)
    for (const bool ellipsoids : {true, false})
    {
      auto settings = AllGlyphs();
      settings.tensorEllipsoids = ellipsoids;
      settings.pixelsPerMm = PixelsPerMm(*slice, 1, width);
      auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
      const vtkIdType numTriangles = 2 * n - 4;
      CPPUNIT_ASSERT_EQUAL(3 * n, glyphs->GetNumberOfPoints());
      CPPUNIT_ASSERT_EQUAL(3 * numTriangles, glyphs->GetNumberOfPolys());
      CPPUNIT_ASSERT(NormalsDeviation(glyphs) <= 1e-5);

      // Glyph 2, the rotated diag(3, 1, 1), has positive values in every direction.
      const double center[3] = {2, 0, 0};
      for (vtkIdType j = 0; j < n; ++j)
      {
        mitk::Vector3D direction;
        for (unsigned int c = 0; c < 3; ++c)
          direction[c] = glyphs->GetPoint(2 * n + j)[c] - center[c];
        direction.Normalize();
        CPPUNIT_ASSERT((direction - mitk::Vector3D(MeshDirection(n, j))).GetNorm() < 1e-5);
      }

      std::map<std::pair<vtkIdType, vtkIdType>, int> directedEdges;
      auto ids = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType t = 2 * numTriangles; t < 3 * numTriangles; ++t)
      {
        glyphs->GetPolys()->GetCellAtId(t, ids);
        CPPUNIT_ASSERT_EQUAL(vtkIdType(3), ids->GetNumberOfIds());
        double p[3][3];
        for (vtkIdType k = 0; k < 3; ++k)
        {
          CPPUNIT_ASSERT(ids->GetId(k) >= 2 * n && ids->GetId(k) < 3 * n);
          glyphs->GetPoint(ids->GetId(k), p[k]);
          ++directedEdges[{ids->GetId(k), ids->GetId((k + 1) % 3)}];
        }
        double normal[3];
        const double ab[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
        const double ac[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
        vtkMath::Cross(ab, ac, normal);
        double outward = 0;
        for (unsigned int c = 0; c < 3; ++c)
          outward += normal[c] * ((p[0][c] + p[1][c] + p[2][c]) / 3 - center[c]);
        CPPUNIT_ASSERT(outward > 0);
      }
      // Closed and consistently oriented: every edge once in each direction.
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(3 * numTriangles), directedEdges.size());
      for (const auto &[edge, count] : directedEdges)
      {
        CPPUNIT_ASSERT_EQUAL(1, count);
        CPPUNIT_ASSERT(directedEdges.count({edge.second, edge.first}) == 1);
      }
    }
  }

  // At a coarse level the glyph point j lies in mesh direction u_j at the tensor's value there,
  // the ellipsoid's radius 1 / sqrt(u^T D^-2 u) or the tensor ODF u^T D u, limited to the range
  // of the values at the ODF directions and normalized by the largest of them. The colour by
  // direction is that of the full mesh.
  void GlyphGenerator_CoarseMesh_TensorValuesAtItsDirections()
  {
    auto image = MakeDegenerateTensorImage();
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));
    // Voxel 2: diag(3, 1, 1) rotated by 30 degrees about z.
    const double c = std::cos(itk::Math::pi / 6);
    const double s = std::sin(itk::Math::pi / 6);
    vnl_matrix_fixed<double, 3, 3> tensor(0.0);
    tensor(0, 0) = 3 * c * c + s * s;
    tensor(0, 1) = tensor(1, 0) = 2 * c * s;
    tensor(1, 1) = 3 * s * s + c * c;
    tensor(2, 2) = 1;
    const vnl_matrix_fixed<double, 3, 3> inverse = vnl_inverse(tensor);
    const vnl_matrix_fixed<double, 3, 3> inverseSquared = inverse * inverse;
    const auto value = [&](bool ellipsoids, const vnl_vector_fixed<double, 3> &u) {
      return ellipsoids ? 1 / std::sqrt(dot_product(u, inverseSquared * u)) : dot_product(u, tensor * u);
    };

    for (const auto &[width, n] : coarseLevels)
    for (const bool ellipsoids : {true, false})
    {
      double largest = 0;
      double smallest = std::numeric_limits<double>::max();
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      {
        largest = std::max(largest, value(ellipsoids, OdfType::GetDirection(i)));
        smallest = std::min(smallest, value(ellipsoids, OdfType::GetDirection(i)));
      }

      auto settings = AllGlyphs();
      settings.tensorEllipsoids = ellipsoids;
      settings.colourByDirection = true;
      auto full = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
      settings.pixelsPerMm = PixelsPerMm(*slice, 1, width);
      auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
      CPPUNIT_ASSERT_EQUAL(3 * n, glyphs->GetNumberOfPoints());
      for (vtkIdType j = 0; j < n; ++j)
      {
        const mitk::Vector3D u = MeshDirection(n, j);
        const double radius = std::clamp(value(ellipsoids, MeshDirection(n, j)), smallest, largest) / largest * 0.5;
        mitk::Point3D expected;
        expected.Fill(0);
        expected[0] = 2;
        expected += u * radius;
        CPPUNIT_ASSERT(Distance(glyphs->GetPoint(2 * n + j), expected) < 1e-5);
        // The null tensor and the one with a negative eigenvalue have no ellipsoid.
        if (ellipsoids)
        {
          CPPUNIT_ASSERT(Distance(glyphs->GetPoint(j), mitk::Point3D(0.0)) == 0);
          CPPUNIT_ASSERT(Distance(glyphs->GetPoint(n + j), [] { mitk::Point3D p(0.0); p[0] = 1; return p; }()) == 0);
        }
      }
      auto colours = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
      auto fullColours = vtkUnsignedCharArray::SafeDownCast(full->GetPointData()->GetArray("ODF_COLORS"));
      for (int k = 0; k < 4; ++k)
        CPPUNIT_ASSERT_EQUAL(int(fullColours->GetTypedComponent(2 * ODF_SAMPLING_SIZE, k)), int(colours->GetTypedComponent(2 * n, k)));
    }
  }

  // At a coarse level the glyph point j of an SH image lies in mesh direction u_j at the SH
  // series there, limited to the range of the values at the ODF directions, normalized by the
  // largest of them and scaled by the GFA of the ODF there; its colour is that of the normalized
  // value.
  void GlyphGenerator_CoarseMesh_ShValuesAtItsDirections()
  {
    auto shImage = LoadShImage();
    const unsigned int numCoefficients = shImage->GetPixelType().GetNumberOfComponents();
    const int order = static_cast<int>(mitk::sh::ShOrder(static_cast<int>(numCoefficients)));
    mitk::Point3D layer1Index;
    layer1Index.Fill(0);
    layer1Index[2] = 1;
    mitk::Point3D origin;
    shImage->GetGeometry()->IndexToWorld(layer1Index, origin);
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*shImage, 0, *Plane(origin[0], origin[1], origin[2], 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value() && slice->normalAxis == 2 && slice->index == 1);

    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetRange(0, 1);
    lut->Build();
    const auto colour = [&lut](double normalizedValue) {
      double rgb[3];
      lut->GetColor(1 - normalizedValue, rgb);
      std::array<int, 3> result;
      for (unsigned int k = 0; k < 3; ++k)
        result[k] = static_cast<unsigned char>(255.0 * rgb[k]);
      return result;
    };

    mitk::ImageReadAccessor access(shImage.GetPointer(), shImage->GetVolumeData(0));
    const auto *values = static_cast<const float *>(access.GetData());
    for (const auto &[width, n] : coarseLevels)
    for (const auto convention : {mitk::ShImage::MRTRIX, mitk::ShImage::FSL})
    {
      shImage->SetShConvention(convention);
      const bool mrtrix = convention == mitk::ShImage::MRTRIX;
      auto settings = AllGlyphs();
      settings.scaling = 1.5;
      settings.scaleBy = mitk::ODFSB_GFA;
      settings.pixelsPerMm = PixelsPerMm(*slice, settings.scaling, width);
      auto glyphs = mitk::OdfGlyphGenerator::Generate(*shImage, 0, *slice, nullptr, settings);
      CPPUNIT_ASSERT_EQUAL(100 * n, glyphs->GetNumberOfPoints());
      auto colours = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
      int checkedGlyphs = 0;
      for (vtkIdType g = 0; g < 100; ++g)
      {
        const vnl_vector<float> coefficients(values + ((g % 10) + 10 * ((g / 10) + 10 * 1)) * numCoefficients, numCoefficients);
        OdfType odf;
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = mitk::sh::GetValue(coefficients, order, OdfType::GetDirection(i), mrtrix);
        const double largest = odf.GetMaxValue();
        if (!(largest > 0))
          continue;
        const double smallest = *std::min_element(odf.Begin(), odf.End());
        ++checkedGlyphs;
        const double scale = 1.5 * odf.GetGeneralizedFractionalAnisotropy() * slice->glyphSpacing * 0.5;
        mitk::Point3D center;
        mitk::Point3D index;
        index[0] = static_cast<double>(g % 10);
        index[1] = static_cast<double>(g / 10);
        index[2] = 1;
        shImage->GetGeometry()->IndexToWorld(index, center);
        for (vtkIdType j = 0; j < n; ++j)
        {
          const mitk::Vector3D u = MeshDirection(n, j);
          const double meshValue = mitk::sh::GetValue(coefficients, order, MeshDirection(n, j), mrtrix);
          const double normalized = std::clamp(meshValue, smallest, largest) / largest;
          const std::string name = "glyph " + std::to_string(g) + " point " + std::to_string(j);
          CPPUNIT_ASSERT_MESSAGE(name, Distance(glyphs->GetPoint(g * n + j), center + u * (normalized * scale)) < 1e-4);
          // The colour of the normalized value, or of a value next to it when it lies at the
          // border of two colours (float sums in another order).
          std::array<int, 3> actual;
          for (unsigned int k = 0; k < 3; ++k)
            actual[k] = colours->GetTypedComponent(g * n + j, k);
          CPPUNIT_ASSERT_MESSAGE(name, actual == colour(normalized) || actual == colour(normalized - 1e-5) || actual == colour(normalized + 1e-5));
        }
      }
      CPPUNIT_ASSERT(checkedGlyphs > 50);
    }
  }

  // ODF images have values only at the ODF directions: at a coarse level, glyph point j lies in
  // mesh direction u_j at the value of the ODF direction nearest to u_j, normalized by the ODF's
  // range, and has the colour of that normalized value.
  void GlyphGenerator_CoarseMesh_OdfValuesOfNearestDirection()
  {
    auto image = MakeOdfImage(1, 1, 1, 0);
    OdfType odf;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = 1 + static_cast<float>(i) / ODF_SAMPLING_SIZE;
    {
      mitk::ImagePixelWriteAccessor<itk::Vector<float, ODF_SAMPLING_SIZE>, 3> access(image);
      itk::Vector<float, ODF_SAMPLING_SIZE> values;
      std::copy(odf.Begin(), odf.End(), values.Begin());
      access.SetPixelByIndex({{0, 0, 0}}, values);
    }
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 30, 0, 0, 1));
    mitk::Point3D center;
    center[0] = 10;
    center[1] = 20;
    center[2] = 30;
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetRange(0, 1);
    lut->Build();

    for (const auto &[width, n] : coarseLevels)
    {
      for (const int normalization : {mitk::ODFN_MINMAX, mitk::ODFN_MAX})
      {
        auto settings = AllGlyphs();
        settings.normalization = normalization;
        settings.pixelsPerMm = PixelsPerMm(*slice, 1, width);
        auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
        CPPUNIT_ASSERT_EQUAL(n, glyphs->GetNumberOfPoints());
        const OdfType normalized = normalization == mitk::ODFN_MINMAX ? odf.MinMaxNormalize() : odf.MaxNormalize();
        auto colours = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
        for (vtkIdType j = 0; j < n; ++j)
        {
          const auto u = MeshDirection(n, j);
          unsigned int nearest = 0;
          for (unsigned int i = 1; i < ODF_SAMPLING_SIZE; ++i)
          {
            if (dot_product(u, OdfType::GetDirection(i)) > dot_product(u, OdfType::GetDirection(nearest)))
              nearest = i;
          }
          const std::string name = std::to_string(n) + " points, normalization " + std::to_string(normalization) + ", point " + std::to_string(j);
          CPPUNIT_ASSERT_MESSAGE(name, Distance(glyphs->GetPoint(j), center + mitk::Vector3D(u) * (normalized[nearest] * 0.5)) < 1e-5);
          double rgb[3];
          lut->GetColor(1 - normalized[nearest], rgb);
          for (unsigned int c = 0; c < 3; ++c)
            CPPUNIT_ASSERT_EQUAL_MESSAGE(name, int(static_cast<unsigned char>(255.0 * rgb[c])), int(colours->GetTypedComponent(j, c)));
        }
      }
    }
  }

  // The ellipsoid of a needle-like tensor is far larger along its axis than at the nearest ODF
  // direction; a coarse direction closer to the axis must not make the glyph larger than the full
  // one, which reaches half the glyph spacing times the scaling.
  void GlyphGenerator_CoarseMesh_NoLargerThanFullGlyph()
  {
    using ItkTensorImageType = mitk::TensorImage::ItkTensorImageType;
    auto itkImage = ItkTensorImageType::New();
    ItkTensorImageType::RegionType region;
    region.SetSize({{8, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    // Needles of eigenvalues 1, 1e-3, 1e-3 along directions between those of the meshes.
    for (itk::IndexValueType x = 0; x < 8; ++x)
    {
      vnl_vector_fixed<double, 3> axis(std::cos(0.37 * x + 0.1), std::sin(0.37 * x + 0.1), 0.3 * x - 1.1);
      axis.normalize();
      ItkTensorImageType::PixelType tensor;
      tensor.Fill(0);
      for (unsigned int r = 0; r < 3; ++r)
        for (unsigned int c = r; c < 3; ++c)
          tensor(r, c) = (r == c ? 1e-3 : 0) + (1 - 1e-3) * axis[r] * axis[c];
      itkImage->SetPixel({{x, 0, 0}}, tensor);
    }
    auto image = mitk::TensorImage::New();
    image->InitializeByItk(itkImage.GetPointer());
    image->SetVolume(itkImage->GetBufferPointer());
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));

    for (const bool ellipsoids : {true, false})
    {
      for (const double width : {0.0, 5.0, 15.0})
      {
        auto settings = AllGlyphs();
        settings.tensorEllipsoids = ellipsoids;
        settings.pixelsPerMm = PixelsPerMm(*slice, 1, width);
        auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
        const vtkIdType n = glyphs->GetNumberOfPoints() / 8;
        for (vtkIdType g = 0; g < 8; ++g)
        {
          mitk::Point3D center(0.0);
          center[0] = static_cast<double>(g);
          double largest = 0;
          for (vtkIdType j = 0; j < n; ++j)
            largest = std::max(largest, Distance(glyphs->GetPoint(g * n + j), center));
          const std::string name = (ellipsoids ? "ellipsoid " : "tensor ODF ") + std::to_string(g) + ", " + std::to_string(n) + " points";
          CPPUNIT_ASSERT_MESSAGE(name + ": " + std::to_string(largest), largest <= 0.5 + 1e-6);
          if (n == ODF_SAMPLING_SIZE)
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(name, 0.5, largest, 1e-6);
        }
      }
    }
  }

  /** R diag(semi-axes) R^T of instance g: its ellipsoid as a symmetric matrix. */
  static vnl_matrix_fixed<double, 3, 3> InstanceEllipsoid(vtkPolyData *instances, vtkIdType g)
  {
    double quaternion[4];
    instances->GetPointData()->GetArray("Orientation")->GetTuple(g, quaternion);
    double rotation[3][3];
    vtkMath::QuaternionToMatrix3x3(quaternion, rotation);
    double axes[3];
    instances->GetPointData()->GetArray("Scale")->GetTuple(g, axes);
    vnl_matrix_fixed<double, 3, 3> ellipsoid;
    for (int r = 0; r < 3; ++r)
    {
      for (int c = 0; c < 3; ++c)
        ellipsoid(r, c) = rotation[r][0] * axes[0] * rotation[c][0] + rotation[r][1] * axes[1] * rotation[c][1] + rotation[r][2] * axes[2] * rotation[c][2];
    }
    return ellipsoid;
  }

  /** Hash of the raw bytes of the instances' points and point arrays. */
  static std::uint64_t InstancesFingerprint(vtkPolyData *instances)
  {
    Fingerprint fingerprint;
    for (vtkDataArray *array : {instances->GetPoints()->GetData(), instances->GetPointData()->GetArray("Orientation"),
                                instances->GetPointData()->GetArray("Scale"), instances->GetPointData()->GetArray("ODF_COLORS")})
    {
      CPPUNIT_ASSERT(array != nullptr);
      fingerprint.Add(array->GetVoidPointer(0), array->GetNumberOfValues() * array->GetDataTypeSize());
    }
    return fingerprint.Get();
  }

  // The instanced ellipsoids are drawn for the voxels of the glyph meshes, in their colours, each
  // the exact ellipsoid of its tensor D with the largest radius of the glyph mesh:
  // R diag(semi-axes) R^T = D x largest radius / largest eigenvalue, semi-axes in ascending order.
  // This holds for any choice of eigenvectors, also of equal eigenvalues.
  void GlyphGenerator_EllipsoidInstances_MatchTheTensors()
  {
    // Three different eigenvalues fix the eigenvectors up to their signs, about half of which make
    // a reflection; at scaling 2 and spacing 1 the largest semi-axis is 1.
    using ItkTensorImageType = mitk::TensorImage::ItkTensorImageType;
    auto itkImage = ItkTensorImageType::New();
    ItkTensorImageType::RegionType region;
    region.SetSize({{16, 1, 1}});
    itkImage->SetRegions(region);
    itkImage->Allocate();
    for (itk::IndexValueType x = 0; x < 16; ++x)
    {
      auto transform = vtkSmartPointer<vtkTransform>::New();
      transform->RotateZ(40.0 * x + 10);
      transform->RotateX(25.0 * x + 5);
      transform->RotateY(15.0 * x);
      ItkTensorImageType::PixelType tensor;
      for (unsigned int r = 0; r < 3; ++r)
      {
        for (unsigned int c = r; c < 3; ++c)
        {
          double value = 0;
          for (unsigned int k = 0; k < 3; ++k)
            value += transform->GetMatrix()->GetElement(r, k) * (3.0 - k) * transform->GetMatrix()->GetElement(c, k);
          tensor(r, c) = value;
        }
      }
      itkImage->SetPixel({{x, 0, 0}}, tensor);
    }
    auto rotatedTensors = mitk::TensorImage::New();
    rotatedTensors->InitializeByItk(itkImage.GetPointer());
    rotatedTensors->SetVolume(itkImage->GetBufferPointer());
    const auto rotatedSlice = mitk::OdfGlyphGenerator::GetSlice(*rotatedTensors, 0, *Plane(0, 0, 0, 0, 0, 1));
    CPPUNIT_ASSERT(rotatedSlice.has_value());
    auto settings = AllGlyphs();
    settings.tensorEllipsoids = true;
    settings.colourByDirection = true;
    settings.scaling = 2;
    auto rotatedInstances = mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(*rotatedTensors, 0, *rotatedSlice, nullptr, settings);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(16), rotatedInstances->GetNumberOfPoints());
    int reflections = 0;
    for (itk::IndexValueType x = 0; x < 16; ++x)
    {
      const auto tensor = itkImage->GetPixel({{x, 0, 0}});
      ItkTensorImageType::PixelType::EigenValuesArrayType eigenValues;
      ItkTensorImageType::PixelType::EigenVectorsMatrixType eigenVectors;
      tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
      vnl_matrix_fixed<double, 3, 3> eigenBase;
      for (int r = 0; r < 3; ++r)
      {
        for (int c = 0; c < 3; ++c)
          eigenBase(r, c) = eigenVectors(c, r);
      }
      reflections += vnl_det(eigenBase) < 0 ? 1 : 0;
      const auto ellipsoid = InstanceEllipsoid(rotatedInstances, x);
      for (int r = 0; r < 3; ++r)
      {
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("rotated diag(3, 2, 1) " + std::to_string(x), tensor(r, c) / 3, ellipsoid(r, c), 1e-6);
      }
    }
    CPPUNIT_ASSERT(reflections > 0 && reflections < 16);

    auto cases = GlyphCases();
    std::vector<GlyphCase> instanceCases;
    for (const auto &glyphCase : cases)
    {
      if (glyphCase.name == "tensor ellipsoids" || glyphCase.name == "degenerate tensor ellipsoids")
        instanceCases.push_back(glyphCase);
    }
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), instanceCases.size());
    // The stable random subset of the whole layer.
    auto random = instanceCases.front();
    random.name = "tensor ellipsoids, random 300 of the layer";
    random.window.reset();
    random.settings.maxGlyphs = 300;
    random.settings.randomPlacement = true;
    instanceCases.push_back(random);

    for (const auto &glyphCase : instanceCases)
    {
      const mitk::OdfGlyphGenerator::Window *window = glyphCase.window ? &*glyphCase.window : nullptr;
      CPPUNIT_ASSERT_MESSAGE(glyphCase.name, mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(*glyphCase.image, glyphCase.settings));
      auto instances = mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(*glyphCase.image, 0, glyphCase.slice, window, glyphCase.settings);
      auto glyphs = glyphCase.Generate();
      const vtkIdType numGlyphs = glyphs->GetNumberOfPoints() / ODF_SAMPLING_SIZE;
      CPPUNIT_ASSERT_MESSAGE(glyphCase.name, numGlyphs > 0);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(glyphCase.name, numGlyphs, instances->GetNumberOfPoints());
      auto instanceColours = vtkUnsignedCharArray::SafeDownCast(instances->GetPointData()->GetArray("ODF_COLORS"));
      auto glyphColours = vtkUnsignedCharArray::SafeDownCast(glyphs->GetPointData()->GetArray("ODF_COLORS"));
      CPPUNIT_ASSERT(instanceColours != nullptr && glyphColours != nullptr);

      const mitk::ImageReadAccessor access(glyphCase.image, glyphCase.image->GetVolumeData(0));
      const auto *values = static_cast<const float *>(access.GetData());
      for (vtkIdType g = 0; g < numGlyphs; ++g)
      {
        const std::string name = glyphCase.name + ", glyph " + std::to_string(g);
        double center[3];
        instances->GetPoint(g, center);
        CPPUNIT_ASSERT_MESSAGE(name, Distance(center, GlyphCenter(glyphs, g)) < 1e-3);
        for (int c = 0; c < 4; ++c)
          CPPUNIT_ASSERT_EQUAL_MESSAGE(name, glyphColours->GetTypedComponent(g * ODF_SAMPLING_SIZE, c), instanceColours->GetTypedComponent(g, c));

        mitk::Point3D centerPoint;
        centerPoint[0] = center[0];
        centerPoint[1] = center[1];
        centerPoint[2] = center[2];
        double largestRadius = 0;
        for (vtkIdType j = g * ODF_SAMPLING_SIZE; j < (g + 1) * ODF_SAMPLING_SIZE; ++j)
          largestRadius = std::max(largestRadius, Distance(glyphs->GetPoint(j), centerPoint));

        itk::Index<3> index;
        glyphCase.image->GetGeometry()->WorldToIndex(centerPoint, index);
        const itk::IndexValueType offset = index[0] + glyphCase.image->GetDimension(0) * (index[1] + glyphCase.image->GetDimension(1) * index[2]);
        const itk::DiffusionTensor3D<float> tensor(values + 6 * offset);
        itk::DiffusionTensor3D<float>::EigenValuesArrayType eigenValues;
        itk::DiffusionTensor3D<float>::EigenVectorsMatrixType eigenVectors;
        tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
        const bool hasEllipsoid = eigenValues[0] > 0 && eigenValues[1] > 0 && eigenValues[2] > 0;

        double axes[3];
        instances->GetPointData()->GetArray("Scale")->GetTuple(g, axes);
        CPPUNIT_ASSERT_MESSAGE(name, axes[0] <= axes[1] && axes[1] <= axes[2]);
        // The glyph mesh's points and the centre are rounded to float.
        const double tolerance = 1e-4 + 1e-5 * largestRadius;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(name, hasEllipsoid ? largestRadius : 0.0, axes[2], tolerance);
        const auto ellipsoid = InstanceEllipsoid(instances, g);
        for (int r = 0; r < 3; ++r)
        {
          for (int c = 0; c < 3; ++c)
          {
            const double expected = hasEllipsoid ? tensor(r, c) * largestRadius / eigenValues[2] : 0;
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(name, expected, ellipsoid(r, c), tolerance);
          }
        }
      }
    }

    // Each instance is computed on its own, so the output does not depend on the number of threads.
    auto wholeLayer = instanceCases.front();
    wholeLayer.window.reset();
    mitk::OmpThreadGuard threads(1);
    const auto reference = InstancesFingerprint(mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(*wholeLayer.image, 0, wholeLayer.slice, nullptr, wholeLayer.settings));
    for (const int numThreads : {3, 8})
    {
      threads.SetNumberOfThreads(numThreads);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(std::to_string(numThreads) + " threads", reference,
                                   InstancesFingerprint(mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(*wholeLayer.image, 0, wholeLayer.slice, nullptr, wholeLayer.settings)));
    }
  }

  // Colours by value change along a glyph and the other normalizations change its shape: only
  // tensor ellipsoids coloured by direction and normalized by their maximum are instances.
  void GlyphGenerator_EllipsoidInstances_OnlyForTheirSettings()
  {
    auto tensors = MakeDegenerateTensorImage();
    auto odfs = MakeOdfImage(2, 1, 1, 1);
    auto settings = AllGlyphs();
    settings.tensorEllipsoids = true;
    settings.colourByDirection = true;
    settings.normalization = mitk::ODFN_MAX;
    CPPUNIT_ASSERT(mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(*tensors, settings));
    CPPUNIT_ASSERT(!mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(*odfs, settings));
    CPPUNIT_ASSERT(!mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(*LoadShImage(), settings));
    const std::pair<std::string, std::function<void(mitk::OdfGlyphGenerator::Settings &)>> changes[] = {
      {"global maximum", [](mitk::OdfGlyphGenerator::Settings &s) { s.normalization = mitk::ODFN_GLOBAL_MAX; }},
      {"tensor ODFs", [](mitk::OdfGlyphGenerator::Settings &s) { s.tensorEllipsoids = false; }},
      {"colours by value", [](mitk::OdfGlyphGenerator::Settings &s) { s.colourByDirection = false; }},
      {"min-max", [](mitk::OdfGlyphGenerator::Settings &s) { s.normalization = mitk::ODFN_MINMAX; }},
      {"no normalization", [](mitk::OdfGlyphGenerator::Settings &s) { s.normalization = mitk::ODFN_NONE; }}};
    for (const auto &[name, change] : changes)
    {
      auto changed = settings;
      change(changed);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, name == "global maximum", mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(*tensors, changed));
    }

    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*odfs, 0, *Plane(10, 20, 30, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    CPPUNIT_ASSERT_THROW(mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(*odfs, 0, *slice, nullptr, settings), mitk::Exception);
  }

  // The instances' sphere is the mesh of the glyph meshes at the same zoom, on the unit sphere,
  // with the sphere's normals.
  void GlyphGenerator_InstanceSphere_MeshOfTheLevel()
  {
    auto image = MakeDegenerateTensorImage();
    const auto slice = mitk::OdfGlyphGenerator::GetSlice(*image, 0, *Plane(0, 0, 0, 0, 0, 1));
    CPPUNIT_ASSERT(slice.has_value());
    const std::pair<double, vtkIdType> levels[] = {{0, ODF_SAMPLING_SIZE}, {5, 42}, {15, 162}, {100, ODF_SAMPLING_SIZE}};
    for (const auto &[width, n] : levels)
    {
      const std::string name = std::to_string(width) + " pixels";
      auto settings = AllGlyphs();
      settings.tensorEllipsoids = true;
      settings.pixelsPerMm = PixelsPerMm(*slice, settings.scaling, width);
      auto sphere = mitk::OdfGlyphGenerator::GetInstanceSphere(*slice, settings);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, n, sphere->GetNumberOfPoints());
      auto glyphs = mitk::OdfGlyphGenerator::Generate(*image, 0, *slice, nullptr, settings);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 3 * n, glyphs->GetNumberOfPoints());

      // The triangles of the first glyph, whose points start at 0.
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 2 * n - 4, sphere->GetNumberOfPolys());
      auto sphereIds = vtkSmartPointer<vtkIdList>::New();
      auto glyphIds = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType c = 0; c < sphere->GetNumberOfPolys(); ++c)
      {
        sphere->GetPolys()->GetCellAtId(c, sphereIds);
        glyphs->GetPolys()->GetCellAtId(c, glyphIds);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name, vtkIdType(3), sphereIds->GetNumberOfIds());
        for (vtkIdType k = 0; k < 3; ++k)
          CPPUNIT_ASSERT_EQUAL_MESSAGE(name, glyphIds->GetId(k), sphereIds->GetId(k));
      }

      vtkDataArray *normals = sphere->GetPointData()->GetNormals();
      CPPUNIT_ASSERT(normals != nullptr);
      for (vtkIdType j = 0; j < n; ++j)
      {
        double expected[3];
        if (n == ODF_SAMPLING_SIZE)
        {
          OdfType::GetBaseMesh()->GetPoint(j, expected);
        }
        else
        {
          const auto direction = MeshDirection(static_cast<unsigned int>(n), static_cast<unsigned int>(j));
          std::copy(direction.begin(), direction.end(), expected);
        }
        vtkMath::Normalize(expected);
        double point[3];
        sphere->GetPoint(j, point);
        double normal[3];
        normals->GetTuple(j, normal);
        for (int c = 0; c < 3; ++c)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(name, expected[c], point[c], 1e-6);
          CPPUNIT_ASSERT_EQUAL_MESSAGE(name, point[c], normal[c]);
        }
      }
    }
  }

  // Defaults no longer overwrite what the node already has, unless asked to.
  void OdfMapper_DefaultProperties_KeepExistingValues()
  {
    auto node = mitk::DataNode::New();
    node->SetData(mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("ImageReconstruction/test_dti.dti")));
    node->SetBoolProperty("VisibleOdfs_T", true);
    node->SetFloatProperty("Scaling", 7);
    node->SetIntProperty("layer", 3);

    mitk::OdfVtkMapper2D::SetDefaultProperties(node);
    bool visibleOdfs = false;
    float scaling = 0;
    int layer = 0;
    CPPUNIT_ASSERT(node->GetBoolProperty("VisibleOdfs_T", visibleOdfs) && visibleOdfs);
    CPPUNIT_ASSERT(node->GetFloatProperty("Scaling", scaling) && scaling == 7);
    CPPUNIT_ASSERT(node->GetIntProperty("layer", layer) && layer == 3);
    CPPUNIT_ASSERT(node->GetProperty("DoRefresh") == nullptr);

    mitk::OdfVtkMapper2D::SetDefaultProperties(node, nullptr, true);
    CPPUNIT_ASSERT(node->GetBoolProperty("VisibleOdfs_T", visibleOdfs) && !visibleOdfs);
    CPPUNIT_ASSERT(node->GetFloatProperty("Scaling", scaling) && scaling == 3);
    CPPUNIT_ASSERT(node->GetIntProperty("layer", layer) && layer == 100);
  }

  /** A view of layer 1 of an image, with a window of 20 x 10 mm around (5, 6, 3). */
  static mitk::OdfVtkMapper2D::GlyphView MakeGlyphView()
  {
    mitk::OdfVtkMapper2D::GlyphView view;
    view.slice.normalAxis = 2;
    view.slice.index = 1;
    view.slice.planeIndex = 1;
    view.planeAxis0.Fill(0);
    view.planeAxis0[0] = 100;
    view.planeAxis1.Fill(0);
    view.planeAxis1[1] = 80;
    view.window.center[0] = 5;
    view.window.center[1] = 6;
    view.window.center[2] = 3;
    view.window.axis1.Fill(0);
    view.window.axis1[0] = 1;
    view.window.axis2.Fill(0);
    view.window.axis2[1] = 1;
    view.window.halfSize1 = 10;
    view.window.halfSize2 = 5;
    view.pixelsPerMm = 2;
    return view;
  }

  // The rule of the 2D mapper: glyphs are generated the first time, for new content and for
  // another slice or plane at once; for a pan or zoom only at the render of LOD 1 that follows
  // an interaction. A render at LOD 1 without a change keeps them.
  void OdfMapper_GenerationRequired()
  {
    using GlyphView = mitk::OdfVtkMapper2D::GlyphView;
    const GlyphView generated = MakeGlyphView();
    for (const int lod : {0, 1})
    {
      const std::string level = ", LOD " + std::to_string(lod);
      CPPUNIT_ASSERT_MESSAGE("first" + level, mitk::OdfVtkMapper2D::IsGenerationRequired(std::nullopt, generated, false, lod));
      CPPUNIT_ASSERT_MESSAGE("unchanged" + level, !mitk::OdfVtkMapper2D::IsGenerationRequired(generated, generated, false, lod));
      CPPUNIT_ASSERT_MESSAGE("content" + level, mitk::OdfVtkMapper2D::IsGenerationRequired(generated, generated, true, lod));
    }

    const std::pair<std::string, std::function<void(GlyphView &)>> sliceChanges[] = {
      {"normal axis", [](GlyphView &v) { v.slice.normalAxis = 1; }},
      {"slice index", [](GlyphView &v) { v.slice.index = 2; }},
      {"plane index", [](GlyphView &v) { v.slice.planeIndex = 1.2; }},
      {"plane axis 0", [](GlyphView &v) { v.planeAxis0[1] = 1; }},
      {"plane axis 1", [](GlyphView &v) { v.planeAxis1[0] = 1; }}};
    const std::pair<std::string, std::function<void(GlyphView &)>> viewChanges[] = {
      {"pan", [](GlyphView &v) { v.window.center[0] += 0.5; }},
      {"zoom", [](GlyphView &v) { v.window.halfSize1 *= 1.1; v.window.halfSize2 *= 1.1; }},
      {"window axes", [](GlyphView &v) { std::swap(v.window.axis1, v.window.axis2); }},
      {"pixels per mm", [](GlyphView &v) { v.pixelsPerMm *= 2; }}};
    for (const int lod : {0, 1})
    {
      for (const auto &[name, change] : sliceChanges)
      {
        GlyphView current = generated;
        change(current);
        CPPUNIT_ASSERT_MESSAGE(name + ", LOD " + std::to_string(lod), mitk::OdfVtkMapper2D::IsGenerationRequired(generated, current, false, lod));
      }
      for (const auto &[name, change] : viewChanges)
      {
        GlyphView current = generated;
        change(current);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ", LOD " + std::to_string(lod), lod > 0, mitk::OdfVtkMapper2D::IsGenerationRequired(generated, current, false, lod));
      }
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingRendering)
