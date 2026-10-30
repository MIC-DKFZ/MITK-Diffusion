/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkIOUtil.h>
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
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
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
  MITK_TEST(GlyphGenerator_NormalizationScaleAndColours);
  MITK_TEST(GlyphGenerator_ShMatchesOdfOfBothConventions);
  MITK_TEST(GlyphGenerator_Order2ShIsNoTensor);
  MITK_TEST(GlyphGenerator_TensorOdfs_MatchOdfImage);
  MITK_TEST(GlyphGenerator_TensorEllipsoid_Semiaxes);
  MITK_TEST(GlyphGenerator_UnsupportedImage_Throws);
  MITK_TEST(OdfMapper_DefaultProperties_KeepExistingValues);
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingRendering)
