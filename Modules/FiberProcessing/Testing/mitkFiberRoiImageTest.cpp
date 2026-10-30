/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkPlanarFigureRoiImage.h>

#include <mitkImageCast.h>
#include <mitkPlanarCircle.h>
#include <mitkPlanarPolygon.h>
#include <mitkSliceNavigationController.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkEuler3DTransform.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <cmath>
#include <map>
#include <string>

/** The ROI image of planar figures (Fiber Processing view, "Generate ROI Image"), on synthetic
 * images and figures whose covered voxels are known. */
class mitkFiberRoiImageTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberRoiImageTestSuite);
  MITK_TEST(SlicePlanes_SetExactlyTheVoxelsOfTheSlice);
  MITK_TEST(Circle_NoHolesAndAreaOfTheCircle);
  MITK_TEST(ObliquePlane_NoHolesAndAreaOfTheFigure);
  MITK_TEST(PlaneOnSliceBorder_SetsOneLayer);
  MITK_TEST(UnplacedFigure_SetsNothing);
  MITK_TEST(FigurePartlyOutside_SetsTheVoxelsInside);
  MITK_TEST(TwoFigures_UnionAndOtherVoxelsKept);
  CPPUNIT_TEST_SUITE_END();

  using ImageType = itk::Image<unsigned char, 3>;

  /** Empty image rotated by rx degrees about x and rz degrees about z. */
  static ImageType::Pointer CreateImage(const itk::Vector<double, 3>& spacing,
                                        const itk::Point<double, 3>& origin,
                                        double rx = 0,
                                        double rz = 0)
  {
    auto image = ImageType::New();
    ImageType::SizeType size = {{40, 40, 30}};
    image->SetRegions(ImageType::RegionType(size));
    image->SetSpacing(spacing);
    image->SetOrigin(origin);
    auto rotation = itk::Euler3DTransform<double>::New();
    rotation->SetRotation(rx * itk::Math::pi / 180, 0, rz * itk::Math::pi / 180);
    ImageType::DirectionType direction;
    direction = rotation->GetMatrix();
    image->SetDirection(direction);
    image->Allocate();
    image->FillBuffer(0);
    return image;
  }

  static itk::Vector<double, 3> Spacing(double x, double y, double z)
  {
    itk::Vector<double, 3> spacing;
    spacing[0] = x;
    spacing[1] = y;
    spacing[2] = z;
    return spacing;
  }

  static itk::Point<double, 3> Point(double x, double y, double z)
  {
    itk::Point<double, 3> point;
    point[0] = x;
    point[1] = y;
    point[2] = z;
    return point;
  }

  /** The plane a render window shows for the given view, through the voxel centre c, as the
   * Workbench gives it to a figure drawn there. */
  static mitk::PlaneGeometry::Pointer SlicePlane(ImageType* image, mitk::AnatomicalPlane view, mitk::Point3D& c)
  {
    mitk::Image::Pointer mitkImage;
    mitk::CastToMitkImage(ImageType::Pointer(image), mitkImage);
    auto controller = mitk::SliceNavigationController::New();
    controller->SetInputWorldTimeGeometry(mitkImage->GetTimeGeometry());
    controller->Update(view);
    itk::Index<3> centre;
    for (unsigned int i = 0; i < 3; ++i)
      centre[i] = image->GetLargestPossibleRegion().GetSize(i) / 2;
    image->TransformIndexToPhysicalPoint(centre, c);
    controller->SelectSliceByPoint(c);
    return controller->GetCurrentPlaneGeometry()->Clone();
  }

  /** Oblique plane through c with the orthonormal in-plane axes u and v. */
  static mitk::PlaneGeometry::Pointer ObliquePlane(const mitk::Point3D& c, const mitk::Vector3D& u, const mitk::Vector3D& v)
  {
    auto plane = mitk::PlaneGeometry::New();
    plane->InitializeStandardPlane(u * 100.0, v * 100.0);
    plane->SetOrigin(c - u * 50.0 - v * 50.0);
    return plane;
  }

  /** Rectangle with half sizes a and b along the plane's axes u and v around the point c. */
  static mitk::PlanarPolygon::Pointer Rectangle(mitk::PlaneGeometry* plane, const mitk::Point3D& c, double a, double b)
  {
    mitk::Vector3D u = plane->GetAxisVector(0);
    mitk::Vector3D v = plane->GetAxisVector(1);
    u.Normalize();
    v.Normalize();
    auto figure = mitk::PlanarPolygon::New();
    figure->SetPlaneGeometry(plane);
    const double corners[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    for (int k = 0; k < 4; ++k)
    {
      mitk::Point2D point;
      plane->Map(c + u * (corners[k][0] * a) + v * (corners[k][1] * b), point);
      if (k == 0)
        figure->PlaceFigure(point);
      else
        figure->AddControlPoint(point);
    }
    figure->SetClosed(true);
    return figure;
  }

  static unsigned int CountSet(const ImageType* image)
  {
    unsigned int count = 0;
    for (itk::ImageRegionConstIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
      count += it.Get() > 0 ? 1 : 0;
    return count;
  }

  /** Thickness of a voxel of the image along the unit vector n. */
  static double Thickness(const ImageType* image, const mitk::Vector3D& n)
  {
    double thickness = 0;
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      double cosine = 0;
      for (unsigned int row = 0; row < 3; ++row)
        cosine += n[row] * image->GetDirection()[row][axis];
      thickness += image->GetSpacing()[axis] * std::abs(cosine);
    }
    return thickness;
  }

public:
  /** The cases of issue #3: a rectangle drawn in an axial, sagittal or coronal render window on
   * images with anisotropic or rotated voxels. Expected are the voxels of the figure's slice whose
   * centres lie inside the rectangle; the former 2D-grid rasterization left holes in sagittal and
   * coronal slices (72 of 135 voxels for 2 x 2 x 1 mm) and shifted the area by up to one voxel. */
  void SlicePlanes_SetExactlyTheVoxelsOfTheSlice()
  {
    const auto origin = Point(-37.3, 12.1, -5.7);
    const std::map<std::string, ImageType::Pointer> images = {
      {"iso2_origin0", CreateImage(Spacing(2, 2, 2), Point(0, 0, 0))},
      {"iso2_origin_frac", CreateImage(Spacing(2, 2, 2), origin)},
      {"1.5x1.5x3", CreateImage(Spacing(1.5, 1.5, 3), origin)},
      {"2x2x1", CreateImage(Spacing(2, 2, 1), origin)},
      {"iso2_rotated", CreateImage(Spacing(2, 2, 2), origin, 20, 10)}};
    const std::map<std::string, mitk::AnatomicalPlane> views = {
      {"axial", mitk::AnatomicalPlane::Axial},
      {"sagittal", mitk::AnatomicalPlane::Sagittal},
      {"coronal", mitk::AnatomicalPlane::Coronal}};
    const double a = 9.3;
    const double b = 7.1;

    for (const auto& [imageName, image] : images)
    {
      for (const auto& [viewName, view] : views)
      {
        const std::string name = imageName + "_" + viewName;
        image->FillBuffer(0);
        mitk::Point3D c;
        auto plane = SlicePlane(image, view, c);
        mitk::Vector3D u = plane->GetAxisVector(0);
        mitk::Vector3D v = plane->GetAxisVector(1);
        mitk::Vector3D n = plane->GetNormal();
        u.Normalize();
        v.Normalize();
        n.Normalize();
        const mitk::Point3D centre = c - n * plane->SignedDistanceFromPlane(c);

        auto figure = Rectangle(plane, centre, a, b);
        mitk::AddPlanarFigureToRoiImage(*figure, *image);

        const double halfThickness = 0.5 * Thickness(image, n);
        unsigned int expected = 0;
        unsigned int missing = 0;
        unsigned int extra = 0;
        for (itk::ImageRegionConstIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
        {
          mitk::Point3D p;
          image->TransformIndexToPhysicalPoint(it.GetIndex(), p);
          const mitk::Vector3D d = p - centre;
          const bool inside = std::abs(d * n) < halfThickness && std::abs(d * u) < a && std::abs(d * v) < b;
          const bool set = it.Get() > 0;
          expected += inside ? 1 : 0;
          missing += inside && !set ? 1 : 0;
          extra += set && !inside ? 1 : 0;
        }
        CPPUNIT_ASSERT_MESSAGE(name + ": the rectangle covers voxels", expected > 0);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": missing voxels (holes)", 0u, missing);
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": voxels outside the rectangle", 0u, extra);
      }
    }
  }

  /** A circle (a polygon of 64 segments) on an axial slice with anisotropic in-plane voxels:
   * every voxel of the slice within the polygon's inscribed radius is set, none beyond the radius,
   * no voxel outside the slice, and the set area differs from pi r^2 by less than one voxel row
   * along the circumference. */
  void Circle_NoHolesAndAreaOfTheCircle()
  {
    auto image = CreateImage(Spacing(1, 1.5, 2), Point(-37.3, 12.1, -5.7));
    mitk::Point3D c;
    auto plane = SlicePlane(image, mitk::AnatomicalPlane::Axial, c);
    mitk::Vector3D n = plane->GetNormal();
    n.Normalize();
    mitk::Vector3D offset;
    offset[0] = 0.3;
    offset[1] = 0.2;
    offset[2] = 0;
    const mitk::Point3D centre = c - n * plane->SignedDistanceFromPlane(c) + offset;
    const double radius = 9.7;

    auto figure = mitk::PlanarCircle::New();
    figure->SetPlaneGeometry(plane);
    mitk::Point2D centre2D;
    mitk::Point2D boundary2D;
    mitk::Vector3D u = plane->GetAxisVector(0);
    u.Normalize();
    plane->Map(centre, centre2D);
    plane->Map(centre + u * radius, boundary2D);
    figure->PlaceFigure(centre2D);
    figure->SetControlPoint(1, boundary2D);

    mitk::AddPlanarFigureToRoiImage(*figure, *image);

    const double inscribed = radius * std::cos(itk::Math::pi / 64);
    unsigned int count = 0;
    for (itk::ImageRegionConstIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      mitk::Point3D p;
      image->TransformIndexToPhysicalPoint(it.GetIndex(), p);
      const mitk::Vector3D d = p - centre;
      const double height = d * n;
      const double inPlane = (d - n * height).GetNorm();
      const bool inSlice = std::abs(height) < 1.0;
      const bool set = it.Get() > 0;
      count += set ? 1 : 0;
      if (inSlice && inPlane < inscribed - 1e-6)
        CPPUNIT_ASSERT_MESSAGE("voxel inside the circle not set", set);
      if (!inSlice || inPlane > radius)
        CPPUNIT_ASSERT_MESSAGE("voxel outside the circle set", !set);
    }
    const double area = count * 1.0 * 1.5;
    CPPUNIT_ASSERT_MESSAGE("area of the circle", std::abs(area - itk::Math::pi * radius * radius) < 2 * itk::Math::pi * radius * 1.5);
  }

  /** A rectangle on a plane oblique to an isotropic image: every voxel that contains a point of
   * the figure at least half a voxel diagonal inside its outline is set (no holes), every set
   * voxel is cut by the plane, and the set voxels, as a layer of the voxel thickness along the
   * normal, have the figure's area within one voxel row along its outline. */
  void ObliquePlane_NoHolesAndAreaOfTheFigure()
  {
    auto image = CreateImage(Spacing(1, 1, 1), Point(-17.3, 2.1, -5.7));
    mitk::Vector3D n;
    n[0] = 1;
    n[1] = 2;
    n[2] = 3;
    n.Normalize();
    mitk::Vector3D u;
    u[0] = 2;
    u[1] = -1;
    u[2] = 0;
    u.Normalize();
    const mitk::Vector3D v = itk::CrossProduct(n, u);
    const mitk::Point3D centre = Point(-17.3 + 19.6, 2.1 + 20.15, -5.7 + 14.4);
    auto plane = ObliquePlane(centre, u, v);
    const double a = 11.3;
    const double b = 7.7;
    auto figure = Rectangle(plane, centre, a, b);

    mitk::AddPlanarFigureToRoiImage(*figure, *image);

    const double margin = 0.5 * std::sqrt(3.0);
    for (double s = -a + margin; s <= a - margin; s += 0.1)
    {
      for (double t = -b + margin; t <= b - margin; t += 0.1)
      {
        itk::Index<3> index;
        CPPUNIT_ASSERT(image->TransformPhysicalPointToIndex(centre + u * s + v * t, index));
        CPPUNIT_ASSERT_MESSAGE("voxel cut by the figure not set", image->GetPixel(index) > 0);
      }
    }

    const double thickness = Thickness(image, n);
    unsigned int count = 0;
    for (itk::ImageRegionConstIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      if (it.Get() == 0)
        continue;
      ++count;
      mitk::Point3D p;
      image->TransformIndexToPhysicalPoint(it.GetIndex(), p);
      CPPUNIT_ASSERT_MESSAGE("set voxel not cut by the plane", std::abs((p - centre) * n) <= 0.5 * thickness);
    }
    const double area = count / thickness;
    CPPUNIT_ASSERT_MESSAGE("area of the figure", std::abs(area - 4 * a * b) < 4 * (a + b));
  }

  /** A plane exactly between two slices selects one of them, not both. */
  void PlaneOnSliceBorder_SetsOneLayer()
  {
    auto image = CreateImage(Spacing(1, 1, 1), Point(0, 0, 0));
    mitk::Vector3D u;
    u.Fill(0);
    u[0] = 1;
    mitk::Vector3D v;
    v.Fill(0);
    v[1] = 1;
    const mitk::Point3D centre = Point(20, 20, 10.5);
    auto plane = ObliquePlane(centre, u, v);
    auto figure = Rectangle(plane, centre, 5.2, 4.2);

    mitk::AddPlanarFigureToRoiImage(*figure, *image);

    std::map<itk::IndexValueType, unsigned int> perSlice;
    for (itk::ImageRegionConstIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      if (it.Get() > 0)
        ++perSlice[it.GetIndex()[2]];
    }
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), perSlice.size());
    CPPUNIT_ASSERT_EQUAL(11u * 9u, perSlice.begin()->second);
  }

  /** A figure without control points, or without a plane, sets nothing. */
  void UnplacedFigure_SetsNothing()
  {
    auto image = CreateImage(Spacing(1, 1, 1), Point(0, 0, 0));
    mitk::Vector3D u;
    u.Fill(0);
    u[0] = 1;
    mitk::Vector3D v;
    v.Fill(0);
    v[1] = 1;
    auto unplaced = mitk::PlanarPolygon::New();
    unplaced->SetPlaneGeometry(ObliquePlane(Point(20, 20, 10), u, v));
    mitk::AddPlanarFigureToRoiImage(*unplaced, *image);
    CPPUNIT_ASSERT_EQUAL(0u, CountSet(image));

    auto withoutPlane = mitk::PlanarPolygon::New();
    mitk::Point2D point;
    point.Fill(0);
    withoutPlane->PlaceFigure(point);
    mitk::AddPlanarFigureToRoiImage(*withoutPlane, *image);
    CPPUNIT_ASSERT_EQUAL(0u, CountSet(image));
  }

  /** A figure reaching beyond the image sets the voxels of its part inside the image. */
  void FigurePartlyOutside_SetsTheVoxelsInside()
  {
    auto image = CreateImage(Spacing(1, 1, 1), Point(0, 0, 0));
    mitk::Vector3D u;
    u.Fill(0);
    u[0] = 1;
    mitk::Vector3D v;
    v.Fill(0);
    v[1] = 1;
    // Rectangle from x = -4.5 to 5.5 and y = 34.5 to 45.5 on slice z = 3: inside the image are
    // x = 0..5 and y = 35..39.
    const mitk::Point3D centre = Point(0.5, 40, 3);
    auto plane = ObliquePlane(centre, u, v);
    auto figure = Rectangle(plane, centre, 5, 5.5);

    mitk::AddPlanarFigureToRoiImage(*figure, *image);

    CPPUNIT_ASSERT_EQUAL(6u * 5u, CountSet(image));
    CPPUNIT_ASSERT(image->GetPixel({{0, 35, 3}}) > 0);
    CPPUNIT_ASSERT(image->GetPixel({{5, 39, 3}}) > 0);
  }

  /** Each figure adds its voxels; voxels set before stay set. */
  void TwoFigures_UnionAndOtherVoxelsKept()
  {
    auto image = CreateImage(Spacing(1, 1, 1), Point(0, 0, 0));
    image->SetPixel({{1, 1, 1}}, 1);
    mitk::Vector3D u;
    u.Fill(0);
    u[0] = 1;
    mitk::Vector3D v;
    v.Fill(0);
    v[1] = 1;
    // Rectangles of 4 x 4 voxels (x = 10..13, 12..15, y = 10..13) on slice z = 5, overlapping in
    // 2 x 4 voxels.
    const mitk::Point3D first = Point(11.5, 11.5, 5);
    const mitk::Point3D second = Point(13.5, 11.5, 5);
    mitk::AddPlanarFigureToRoiImage(*Rectangle(ObliquePlane(first, u, v), first, 2, 2), *image);
    mitk::AddPlanarFigureToRoiImage(*Rectangle(ObliquePlane(second, u, v), second, 2, 2), *image);

    CPPUNIT_ASSERT_EQUAL(1u + 6u * 4u, CountSet(image));
    CPPUNIT_ASSERT(image->GetPixel({{1, 1, 1}}) > 0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberRoiImage)
