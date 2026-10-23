/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkImageGeometry4D_h
#define mitkImageGeometry4D_h

#include <itkImageBase.h>

namespace mitk
{
  /**
   * SH coefficient and peak images exist in two layouts: a 3D image with one vector per voxel and a
   * 4D image with one volume per component (the file layout, and mitk::PeakImage). These functions
   * carry the geometry from one layout to the other in double precision. The fourth axis has
   * spacing 1, origin 0 and no rotation. They set spacing, origin, direction and regions; the
   * caller allocates.
   */

  /** Geometry of image3 plus a fourth axis with numberOfVolumes entries. */
  inline void CopyGeometryTo4D(const itk::ImageBase<3> *image3, itk::ImageBase<4> *image4, itk::SizeValueType numberOfVolumes)
  {
    const auto &region3 = image3->GetLargestPossibleRegion();
    itk::ImageBase<4>::RegionType region4;
    itk::ImageBase<4>::SpacingType spacing4;
    itk::ImageBase<4>::PointType origin4;
    itk::ImageBase<4>::DirectionType direction4;
    spacing4.Fill(1);
    origin4.Fill(0);
    direction4.SetIdentity();
    for (unsigned int r = 0; r < 3; ++r)
    {
      region4.SetIndex(r, region3.GetIndex(r));
      region4.SetSize(r, region3.GetSize(r));
      spacing4[r] = image3->GetSpacing()[r];
      origin4[r] = image3->GetOrigin()[r];
      for (unsigned int c = 0; c < 3; ++c)
        direction4[r][c] = image3->GetDirection()[r][c];
    }
    region4.SetIndex(3, 0);
    region4.SetSize(3, numberOfVolumes);

    image4->SetSpacing(spacing4);
    image4->SetOrigin(origin4);
    image4->SetDirection(direction4);
    image4->SetRegions(region4);
  }

  /** Geometry of the first three axes of image4. */
  inline void CopyGeometryTo3D(const itk::ImageBase<4> *image4, itk::ImageBase<3> *image3)
  {
    const auto &region4 = image4->GetLargestPossibleRegion();
    itk::ImageBase<3>::RegionType region3;
    itk::ImageBase<3>::SpacingType spacing3;
    itk::ImageBase<3>::PointType origin3;
    itk::ImageBase<3>::DirectionType direction3;
    for (unsigned int r = 0; r < 3; ++r)
    {
      region3.SetIndex(r, region4.GetIndex(r));
      region3.SetSize(r, region4.GetSize(r));
      spacing3[r] = image4->GetSpacing()[r];
      origin3[r] = image4->GetOrigin()[r];
      for (unsigned int c = 0; c < 3; ++c)
        direction3[r][c] = image4->GetDirection()[r][c];
    }

    image3->SetSpacing(spacing3);
    image3->SetOrigin(origin3);
    image3->SetDirection(direction3);
    image3->SetRegions(region3);
  }
}

#endif
