/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberfoxMotion_h
#define mitkFiberfoxMotion_h

#include <itkMatrix.h>
#include <itkPoint.h>
#include <itkScalableAffineTransform.h>
#include <itkVector.h>
#include <mitkDiffusionImageHelperFunctions.h>

namespace mitk
{
  /**
  * \brief Rigid head motion of one simulated volume: x' = R (x - c) + c + t.
  *
  * R = Rz Ry Rx of the rotation angles in degrees (mitk::imv::GetRotationMatrixItk), t the
  * translation in mm and c a centre that is the same for all volumes of a simulation (the centre of
  * the unmoved fiber bundle). The simulation moves the fibers and the tissue mask forward and reads
  * the maps defined on the unmoved head (field map, volume fractions) through the exact inverse
  * x = R^T (x' - c - t) + c, so fibers, mask and maps stay aligned in every volume.
  */
  class FiberfoxMotion
  {
  public:
    using VectorType = itk::Vector<double, 3>;
    using PointType = itk::Point<double, 3>;
    using MatrixType = itk::Matrix<double, 3, 3>;
    using TransformType = itk::ScalableAffineTransform<double, 3>;

    /** No motion. */
    FiberfoxMotion()
    {
      m_Rotation.Fill(0.0);
      m_Translation.Fill(0.0);
      m_Center.Fill(0.0);
      m_Matrix.SetIdentity();
      m_InverseMatrix.SetIdentity();
    }

    FiberfoxMotion(const VectorType& rotation, const VectorType& translation, const PointType& center)
      : m_Rotation(rotation)
      , m_Translation(translation)
      , m_Center(center)
      , m_Matrix(mitk::imv::GetRotationMatrixItk<double>(rotation[0], rotation[1], rotation[2]))
      , m_InverseMatrix(m_Matrix.GetTranspose())
    {
    }

    /** Position of an unmoved point after the motion. */
    template <class TPoint>
    TPoint Forward(const TPoint& point) const
    {
      VectorType d;
      for (unsigned int i = 0; i < 3; ++i)
        d[i] = point[i] - m_Center[i];
      d = m_Matrix * d;
      TPoint moved;
      for (unsigned int i = 0; i < 3; ++i)
        moved[i] = d[i] + m_Center[i] + m_Translation[i];
      return moved;
    }

    /** Unmoved position of a point of the moved head. */
    template <class TPoint>
    TPoint Inverse(const TPoint& point) const
    {
      VectorType d;
      for (unsigned int i = 0; i < 3; ++i)
        d[i] = point[i] - m_Center[i] - m_Translation[i];
      d = m_InverseMatrix * d;
      TPoint unmoved;
      for (unsigned int i = 0; i < 3; ++i)
        unmoved[i] = d[i] + m_Center[i];
      return unmoved;
    }

    /** The forward motion as a transform, e.g. for mitk::FiberBundle::TransformFibers. */
    TransformType::Pointer GetTransform() const
    {
      auto transform = TransformType::New();
      transform->SetCenter(m_Center);
      transform->SetMatrix(m_Matrix);
      transform->SetTranslation(m_Translation);
      return transform;
    }

    bool IsIdentity() const
    {
      return m_Rotation.GetSquaredNorm() == 0 && m_Translation.GetSquaredNorm() == 0;
    }

    const VectorType& GetRotation() const { return m_Rotation; }        ///< degrees about x, y and z
    const VectorType& GetTranslation() const { return m_Translation; }  ///< mm
    const PointType& GetCenter() const { return m_Center; }
    const MatrixType& GetMatrix() const { return m_Matrix; }

  private:
    VectorType m_Rotation;
    VectorType m_Translation;
    PointType m_Center;
    MatrixType m_Matrix;
    MatrixType m_InverseMatrix;
  };
}

#endif
