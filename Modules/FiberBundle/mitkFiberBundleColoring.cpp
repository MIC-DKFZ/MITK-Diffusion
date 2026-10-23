/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundle.h"
#include "mitkFiberPolyDataBuilder.h"

#include <boost/timer/progress_display.hpp>
#include <itkLinearInterpolateImageFunction.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkLookupTable.h>
#include <vtkCellArray.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkUnsignedCharArray.h>

// Colour codings of mitk::FiberBundle: one RGBA tuple per point.

namespace
{
  /** Opacity in 0..255; values outside (scalar maps beyond 1, NaN from a zero maximum) would be
   *  undefined when cast to unsigned char. */
  unsigned char ClampToByte(double value)
  {
    if (!(value > 0.0))
      return 0;
    if (value >= 255.0)
      return 255;
    return static_cast<unsigned char>(value);
  }
}

void mitk::FiberBundle::ColorFibersByLength(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type)
{
  if (m_MaxFiberLength<=0)
    return;

  auto numOfPoints = this->GetNumberOfPoints();

  //colors and alpha value for each single point, RGBA = 4 components
  unsigned char rgba[4] = {0,0,0,0};
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(numOfPoints * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  auto numOfFibers = m_FiberPolyData->GetNumberOfLines();
  if (numOfFibers < 1)
    return;

  mitk::LookupTable::Pointer mitkLookup = mitk::LookupTable::New();
  mitkLookup->SetType(type);
  if (type!=mitk::LookupTable::MULTILABEL)
    mitkLookup->GetVtkLookupTable()->SetTableRange(m_MinFiberLength, m_MaxFiberLength);

  vtkNew<vtkIdList> idBuffer;
  for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);

    float l = m_FiberLengths.at(i)/m_MaxFiberLength;
    double color[3];
    mitkLookup->GetColor(m_FiberLengths.at(i), color);

    for (const vtkIdType id : ids)
    {

      rgba[0] = static_cast<unsigned char>(255.0 * color[0]);
      rgba[1] = static_cast<unsigned char>(255.0 * color[1]);
      rgba[2] = static_cast<unsigned char>(255.0 * color[2]);
      if (opacity)
        rgba[3] = ClampToByte(255.0f * l);
      else
        rgba[3] = static_cast<unsigned char>(255.0);
      m_FiberColors->InsertTypedTuple(id, rgba);
    }

    if (weight_fibers)
      this->SetFiberWeight(i, m_FiberLengths.at(i));
  }

  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::ColorFibersByOrientation()
{
  vtkPoints* extrPoints = m_FiberPolyData->GetPoints();
  vtkIdType numOfPoints = 0;
  if (extrPoints!=nullptr)
    numOfPoints = extrPoints->GetNumberOfPoints();

  //colors and alpha value for each single point, RGBA = 4 components
  unsigned char rgba[4] = {0,0,0,0};
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");
  // Points on no line (vertices of one-point fibers) are white instead of uninitialised.
  m_FiberColors->SetNumberOfTuples(numOfPoints);
  const unsigned char white[4] = {255, 255, 255, 255};
  for (vtkIdType i=0; i<numOfPoints; ++i)
    m_FiberColors->SetTypedTuple(i, white);

  auto numOfFibers = m_FiberPolyData->GetNumberOfLines();
  if (numOfFibers < 1)
    return;

  /* extract single fibers of fiberBundle */
  vtkCellArray* fiberList = m_FiberPolyData->GetLines();
  fiberList->InitTraversal();
  for (int fi=0; fi<numOfFibers; ++fi) {

    vtkIdType const* idList; // contains the point id's of the line
    vtkIdType pointsPerFiber; // number of points for current line
    fiberList->GetNextCell(pointsPerFiber, idList);

    /* single fiber checkpoints: is number of points valid */
    if (pointsPerFiber > 1)
    {
      /* operate on points of single fiber */
      for (int i=0; i <pointsPerFiber; ++i)
      {
        /* process all points elastV[0]ept starting and endpoint for calculating color value take current point, previous point and next point */
        if (i<pointsPerFiber-1 && i > 0)
        {
          /* The color value of the current point is influenced by the previous point and next point. */
          vnl_vector_fixed< double, 3 > currentPntvtk(extrPoints->GetPoint(idList[i])[0], extrPoints->GetPoint(idList[i])[1],extrPoints->GetPoint(idList[i])[2]);
          vnl_vector_fixed< double, 3 > nextPntvtk(extrPoints->GetPoint(idList[i+1])[0], extrPoints->GetPoint(idList[i+1])[1], extrPoints->GetPoint(idList[i+1])[2]);
          vnl_vector_fixed< double, 3 > prevPntvtk(extrPoints->GetPoint(idList[i-1])[0], extrPoints->GetPoint(idList[i-1])[1], extrPoints->GetPoint(idList[i-1])[2]);

          vnl_vector_fixed< double, 3 > diff1;
          diff1 = currentPntvtk - nextPntvtk;

          vnl_vector_fixed< double, 3 > diff2;
          diff2 = currentPntvtk - prevPntvtk;

          vnl_vector_fixed< double, 3 > diff;
          diff = (diff1 - diff2) / 2.0;
          diff.normalize();

          rgba[0] = static_cast<unsigned char>(255.0 * std::fabs(diff[0]));
          rgba[1] = static_cast<unsigned char>(255.0 * std::fabs(diff[1]));
          rgba[2] = static_cast<unsigned char>(255.0 * std::fabs(diff[2]));
          rgba[3] = static_cast<unsigned char>(255.0);
        }
        else if (i==0)
        {
          /* First point has no previous point, therefore only diff1 is taken */

          vnl_vector_fixed< double, 3 > currentPntvtk(extrPoints->GetPoint(idList[i])[0], extrPoints->GetPoint(idList[i])[1],extrPoints->GetPoint(idList[i])[2]);
          vnl_vector_fixed< double, 3 > nextPntvtk(extrPoints->GetPoint(idList[i+1])[0], extrPoints->GetPoint(idList[i+1])[1], extrPoints->GetPoint(idList[i+1])[2]);

          vnl_vector_fixed< double, 3 > diff1;
          diff1 = currentPntvtk - nextPntvtk;
          diff1.normalize();

          rgba[0] = static_cast<unsigned char>(255.0 * std::fabs(diff1[0]));
          rgba[1] = static_cast<unsigned char>(255.0 * std::fabs(diff1[1]));
          rgba[2] = static_cast<unsigned char>(255.0 * std::fabs(diff1[2]));
          rgba[3] = static_cast<unsigned char>(255.0);
        }
        else if (i==pointsPerFiber-1)
        {
          /* Last point has no next point, therefore only diff2 is taken */
          vnl_vector_fixed< double, 3 > currentPntvtk(extrPoints->GetPoint(idList[i])[0], extrPoints->GetPoint(idList[i])[1],extrPoints->GetPoint(idList[i])[2]);
          vnl_vector_fixed< double, 3 > prevPntvtk(extrPoints->GetPoint(idList[i-1])[0], extrPoints->GetPoint(idList[i-1])[1], extrPoints->GetPoint(idList[i-1])[2]);

          vnl_vector_fixed< double, 3 > diff2;
          diff2 = currentPntvtk - prevPntvtk;
          diff2.normalize();

          rgba[0] = static_cast<unsigned char>(255.0 * std::fabs(diff2[0]));
          rgba[1] = static_cast<unsigned char>(255.0 * std::fabs(diff2[1]));
          rgba[2] = static_cast<unsigned char>(255.0 * std::fabs(diff2[2]));
          rgba[3] = static_cast<unsigned char>(255.0);
        }
        m_FiberColors->InsertTypedTuple(idList[i], rgba);
      }
    }
    else if (pointsPerFiber == 1)
    {
      /* a single point does not define a fiber (use vertex mechanisms instead */
      continue;
    }
    else
    {
      MITK_DEBUG << "Fiber with 0 points detected... please check your tractography algorithm!" ;
      continue;
    }
  }
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::ColorFibersByCurvature(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type)
{
  double window = 5;

  //colors and alpha value for each single point, RGBA = 4 components
  unsigned char rgba[4] = {0,0,0,0};
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(m_FiberPolyData->GetNumberOfPoints() * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  std::vector< double > values;
  double min = 1;
  double max = 0;
  MITK_INFO << "Coloring fibers by curvature";
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));

  vtkNew<vtkIdList> idBuffer;
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    ++disp;
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    const auto numPoints = static_cast<int>(ids.size());
    double mean_curv = 0;

    // calculate curvatures
    for (int j=0; j<numPoints; j++)
    {
      double dist = 0;
      int c = j;
      std::vector< vnl_vector_fixed< double, 3 > > vectors;
      vnl_vector_fixed< double, 3 > meanV; meanV.fill(0.0);
      // backwards down to the fiber's first segment (points c-1 and c)
      while(dist<window/2 && c>0)
      {
        double p1[3];
        m_FiberPolyData->GetPoint(ids[c-1], p1);
        double p2[3];
        m_FiberPolyData->GetPoint(ids[c], p2);

        vnl_vector_fixed< double, 3 > v;
        v[0] = p2[0]-p1[0];
        v[1] = p2[1]-p1[1];
        v[2] = p2[2]-p1[2];
        dist += v.magnitude();
        v.normalize();
        vectors.push_back(v);
        meanV += v;
        c--;
      }
      c = j;
      dist = 0;
      while(dist<window/2 && c<numPoints-1)
      {
        double p1[3];
        m_FiberPolyData->GetPoint(ids[c], p1);
        double p2[3];
        m_FiberPolyData->GetPoint(ids[c+1], p2);

        vnl_vector_fixed< double, 3 > v;
        v[0] = p2[0]-p1[0];
        v[1] = p2[1]-p1[1];
        v[2] = p2[2]-p1[2];
        dist += v.magnitude();
        v.normalize();
        vectors.push_back(v);
        meanV += v;
        c++;
      }
      meanV.normalize();

      double dev = 0;
      for (unsigned int c=0; c<vectors.size(); c++)
      {
        double angle = dot_product(meanV, vectors.at(c));
        if (angle>1.0)
          angle = 1.0;
        if (angle<-1.0)
          angle = -1.0;
        dev += acos(angle)*180/itk::Math::pi;
      }
      if (vectors.size()>0)
        dev /= vectors.size();

      if (weight_fibers)
        mean_curv += dev;
      dev = 1.0-dev/180.0;
      values.push_back(dev);
      if (dev<min)
        min = dev;
      if (dev>max)
        max = dev;
    }

    if (weight_fibers)
      this->SetFiberWeight(i, mean_curv/numPoints);
  }

  mitk::LookupTable::Pointer mitkLookup = mitk::LookupTable::New();
  mitkLookup->SetType(type);
  if (type!=mitk::LookupTable::MULTILABEL)
    mitkLookup->GetVtkLookupTable()->SetTableRange(min, max);

  unsigned int count = 0;
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    for (const vtkIdType id : FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer))
    {
      double color[3];
      double dev = values.at(count);
      mitkLookup->GetColor(dev, color);

      rgba[0] = static_cast<unsigned char>(255.0 * color[0]);
      rgba[1] = static_cast<unsigned char>(255.0 * color[1]);
      rgba[2] = static_cast<unsigned char>(255.0 * color[2]);

      if (opacity)
        rgba[3] = ClampToByte(255.0f * dev/max);
      else
        rgba[3] = static_cast<unsigned char>(255.0);

      m_FiberColors->InsertTypedTuple(id, rgba);
      count++;
    }
  }
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::ColorFibersByScalarMap(const mitk::Image* FAimage, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate)
{
  // Every pixel type is sampled as float: an integer interpolator would truncate interpolated
  // values (a point halfway across a 0/1 mask edge would count as 0). Integer labels below 2^24
  // are exact in float.
  using ImageType = itk::Image<float, 3>;
  ImageType::Pointer itkImage = ImageType::New();
  CastToItkImage(FAimage, itkImage);

  ColorFibersByScalarMap<float>(itkImage.GetPointer(), opacity, weight_fibers, type, max_cap, interpolate );
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

template <typename TPixel>
void mitk::FiberBundle::ColorFibersByScalarMap(const itk::Image<TPixel, 3>* image, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate)
{
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(m_FiberPolyData->GetNumberOfPoints() * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  unsigned char rgba[4] = {0,0,0,0};
  vtkPoints* pointSet = m_FiberPolyData->GetPoints();

  if (type==mitk::LookupTable::MULTILABEL)
    interpolate = false;

  auto interpolator = itk::LinearInterpolateImageFunction< itk::Image<TPixel, 3>, float >::New();
  interpolator->SetInputImage(image);

  double min = 999999;
  double max = -999999;

  vtkNew<vtkIdList> idBuffer;
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    const auto numPoints = static_cast<vtkIdType>(ids.size());
    double mean_val = 0;

    for (const vtkIdType id : ids)
    {
      double p[3];
      m_FiberPolyData->GetPoint(id, p);
      auto pixelValue = mitk::imv::GetImageValue<TPixel>(mitk::imv::GetItkPoint(p), interpolate, interpolator);

      if (pixelValue>max)
        max = pixelValue;
      if (pixelValue<min)
        min = pixelValue;

      if (weight_fibers)
        mean_val += pixelValue;
    }

    if (weight_fibers)
      this->SetFiberWeight(i, mean_val/numPoints);
  }

  mitk::LookupTable::Pointer mitkLookup = mitk::LookupTable::New();
  mitkLookup->SetType(type);
  if (type!=mitk::LookupTable::MULTILABEL && max_cap > 0)
    mitkLookup->GetVtkLookupTable()->SetTableRange(min, max*max_cap);

  for(long i=0; i<m_FiberPolyData->GetNumberOfPoints(); ++i)
  {
    itk::Point<float, 3> px;
    px[0] = pointSet->GetPoint(i)[0];
    px[1] = pointSet->GetPoint(i)[1];
    px[2] = pointSet->GetPoint(i)[2];
    auto pixelValue = mitk::imv::GetImageValue<TPixel>(px, interpolate, interpolator);

    double color[3];
    mitkLookup->GetColor(pixelValue, color);

    rgba[0] = static_cast<unsigned char>(255.0 * color[0]);
    rgba[1] = static_cast<unsigned char>(255.0 * color[1]);
    rgba[2] = static_cast<unsigned char>(255.0 * color[2]);
    if (opacity)
      rgba[3] = ClampToByte(255.0 * pixelValue);
    else
      rgba[3] = static_cast<unsigned char>(255.0);
    m_FiberColors->InsertTypedTuple(i, rgba);
  }
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}


void mitk::FiberBundle::ColorFibersByFiberWeights(bool opacity, mitk::LookupTable::LookupTableType type)
{
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(m_FiberPolyData->GetNumberOfPoints() * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  unsigned char rgba[4] = {0,0,0,0};
  unsigned int counter = 0;

  float max = -999999;
  float min = 999999;
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    float weight = this->GetFiberWeight(i);
    if (weight>max)
      max = weight;
    if (weight<min)
      min = weight;
  }
  if (fabs(max-min)<0.00001f)
  {
    max = 1;
    min = 0;
  }

  mitk::LookupTable::Pointer mitkLookup = mitk::LookupTable::New();
  mitkLookup->SetType(type);
  if (type!=mitk::LookupTable::MULTILABEL)
    mitkLookup->GetVtkLookupTable()->SetTableRange(min, max);

  vtkNew<vtkIdList> idBuffer;
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    const auto numPoints = static_cast<vtkIdType>(FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer).size());
    auto weight = this->GetFiberWeight(i);

    double color[3];
    mitkLookup->GetColor(weight, color);

    for (int j=0; j<numPoints; j++)
    {
      rgba[0] = static_cast<unsigned char>(255.0 * color[0]);
      rgba[1] = static_cast<unsigned char>(255.0 * color[1]);
      rgba[2] = static_cast<unsigned char>(255.0 * color[2]);
      if (opacity)
        rgba[3] = ClampToByte(255.0f * weight/max);
      else
        rgba[3] = static_cast<unsigned char>(255.0);

      m_FiberColors->InsertTypedTuple(counter, rgba);
      counter++;
    }
  }

  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::SetSingleFiberColor(float r, float g, float b, unsigned int cellId, float alpha)
{
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(m_FiberPolyData->GetNumberOfPoints() * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  // the given colour for fiber cellId, white for all others
  const unsigned char picked[4] = {static_cast<unsigned char>(r), static_cast<unsigned char>(g), static_cast<unsigned char>(b), static_cast<unsigned char>(alpha)};
  const unsigned char other[4] = {255, 255, 255, static_cast<unsigned char>(alpha)};
  unsigned int counter = 0;
  vtkNew<vtkIdList> idBuffer;
  for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    const auto numPoints = static_cast<vtkIdType>(FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer).size());
    for (vtkIdType j=0; j<numPoints; j++)
      m_FiberColors->InsertTypedTuple(counter++, i==cellId ? picked : other);
  }

  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::SetFiberColors(float r, float g, float b, float alpha)
{
  m_FiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  m_FiberColors->ReserveValues(m_FiberPolyData->GetNumberOfPoints() * 4);
  m_FiberColors->SetNumberOfComponents(4);
  m_FiberColors->SetName("FIBER_COLORS");

  unsigned char rgba[4] = {0,0,0,0};
  for(long i=0; i<m_FiberPolyData->GetNumberOfPoints(); ++i)
  {
    rgba[0] = static_cast<unsigned char>(r);
    rgba[1] = static_cast<unsigned char>(g);
    rgba[2] = static_cast<unsigned char>(b);
    rgba[3] = static_cast<unsigned char>(alpha);
    m_FiberColors->InsertTypedTuple(i, rgba);
  }
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}

void mitk::FiberBundle::SetFiberColors(vtkSmartPointer<vtkUnsignedCharArray> fiberColors)
{
  // One RGBA or RGB tuple per point; RGB colours get an opaque alpha.
  const int components = fiberColors==nullptr ? 0 : fiberColors->GetNumberOfComponents();
  if (fiberColors==nullptr || (components!=3 && components!=4) || fiberColors->GetNumberOfTuples()<m_FiberPolyData->GetNumberOfPoints())
    mitkThrow() << "Fiber colors need one RGB or RGBA tuple per point (" << m_FiberPolyData->GetNumberOfPoints() << " points).";

  for(long i=0; i<m_FiberPolyData->GetNumberOfPoints(); ++i)
  {
    unsigned char target[4] = {0,0,0,255};
    for (int c=0; c<components; ++c)
      target[c] = fiberColors->GetValue(i*components+c);
    m_FiberColors->InsertTypedTuple(i, target);
  }
  m_UpdateTime3D.Modified();
  m_UpdateTime2D.Modified();
}
