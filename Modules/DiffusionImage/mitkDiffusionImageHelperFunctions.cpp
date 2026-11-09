/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionImageHelperFunctions.h"

#include "itkVectorContainer.h"
#include <vtkBox.h>
#include <itksys/SystemTools.hxx>
#include <charconv>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <sstream>
#include <vnl/vnl_det.h>

// Intersect a finite line (with end points p0 and p1) with all of the
// cells of a vtkImageData
std::vector< std::pair< itk::Index<3>, double > > mitk::imv::IntersectImage(const itk::Vector<double,3>& spacing, itk::Index<3> si, itk::Index<3> ei, const itk::ContinuousIndex<float, 3>& sf, const itk::ContinuousIndex<float, 3>& ef)
{
  std::vector< std::pair< itk::Index<3>, double > > out;
  if (si == ei)
  {
    double d[3];
    for (int i=0; i<3; ++i)
      d[i] = static_cast<double>(sf[i]-ef[i])*spacing[i];
    double len = std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
    out.push_back(  std::pair< itk::Index<3>, double >(si, len) );
    return out;
  }

  double bounds[6];

  double entrancePoint[3];
  double exitPoint[3];

  double startPoint[3];
  double endPoint[3];

  double t0, t1;
  for (unsigned int i=0; i<3; ++i)
  {
    startPoint[i] = static_cast<double>(sf[i]);
    endPoint[i] = static_cast<double>(ef[i]);

    if (si[i]>ei[i])
    {
      auto t = si[i];
      si[i] = ei[i];
      ei[i] = t;
    }
  }

  for (auto x = si[0]; x<=ei[0]; ++x)
    for (auto y = si[1]; y<=ei[1]; ++y)
      for (auto z = si[2]; z<=ei[2]; ++z)
      {
        bounds[0] = static_cast<double>(x) - 0.5;
        bounds[1] = static_cast<double>(x) + 0.5;
        bounds[2] = static_cast<double>(y) - 0.5;
        bounds[3] = static_cast<double>(y) + 0.5;
        bounds[4] = static_cast<double>(z) - 0.5;
        bounds[5] = static_cast<double>(z) + 0.5;

        int entryPlane;
        int exitPlane;
        int hit = vtkBox::IntersectWithLine(bounds,
                                            startPoint,
                                            endPoint,
                                            t0,
                                            t1,
                                            entrancePoint,
                                            exitPoint,
                                            entryPlane,
                                            exitPlane);
        if (hit)
        {
          if (entryPlane>=0 && exitPlane>=0)
          {
            double d[3];
            for (int i=0; i<3; ++i)
              d[i] = (exitPoint[i] - entrancePoint[i])*spacing[i];
            double len = std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );

            itk::Index<3> idx; idx[0] = x; idx[1] = y; idx[2] = z;
            out.push_back(  std::pair< itk::Index<3>, double >(idx, len) );
          }
          else if (entryPlane>=0)
          {
            double d[3];
            for (int i=0; i<3; ++i)
              d[i] = (static_cast<double>(ef[i]) - entrancePoint[i])*spacing[i];
            double len = std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );

            itk::Index<3> idx; idx[0] = x; idx[1] = y; idx[2] = z;
            out.push_back(  std::pair< itk::Index<3>, double >(idx, len) );
          }
          else if (exitPlane>=0)
          {
            double d[3];
            for (int i=0; i<3; ++i)
              d[i] = (exitPoint[i]-static_cast<double>(sf[i]))*spacing[i];
            double len = std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );

            itk::Index<3> idx; idx[0] = x; idx[1] = y; idx[2] = z;
            out.push_back(  std::pair< itk::Index<3>, double >(idx, len) );
          }
        }
      }
  return out;
}

//------------------------- gradients-function ------------------------------------


namespace
{
  /** The numbers of each non-empty line of a text file. */
  std::vector<std::vector<double>> ReadNumberRows(const std::string& file, const std::string& what)
  {
    if (!itksys::SystemTools::FileExists(file))
      mitkThrow() << what << " file not existing: " << file;
    std::ifstream stream(file);
    if (!stream.is_open())
      mitkThrow() << what << " file could not be opened: " << file;

    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(stream, line))
    {
      std::vector<double> row;
      std::istringstream tokens(line);
      std::string token;
      while (tokens >> token)
      {
        double value = 0;
        const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
        if (error != std::errc() || end != token.data() + token.size())
          mitkThrow() << "Encountered invalid " << what << " file entry >" << token << "< in " << file;
        row.push_back(value);
      }
      if (!row.empty())
        rows.push_back(row);
    }
    return rows;
  }
}

mitk::gradients::GradientDirectionContainerType::ConstPointer mitk::gradients::ReadBvalsBvecs(std::string bvals_file, std::string bvecs_file, double& reference_bval)
{
  std::vector<double> bval_entries;
  for (const auto& row : ReadNumberRows(bvals_file, "bvals"))
    bval_entries.insert(bval_entries.end(), row.begin(), row.end());
  const auto bvec_rows = ReadNumberRows(bvecs_file, "bvecs");

  const auto n = bval_entries.size();
  // FSL writes three rows of n values. A file with one gradient per line is rejected rather than
  // transposed silently, since its frame is not known.
  if (bvec_rows.size() != 3 || bvec_rows[0].size() != n || bvec_rows[1].size() != n || bvec_rows[2].size() != n)
  {
    mitkThrow() << "bvecs file " << bvecs_file << " has " << bvec_rows.size() << " rows; expected 3 rows of "
                << n << " values, one per b-value in " << bvals_file;
  }

  reference_bval = -1;
  for (auto b : bval_entries)
    reference_bval = std::max(reference_bval, b);

  auto directioncontainer = mitk::gradients::GradientDirectionContainerType::New();
  for (unsigned int i = 0; i < n; i++)
  {
    double b_val = bval_entries.at(i);

    mitk::gradients::GradientDirectionType vec;
    vec[0] = bvec_rows[0][i];
    vec[1] = bvec_rows[1][i];
    vec[2] = bvec_rows[2][i];

    if (b_val>0 && vec.magnitude() < 0.001) // could be IVIM data --> set vec to  1,0,0
    {
        vec[0] = 1.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    }

    // Adjust the vector length to encode gradient strength
    if (reference_bval>0)
    {
      double factor = b_val/reference_bval;
      if(vec.magnitude() > 0)
      {
        vec.normalize();
        vec[0] = sqrt(factor)*vec[0];
        vec[1] = sqrt(factor)*vec[1];
        vec[2] = sqrt(factor)*vec[2];
      }
    }

    directioncontainer->InsertElement(i,vec);
  }

  return GradientDirectionContainerType::ConstPointer(directioncontainer);
}

mitk::DiffusionPropertyHelper::MeasurementFrameType mitk::gradients::FslBvecsMeasurementFrame(const itk::Matrix<double, 3, 3>& direction)
{
  // FSL gives bvecs in the voxel frame, with x negated when the voxel-to-world matrix has a positive
  // determinant (RAS and LPS matrices have the same determinant). The world vector is then
  // direction * flip * bvec; measurement frames hold the frame axes as rows.
  vnl_matrix_fixed<double, 3, 3> toWorld = direction.GetVnlMatrix();
  if (vnl_det(toWorld) > 0)
    toWorld.set_column(0, -toWorld.get_column(0));
  return toWorld.transpose();
}

void mitk::gradients::WriteBvalsBvecs(std::string bvals_file, std::string bvecs_file, GradientDirectionContainerType::ConstPointer gradients, double reference_bval)
{
  // The bvecs are normalized, so the b-value is all that is left of the gradient length: it is
  // written as b_ref |g|^2, rounded to 1e-6 s/mm^2 and in its shortest form ("1000", "999.999",
  // "162.5"), not as the integer key of the b-value map, which changed b-values by up to 0.5.
  std::ofstream myfile;
  myfile.open (bvals_file.c_str());
  for(unsigned int i=0; i<gradients->Size(); i++)
  {
    const double b = reference_bval * gradients->ElementAt(i).squared_magnitude();
    myfile << std::format("{}", std::round(b * 1e6) / 1e6) << " ";
  }
  myfile.close();

  std::ofstream myfile2;
  myfile2.open (bvecs_file.c_str());
  myfile2.precision(std::numeric_limits<double>::max_digits10);
  for(int j=0; j<3; j++)
  {
    for(unsigned int i=0; i<gradients->Size(); i++)
    {
      GradientDirectionType direction = gradients->ElementAt(i);
      direction.normalize();
      myfile2 << direction.get(j) << " ";
    }
    myfile2 << std::endl;
  }
}

std::vector<unsigned int> mitk::gradients::GetAllUniqueDirections(const BValueMap & refBValueMap, GradientDirectionContainerType::ConstPointer refGradientsContainer )
{

  IndiciesVector directioncontainer;
  auto mapIterator = refBValueMap.begin();

  if(refBValueMap.find(0) != refBValueMap.end() && refBValueMap.size() > 1)
    mapIterator++; //skip bzero Values

  for( ; mapIterator != refBValueMap.end(); mapIterator++){

    IndiciesVector currentShell = mapIterator->second;

    while(currentShell.size()>0)
    {
      unsigned int wntIndex = currentShell.back();
      currentShell.pop_back();

      auto containerIt = directioncontainer.begin();
      bool directionExist = false;
      while(containerIt != directioncontainer.end())
      {
        // The gradient length encodes the b-value; only the directions are compared.
        if (fabs(dot_product(GradientDirectionType(refGradientsContainer->ElementAt(*containerIt)).normalize(),
                             GradientDirectionType(refGradientsContainer->ElementAt(wntIndex)).normalize())) > 0.9998)
        {
          directionExist = true;
          break;
        }
        containerIt++;
      }
      if(!directionExist)
      {
        directioncontainer.push_back(wntIndex);
      }
    }
  }

  return directioncontainer;
}


void mitk::gradients::Cart2Sph(double x, double y, double z, double *spherical)
{
    double phi, th, rad;
    rad = sqrt(x*x+y*y+z*z);
    if( rad < mitk::eps )
    {
        th = itk::Math::pi/2;
        phi = itk::Math::pi/2;
    }
    else
    {
        th = acos(z/rad);
        phi = atan2(y, x);
    }
    spherical[0] = phi;
    spherical[1] = th;
    spherical[2] = rad;
}


vnl_vector_fixed<double, 3> mitk::gradients::Sph2Cart(const double& theta, const double& phi, const double& rad)
{
    vnl_vector_fixed<double, 3> dir;
    dir[0] = rad * sin(theta) * cos(phi);
    dir[1] = rad * sin(theta) * sin(phi);
    dir[2] = rad * cos(theta);
    return dir;
}

vnl_matrix<double> mitk::gradients::ComputeSphericalFromCartesian(const IndiciesVector & refShell, const GradientDirectionContainerType * refGradientsContainer)
{

  vnl_matrix<double> Q(3, refShell.size());
  Q.fill(0.0);

  for(unsigned int i = 0; i < refShell.size(); i++)
  {
    GradientDirectionType dir = refGradientsContainer->ElementAt(refShell[i]);
    double x = dir.normalize().get(0);
    double y = dir.normalize().get(1);
    double z = dir.normalize().get(2);
    double cart[3];
    mitk::gradients::Cart2Sph(x,y,z,cart);
    Q(0,i) = cart[0];
    Q(1,i) = cart[1];
    Q(2,i) = cart[2];
  }
  return Q;
}

mitk::gradients::GradientDirectionContainerType::Pointer mitk::gradients::CreateNormalizedUniqueGradientDirectionContainer(const mitk::gradients::BValueMap & bValueMap,
                                                                                                                           const GradientDirectionContainerType *origninalGradentcontainer)
{
  mitk::gradients::GradientDirectionContainerType::Pointer directioncontainer = mitk::gradients::GradientDirectionContainerType::New();
  auto mapIterator = bValueMap.begin();

  if(bValueMap.find(0) != bValueMap.end() && bValueMap.size() > 1){
    mapIterator++; //skip bzero Values
    vnl_vector_fixed<double, 3> vec;
    vec.fill(0.0);
    directioncontainer->push_back(vec);
  }

  for( ; mapIterator != bValueMap.end(); mapIterator++){

    IndiciesVector currentShell = mapIterator->second;

    while(currentShell.size()>0)
    {
      unsigned int wntIndex = currentShell.back();
      currentShell.pop_back();

      mitk::gradients::GradientDirectionContainerType::Iterator containerIt = directioncontainer->Begin();
      bool directionExist = false;
      while(containerIt != directioncontainer->End())
      {
        // The gradient length encodes the b-value; only the directions are compared.
        if (fabs(dot_product(containerIt.Value(), GradientDirectionType(origninalGradentcontainer->ElementAt(wntIndex)).normalize())) > 0.9998)
        {
          directionExist = true;
          break;
        }
        containerIt++;
      }
      if(!directionExist)
      {
        GradientDirectionType dir(origninalGradentcontainer->ElementAt(wntIndex));
        directioncontainer->push_back(dir.normalize());
      }
    }
  }

  return directioncontainer;
}
