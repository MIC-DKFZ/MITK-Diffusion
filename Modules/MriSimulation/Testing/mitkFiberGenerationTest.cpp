/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkTestingMacros.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkPlanarEllipse.h>
#include <itkFibersFromPlanarFiguresFilter.h>
#include <omp.h>
#include <vtkCell.h>
#include <vtkPolyData.h>

namespace
{
  /** Largest distance of the bundle's first (last) fiber points from the centre of the first (last)
   *  fiducial, relative to the ellipse: at most 1 when every fiber starts and ends inside them. */
  double MaxEllipseRadiusOfEnds(mitk::FiberBundle* bundle, mitk::PlanarEllipse* first, mitk::PlanarEllipse* last)
  {
    double maxRadius = 0;
    for (unsigned int i = 0; i < bundle->GetNumFibers(); ++i)
    {
      vtkCell* cell = bundle->GetFiberPolyData()->GetCell(i);
      for (auto end : { std::make_pair(first, 0), std::make_pair(last, static_cast<int>(cell->GetNumberOfPoints()) - 1) })
      {
        mitk::Point3D world(cell->GetPoints()->GetPoint(end.second));
        mitk::Point2D point;
        end.first->GetPlaneGeometry()->Map(world, point);
        const mitk::Point2D centre = end.first->GetControlPoint(0);
        const mitk::Vector2D axis1 = end.first->GetControlPoint(1) - centre;
        const mitk::Vector2D axis2 = end.first->GetControlPoint(2) - centre;
        const mitk::Vector2D d = point - centre;
        const double u = (d * axis1) / axis1.GetSquaredNorm();
        const double v = (d * axis2) / axis2.GetSquaredNorm();
        maxRadius = std::max(maxRadius, std::sqrt(u * u + v * v));
      }
    }
    MITK_INFO << "Largest relative ellipse radius of the fiber ends: " << maxRadius;
    return maxRadius;
  }
}

/**Documentation
 *  Test if fiber transfortaiom methods work correctly
 */
int mitkFiberGenerationTest(int argc, char* argv[])
{
    omp_set_num_threads(1);
    MITK_TEST_BEGIN("mitkFiberGenerationTest");

    MITK_TEST_CONDITION_REQUIRED(argc==6,"check for input data")

            try{

        mitk::PlanarEllipse::Pointer pf1 = mitk::IOUtil::Load<mitk::PlanarEllipse>(argv[1]);
        mitk::PlanarEllipse::Pointer pf2 = mitk::IOUtil::Load<mitk::PlanarEllipse>(argv[2]);
        mitk::PlanarEllipse::Pointer pf3 = mitk::IOUtil::Load<mitk::PlanarEllipse>(argv[3]);
        mitk::FiberBundle::Pointer uniform = mitk::IOUtil::Load<mitk::FiberBundle>(argv[4]);
        mitk::FiberBundle::Pointer gaussian = mitk::IOUtil::Load<mitk::FiberBundle>(argv[5]);

        mitk::FiberGenerationParameters parameters;
        std::vector< mitk::PlanarEllipse::Pointer > fid; fid.push_back(pf1); fid.push_back(pf2); fid.push_back(pf3);
        std::vector< unsigned int > flip; flip.push_back(0); flip.push_back(0); flip.push_back(0);
        parameters.m_Fiducials.push_back(fid); parameters.m_FlipList.push_back(flip);
        parameters.m_Density = 50;
        parameters.m_Tension = 0;
        parameters.m_Continuity = 0;
        parameters.m_Bias = 0;
        parameters.m_Sampling = 1;
        parameters.m_Variance = 0.1;

        // check uniform fiber distribution
        {
            itk::FibersFromPlanarFiguresFilter::Pointer filter = itk::FibersFromPlanarFiguresFilter::New();
            parameters.m_Distribution = mitk::FiberGenerationParameters::DISTRIBUTE_UNIFORM;
            filter->SetParameters(parameters);
            filter->SetFixSeed(true);
            filter->Update();
            std::vector< mitk::FiberBundle::Pointer > fiberBundles = filter->GetFiberBundles();
            bool equal = uniform->Equals(fiberBundles.at(0));
            if (!equal)
            {
              std::string out = mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + "uniform.fib";
              MITK_INFO << "Saving test bundle to " << out;
              mitk::IOUtil::Save(fiberBundles.at(0), out);
            }
            MITK_TEST_CONDITION_REQUIRED(equal,"check uniform bundle")
            MITK_TEST_CONDITION_REQUIRED(MaxEllipseRadiusOfEnds(fiberBundles.at(0), pf1, pf3) <= 1.0 + 1e-6, "uniform fibers start and end inside the fiducials")
        }

        // check gaussian fiber distribution
        {
            itk::FibersFromPlanarFiguresFilter::Pointer filter = itk::FibersFromPlanarFiguresFilter::New();
            parameters.m_Distribution = mitk::FiberGenerationParameters::DISTRIBUTE_GAUSSIAN;
            filter->SetParameters(parameters);
            filter->SetFixSeed(true);
            filter->Update();
            std::vector< mitk::FiberBundle::Pointer > fiberBundles = filter->GetFiberBundles();
            bool equal = gaussian->Equals(fiberBundles.at(0));
            if (!equal)
            {
              std::string out = mitk::IOUtil::GetTempPath() + mitk::IOUtil::GetDirectorySeparator() + "gaussian.fib";
              MITK_INFO << "Saving test bundle to " << out;
              mitk::IOUtil::Save(fiberBundles.at(0), out);
            }
            MITK_TEST_CONDITION_REQUIRED(equal,"check gaussian bundle")
            MITK_TEST_CONDITION_REQUIRED(MaxEllipseRadiusOfEnds(fiberBundles.at(0), pf1, pf3) <= 1.0 + 1e-6, "gaussian fibers start and end inside the fiducials")
        }
    }
    catch(...) {
        return EXIT_FAILURE;
    }

    // always end with this!
    MITK_TEST_END();
}
