/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

/** A chain of FiberBundle transformations on fiberBundleX.fib against the MITK-Data reference. */
class mitkFiberTransformationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberTransformationTestSuite);
  MITK_TEST(TransformationChain_EqualsReference);
  CPPUNIT_TEST_SUITE_END();

public:
  void TransformationChain_EqualsReference()
  {
    auto fib = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/fiberBundleX.fib"));
    auto reference = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/fiberBundleX_transformed.fib"));

    fib->RotateAroundAxis(90, 45, 10);
    fib->TranslateFibers(2, 3, 5);
    fib->ScaleFibers(1, 0.1, 1.3);
    fib->RemoveLongFibers(150);
    fib->RemoveShortFibers(20);
    fib->ResampleSpline(1.0);
    fib->ApplyCurvatureThreshold(3.0, true);
    fib->MirrorFibers(0);
    fib->MirrorFibers(1);
    fib->MirrorFibers(2);

    CPPUNIT_ASSERT_MESSAGE("Transformed bundle equals the reference", fib->Equals(reference));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberTransformation)
