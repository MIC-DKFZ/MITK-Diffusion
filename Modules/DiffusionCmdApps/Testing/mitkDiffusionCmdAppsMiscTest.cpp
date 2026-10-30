/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkImageReadAccessor.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>

#include <vnl/vnl_matrix_fixed.h>

using namespace mitk::CmdAppTest;

/** Smoke tests of the apps in Modules/DiffusionCmdApps/Misc; FiberScreenshot renders and has a test of its own. */
class mitkDiffusionCmdAppsMiscTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsMiscTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(FileFormatConverter_DwiToNifti_KeepsGradients);
  MITK_TEST(ImageResampler_Spacing);
  MITK_TEST(CopyGeometry_RotatedReference_TurnsTensors);
  MITK_TEST(CopyGeometry_UnknownInputType_Fails);
  MITK_TEST(CopyGeometry_InputTypeMismatch_Fails);
  MITK_TEST(DiffusionDICOMLoader_MissingFolder_Fails);
  MITK_TEST(ShToOdfImage_Odfs);
  MITK_TEST(PeakExtraction_Peaks);
  MITK_TEST(DReg_SameImage_KeepsGrid);
  CPPUNIT_TEST_SUITE_END();

public:
  void HelpAndXml()
  {
    CheckHelpAndXml("Misc", {"FileFormatConverter", "ImageResampler", "CopyGeometry", "DiffusionDICOMLoader",
                             "ShToOdfImage", "DReg", "FiberScreenshot", "PeakExtraction"});
  }

  void FileFormatConverter_DwiToNifti_KeepsGradients()
  {
    const auto dir = WorkDir("FileFormatConverter");
    RunOk("FileFormatConverter", {"-i", Data("dwi_voxel_b0x3_b1000x192.dwi"), "-o", "dwi.01.nii.gz"}, dir);
    auto dwi = Load<mitk::Image>(dir, "dwi.01.nii.gz");
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    CPPUNIT_ASSERT_EQUAL(195u, static_cast<unsigned int>(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi)->Size()));
  }

  void ImageResampler_Spacing()
  {
    const auto dir = WorkDir("ImageResampler");
    RunOk("ImageResampler", {"-i", Data("ROIIMAGE.nrrd"), "-o", "resampled.nrrd", "-s", "2,2,2"}, dir);
    auto image = Load<mitk::Image>(dir, "resampled.nrrd");
    for (unsigned int d = 0; d < 3; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, image->GetGeometry()->GetSpacing()[d], 1e-6);
  }

  // The reference has the voxel grid of the tensor image turned by 90 degrees about z: the tensors
  // turn with it (#22).
  void CopyGeometry_RotatedReference_TurnsTensors()
  {
    const auto dir = WorkDir("CopyGeometry");
    auto input = mitk::IOUtil::Load<mitk::TensorImage>(Data("tensorImage.dti"));

    auto reference = mitk::Image::New();
    reference->Initialize(mitk::MakeScalarPixelType<float>(), *input->GetGeometry());
    mitk::AffineTransform3D::MatrixType matrix;
    matrix.Fill(0);
    matrix[0][1] = -3;
    matrix[1][0] = 3;
    matrix[2][2] = 3;
    auto transform = mitk::AffineTransform3D::New();
    transform->SetMatrix(matrix);
    reference->GetGeometry()->SetIndexToWorldTransform(transform);
    mitk::IOUtil::Save(reference, Path(dir, "reference.nrrd"));

    RunOk("CopyGeometry", {"-i", Data("tensorImage.dti"), "-r", "reference.nrrd", "-o", "turned.dti"}, dir);
    auto output = Load<mitk::TensorImage>(dir, "turned.dti");

    vnl_matrix_fixed<double, 3, 3> rotation(0.0);
    rotation(0, 1) = -1;
    rotation(1, 0) = 1;
    rotation(2, 2) = 1;
    mitk::ImageReadAccessor inAccessor(input.GetPointer());
    mitk::ImageReadAccessor outAccessor(output.GetPointer());
    const auto *in = static_cast<const float *>(inAccessor.GetData());
    const auto *out = static_cast<const float *>(outAccessor.GetData());
    const std::size_t voxels = static_cast<std::size_t>(input->GetDimension(0)) * input->GetDimension(1) * input->GetDimension(2);
    double maxDifference = 0;
    for (std::size_t v = 0; v < voxels; ++v)
    {
      const float *t = in + 6 * v;
      vnl_matrix_fixed<double, 3, 3> d;
      d(0, 0) = t[0]; d(0, 1) = t[1]; d(0, 2) = t[2];
      d(1, 0) = t[1]; d(1, 1) = t[3]; d(1, 2) = t[4];
      d(2, 0) = t[2]; d(2, 1) = t[4]; d(2, 2) = t[5];
      const auto expected = rotation * d * rotation.transpose();
      const float *o = out + 6 * v;
      const double values[6] = {expected(0, 0), expected(0, 1), expected(0, 2), expected(1, 1), expected(1, 2), expected(2, 2)};
      for (unsigned int c = 0; c < 6; ++c)
        maxDifference = std::max(maxDifference, std::abs(o[c] - values[c]));
    }
    CPPUNIT_ASSERT(maxDifference < 1e-9);
  }

  void CopyGeometry_UnknownInputType_Fails()
  {
    const auto dir = WorkDir("CopyGeometryFails");
    RunFails("CopyGeometry", {"-i", Data("sh_image_test.nii.gz"), "-r", Data("sh_image_test.nii.gz"), "-o", "out.nii.gz", "-t", "tensor"}, dir);
  }

  // -t sh on a tensor file: the tensor reader wins over the preferred SH reader, and the app
  // refuses the input instead of turning it as a tensor image.
  void CopyGeometry_InputTypeMismatch_Fails()
  {
    const auto dir = WorkDir("CopyGeometryMismatch");
    RunFails("CopyGeometry", {"-i", Data("tensorImage.dti"), "-r", Data("tensorImage.dti"), "-o", "out.dti", "-t", "sh"}, dir);
    CPPUNIT_ASSERT(!std::filesystem::exists(Path(dir, "out.dti")));
  }

  // No DICOM DWI in the test data: only the error path.
  void DiffusionDICOMLoader_MissingFolder_Fails()
  {
    const auto dir = WorkDir("DiffusionDICOMLoader");
    RunFails("DiffusionDICOMLoader", {"-i", Path(dir, "missing"), "-o", "dwi.nrrd"}, dir);
  }

  void ShToOdfImage_Odfs()
  {
    const auto dir = WorkDir("ShToOdfImage");
    RunOk("ShToOdfImage", {"-i", Data("sh_image_test.nii.gz"), "-o", "odfs.odf"}, dir);
    auto odfs = Load<mitk::OdfImage>(dir, "odfs.odf");
    CheckSize(odfs, 10, 10, 3);
  }

  void PeakExtraction_Peaks()
  {
    const auto dir = WorkDir("PeakExtraction");
    RunOk("PeakExtraction", {"-i", Data("sh_image_test.nii.gz"), "-o", "sub.01"}, dir);
    auto peaks = Load<mitk::PeakImage>(dir, "sub.01_PEAKS.nrrd", {"Peak Image"});
    CheckSize(peaks, 10, 10, 3);
    CPPUNIT_ASSERT_EQUAL(6u, peaks->GetDimension(3));
  }

  void DReg_SameImage_KeepsGrid()
  {
    const auto dir = WorkDir("DReg");
    const auto gfa = Data("StreamlineTractography/gfa_image.nii.gz");
    RunOk("DReg", {"-f", gfa, "-m", gfa, "-o", "registered.nii.gz"}, dir);
    CheckSize(Load<mitk::Image>(dir, "registered.nii.gz"), 64, 64, 3);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsMisc)
