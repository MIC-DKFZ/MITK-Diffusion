/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionTestData.h>
#include <mitkFileReaderSelector.h>
#include <mitkIOMimeTypes.h>
#include <mitkIOUtil.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>

/**
 * Opens one file of every diffusion data type the way the Workbench does: through mitk::IOUtil,
 * in a process that links MITK's core but no diffusion module. The readers are registered by the
 * service activators of DiffusionImage, DiffusionModelling and FiberBundle, which run when the
 * module is loaded. The Workbench loads a diffusion module only when a diffusion plugin starts, so
 * the readers are there right after its start only because MitkCore auto-loads the modules'
 * autoload/IO modules. The test driver of a diffusion module links the module and would find the
 * readers without that, hence this driver of its own.
 *
 * Each case checks that the reader is offered for the file once and that the file, opened with
 * it as a user chooses it in the open dialog, gives the expected data type. Formats the test
 * data has no file of are written first by the diffusion writers from one it has.
 */
class mitkDiffusionWorkbenchReadersTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionWorkbenchReadersTestSuite);
  MITK_TEST(Dwi_NrrdNiftiDicom_ReadAsDwi);
  MITK_TEST(Tensor_Dti_ReadAsTensorImage);
  MITK_TEST(Odf_OdfAndQbi_ReadAsOdfImage);
  MITK_TEST(Sh_NrrdAndNiftiLegacyAndMrtrix_ReadAsShImage);
  MITK_TEST(Peaks_NrrdAndNifti_ReadAsPeakImage);
  MITK_TEST(Tractograms_FibVtkTrkTckTrxDicom_ReadAsFiberBundle);
  MITK_TEST(PlanarFigureComposite_Pfc_ReadAsComposite);
  CPPUNIT_TEST_SUITE_END();

private:
  std::string m_TempDir;

  /** The reader descriptions offered for the file; a reader registered twice is there twice. */
  static std::multiset<std::string> OfferedReaders(const std::string &path)
  {
    std::multiset<std::string> descriptions;
    for (const auto &item : mitk::FileReaderSelector(path).Get())
      descriptions.insert(item.GetDescription());
    return descriptions;
  }

  /** Opens the file with the named reader and checks that it gives one object of the class. */
  static mitk::BaseData::Pointer Load(const std::string &path, const std::string &reader, const std::string &className)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Reader \"" + reader + "\" offered for " + path,
                                 std::size_t(1),
                                 OfferedReaders(path).count(reader));

    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({reader}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    const auto data = mitk::IOUtil::Load(path, &functor);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(path, std::size_t(1), data.size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE(path, className, std::string(data.front()->GetNameOfClass()));
    return data.front();
  }

  /** A DWI is a plain mitk::Image; the DWI readers give it the gradients as a property. */
  static void LoadDwi(const std::string &path)
  {
    auto data = Load(path, "Diffusion Weighted Images", "Image");
    CPPUNIT_ASSERT_MESSAGE(path, data->GetProperty("DWMRI.GradientDirections").IsNotNull());
  }

  /** Copies a test-data file into the temporary folder under another name (extension). */
  std::string Copy(const std::string &file, const std::string &name) const
  {
    const std::string path = m_TempDir + "/" + name;
    std::filesystem::copy_file(mitk::GetDiffusionTestDataFilePath(file), path);
    return path;
  }

  static std::string FileContent(const std::string &path)
  {
    std::ifstream file(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
  }

public:
  void setUp() override
  {
    m_TempDir = mitk::IOUtil::CreateTemporaryDirectory("DiffusionWorkbenchReadersTest_XXXXXX");
  }

  void tearDown() override
  {
    std::filesystem::remove_all(m_TempDir);
  }

  void Dwi_NrrdNiftiDicom_ReadAsDwi()
  {
    LoadDwi(mitk::GetDiffusionTestDataFilePath("dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"));
    LoadDwi(mitk::GetDiffusionTestDataFilePath("diffusionImage.hdwi"));
    LoadDwi(Copy("dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi", "dwi.nrrd"));
    LoadDwi(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/diffusion_image.nii.gz"));
    // The DICOM reader reads the series of the folder.
    LoadDwi(mitk::GetDiffusionTestDataFilePath("DicomDwi/Siemens/v0_z0.dcm"));
  }

  void Tensor_Dti_ReadAsTensorImage()
  {
    Load(mitk::GetDiffusionTestDataFilePath("tensorImage.dti"), "Diffusion Tensor Image", "TensorImage");
  }

  void Odf_OdfAndQbi_ReadAsOdfImage()
  {
    Load(mitk::GetDiffusionTestDataFilePath("qBallImage.qbi"), "ODF Image", "OdfImage");
    Load(Copy("qBallImage.qbi", "odf.odf"), "ODF Image", "OdfImage");
  }

  void Sh_NrrdAndNiftiLegacyAndMrtrix_ReadAsShImage()
  {
    // Legacy basis: 4D NIfTI without mark, and the coefficients as vector pixels of a 3D NRRD.
    auto sh = Load(mitk::GetDiffusionTestDataFilePath("sh_image_test.nii.gz"), "SH Image", "ShImage");
    Load(mitk::GetDiffusionTestDataFilePath("qBallImage_SHCoeffs.nrrd"), "SH Image", "ShImage");
    // MRtrix 3 basis: NIfTI written by MRtrix, and NRRD written by MITK with the mark.
    Load(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_image_sh.nii.gz"), "SH Image", "ShImage");
    const std::string nrrd = m_TempDir + "/sh.nrrd";
    mitk::IOUtil::Save(sh, "SH_IMAGE", nrrd);
    CPPUNIT_ASSERT(FileContent(nrrd).find("MITK_SH_convention:=MRtrix3") != std::string::npos);
    Load(nrrd, "SH Image", "ShImage");
  }

  void Peaks_NrrdAndNifti_ReadAsPeakImage()
  {
    auto peaks = Load(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), "Peak Image", "PeakImage");
    Load(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/csd_peak_image.nii.gz"), "Peak Image", "PeakImage");
    const std::string nrrd = m_TempDir + "/peaks.nrrd";
    mitk::IOUtil::Save(peaks, "ODF_PEAKS", nrrd);
    Load(nrrd, "Peak Image", "PeakImage");
  }

  void Tractograms_FibVtkTrkTckTrxDicom_ReadAsFiberBundle()
  {
    const std::string vtkReader = "VTK Fiber Bundle Reader";
    auto fibers = Load(mitk::GetDiffusionTestDataFilePath("fiberBundleX.fib"), vtkReader, "FiberBundle");
    Load(Copy("fiberBundleX.fib", "fibers.vtk"), vtkReader, "FiberBundle");
    Load(mitk::GetDiffusionTestDataFilePath("ThirdPartyTractograms/dipy_las.trk"),
         "TrackVis Fiber Bundle Reader",
         "FiberBundle");
    Load(mitk::GetDiffusionTestDataFilePath("ThirdPartyTractograms/mrtrix.tck"),
         "tck Fiber Bundle Reader (MRtrix format)",
         "FiberBundle");
    Load(mitk::GetDiffusionTestDataFilePath("ThirdPartyTractograms/TRX/gs.trx"), "TRX Fiber Bundle Reader", "FiberBundle");
    const std::string dicom = m_TempDir + "/tractogram.dcm";
    mitk::IOUtil::Save(fibers, mitk::IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.dcm", dicom);
    Load(dicom, "DICOM Fiber Bundle Reader", "FiberBundle");
  }

  void PlanarFigureComposite_Pfc_ReadAsComposite()
  {
    const std::string path = m_TempDir + "/composite.pfc";
    std::ofstream(path) << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<comptype>1</comptype>\n";
    Load(path, "Planar Figure Composite Reader", "PlanarFigureComposite");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionWorkbenchReaders)
