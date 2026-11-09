/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkShImageReader.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"

#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "itkNiftiImageIO.h"
#include "mitkITKImageImport.h"
#include "mitkImageDataItem.h"
#include <mitkLocaleSwitch.h>
#include <itkShCoefficientImageImporter.h>

#include <mitkDiffusionImageHelperFunctions.h>
#include <mitkDiffusionModellingHelperFunctions.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <vnl/algo/vnl_svd.h>

#include <cmath>
#include <optional>

namespace
{
  // Constant initialization: the reader is constructed while the library is loaded, possibly
  // before dynamically initialized globals of this file.
  constexpr const char *CONVENTION_OPTION = "SH convention";
  constexpr const char *AUTO = "auto (file marker, else MITK legacy; vector layout: MITK FSL)";
  constexpr const char *MRTRIX3 = "MRtrix 3";
  constexpr const char *MITK_LEGACY = "MITK legacy (MRtrix 0.2)";
  constexpr const char *FSL_DIPY = "FSL/Dipy";

  /** How the coefficients in the file are defined. */
  enum class FileBasis
  {
    MitkLegacy, // MITK's internal basis (MRtrix 0.2), taken over unchanged
    Mrtrix3,    // MRtrix 3: orthonormal, m != 0 larger by sqrt(2) than MITK's
    MitkFsl,    // MITK's FSL definition (ShImage::FSL), taken over unchanged
    Dipy        // dipy's default basis (descoteaux07, legacy) in the FSL bvec frame
  };

  /**
   * Matrix that maps coefficients of dipy's default basis to MITK's basis of the same order.
   * dipy's basis (descoteaux07 with legacy=True, the default of its CSA, CSD and sf_to_sh) is
   * sqrt(2) Re Y_l^|m| for m < 0, Y_l^0 and sqrt(2) Im Y_l^|m| for m > 0, with the Condon-Shortley
   * phase, at directions in the frame of the bvecs dipy was given: the FSL frame, i.e. voxel axes
   * with x negated when the image matrix has a positive determinant. MITK's coefficients are
   * functions of world (LPS) directions. Both bases span the same space, so sampling dipy's
   * function on a dense direction set and fitting MITK's basis is exact, also for oblique images.
   */
  vnl_matrix<double> DipyToMitkBasis(int order, const itk::Matrix<double, 3, 3> &direction)
  {
    // bvec-frame direction v = frame * u for a world direction u (frame axes as rows).
    const vnl_matrix_fixed<double, 3, 3> frame = mitk::gradients::FslBvecsMeasurementFrame(direction);

    // Fibonacci lattice: unit vectors spread evenly over the sphere, far more than coefficients.
    const unsigned int count = 1000;
    vnl_matrix<double> world(3, count);
    const double golden = itk::Math::pi * (3.0 - std::sqrt(5.0));
    for (unsigned int i = 0; i < count; ++i)
    {
      const double z = 1.0 - (2.0 * i + 1.0) / count;
      const double r = std::sqrt(1.0 - z * z);
      world(0, i) = r * std::cos(golden * i);
      world(1, i) = r * std::sin(golden * i);
      world(2, i) = z;
    }

    const unsigned int coefficients = (order + 1) * (order + 2) / 2;
    vnl_matrix<double> dipy(count, coefficients);
    for (unsigned int i = 0; i < count; ++i)
    {
      const vnl_vector_fixed<double, 3> v = frame * world.get_column(i);
      const double theta = std::acos(std::max(-1.0, std::min(1.0, v[2])));
      const double phi = std::atan2(v[1], v[0]);
      for (int l = 0; l <= order; l += 2)
        for (int m = -l; m <= l; ++m)
        {
          const auto y = boost::math::spherical_harmonic<double>(static_cast<unsigned int>(l), std::abs(m), theta, phi);
          dipy(i, l * (l + 1) / 2 + m) = m < 0 ? std::sqrt(2.0) * y.real() : (m > 0 ? std::sqrt(2.0) * y.imag() : y.real());
        }
    }

    const vnl_matrix<float> mitkBasis = mitk::sh::CalcShBasisForDirections(order, world, true);
    vnl_matrix<double> mitk(mitkBasis.rows(), mitkBasis.cols());
    for (unsigned int r = 0; r < mitk.rows(); ++r)
      for (unsigned int c = 0; c < mitk.cols(); ++c)
        mitk(r, c) = mitkBasis(r, c);
    return vnl_svd<double>(mitk).solve(dipy);
  }

  /** Convention marked by ShImageWriter (NRRD key or NIfTI descrip); it overrides the option. */
  std::optional<FileBasis> MitkMarker(const itk::MetaDataDictionary &dictionary)
  {
    std::string value;
    if (itk::ExposeMetaData<std::string>(dictionary, "MITK_SH_convention", value)
        || (itk::ExposeMetaData<std::string>(dictionary, "ITK_FileNotes", value) && value.rfind("MITK_SH_convention=", 0) == 0))
    {
      value = value.substr(value.find('=') + 1);
      if (value == "FSL")
        return FileBasis::MitkFsl;
      if (value == "MRtrix3")
        return FileBasis::Mrtrix3;
    }
    return std::nullopt;
  }

  /** MRtrix writes "MRtrix version: <version>" into the NIfTI descrip; version 3 and later use the orthonormal basis. */
  bool WrittenByMrtrix3(const itk::MetaDataDictionary &dictionary)
  {
    const std::string prefix = "MRtrix version: ";
    std::string notes;
    if (!itk::ExposeMetaData<std::string>(dictionary, "ITK_FileNotes", notes) || notes.rfind(prefix, 0) != 0)
      return false;
    try
    {
      return std::stoi(notes.substr(prefix.size())) >= 3;
    }
    catch (const std::exception &)
    {
      return false;
    }
  }
}

namespace mitk
{
  ShImageReader::ShImageReader(const ShImageReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  ShImageReader::ShImageReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::SH_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::SH_MIMETYPE_DESCRIPTION() )
  {
    Options defaultOptions;
    defaultOptions[CONVENTION_OPTION] = std::string(AUTO);
    defaultOptions[std::string(CONVENTION_OPTION) + ".enum"] = std::vector<std::string>{AUTO, MRTRIX3, MITK_LEGACY, FSL_DIPY};
    this->SetDefaultOptions(defaultOptions);

    this->RegisterService();
  }

  ShImageReader::~ShImageReader() = default;

  IFileReader::ConfidenceLevel ShImageReader::GetConfidenceLevel() const
  {
    // The vector length alone does not make a 3D vector image an SH image, so MITK's image reader
    // stays the default for the vector layout; the SH reader reads it when chosen.
    const auto level = AbstractFileReader::GetConfidenceLevel();
    if (level == Supported && DiffusionIOMimeTypes::IsShVectorLayout(this->GetInputLocation()))
      return PartiallySupported;
    return level;
  }

  template <int shOrder>
  mitk::Image::Pointer ShImageReader::ConvertShImage(ShImage::ShOnDiskType::Pointer img)
  {
    typedef itk::ShCoefficientImageImporter< float, shOrder > ImporterType;
    typename ImporterType::Pointer importer = ImporterType::New();
    importer->SetInputImage(img);
    importer->GenerateData();

    itk::Matrix<double, 3, 3> direction;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        direction[r][c] = img->GetDirection()[r][c];
    return TakeCoefficients<shOrder>(importer->GetCoefficientImage(), img->GetMetaDataDictionary(), direction, false);
  }

  template <int shOrder>
  mitk::Image::Pointer ShImageReader::ConvertVectorShImage(VectorImageType::Pointer img)
  {
    auto coefficients = CoefficientImageType<shOrder>::New();
    coefficients->SetSpacing(img->GetSpacing());
    coefficients->SetOrigin(img->GetOrigin());
    coefficients->SetDirection(img->GetDirection());
    coefficients->SetRegions(img->GetLargestPossibleRegion());
    coefficients->Allocate();
    itk::ImageRegionConstIterator<VectorImageType> itVector(img, img->GetLargestPossibleRegion());
    itk::ImageRegionIterator<CoefficientImageType<shOrder>> it(coefficients, coefficients->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it, ++itVector)
    {
      typename CoefficientImageType<shOrder>::PixelType pixel;
      for (unsigned int j = 0; j < pixel.Size(); ++j)
        pixel[j] = itVector.Get()[j];
      it.Set(pixel);
    }
    return TakeCoefficients<shOrder>(coefficients, img->GetMetaDataDictionary(), img->GetDirection(), true);
  }

  template <int shOrder>
  mitk::Image::Pointer ShImageReader::TakeCoefficients(CoefficientImageType<shOrder> *coefficients,
                                                       const itk::MetaDataDictionary &dictionary,
                                                       const itk::Matrix<double, 3, 3> &direction,
                                                       bool vectorLayout)
  {
    mitk::ShImage::Pointer shImage = mitk::ShImage::New();

    auto basis = MitkMarker(dictionary);
    if (!basis)
    {
      const std::string option = us::any_cast<std::string>(this->GetOptions()[CONVENTION_OPTION]);
      if (option == MRTRIX3 || (option == AUTO && WrittenByMrtrix3(dictionary)))
        basis = FileBasis::Mrtrix3;
      else if (option == FSL_DIPY)
        basis = FileBasis::Dipy;
      else if (option == AUTO && vectorLayout)
        // The vector layout comes from MITK versions whose Q-ball filters used MITK's FSL definition.
        basis = FileBasis::MitkFsl;
      else
        basis = FileBasis::MitkLegacy;
    }

    if (*basis == FileBasis::MitkFsl)
      shImage->SetShConvention(mitk::ShImage::SH_CONVENTION::FSL);

    if (*basis == FileBasis::Dipy)
    {
      const vnl_matrix<double> toMitk = DipyToMitkBasis(shOrder, direction);
      itk::ImageRegionIterator<CoefficientImageType<shOrder>> it(coefficients, coefficients->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
      {
        auto pixel = it.Get();
        vnl_vector<double> dipy(pixel.Size());
        for (unsigned int j = 0; j < pixel.Size(); ++j)
          dipy[j] = pixel[j];
        const vnl_vector<double> mitk = toMitk * dipy;
        for (unsigned int j = 0; j < pixel.Size(); ++j)
          pixel[j] = static_cast<float>(mitk[j]);
        it.Set(pixel);
      }
    }

    if (*basis == FileBasis::Mrtrix3)
    {
      // MITK's internal basis is MRtrix 0.2's, whose functions with m != 0 are MRtrix 3's divided
      // by sqrt(2); m = 0 sits at index l(l+1)/2 of each even order l.
      itk::ImageRegionIterator<CoefficientImageType<shOrder>> it(coefficients, coefficients->GetLargestPossibleRegion());
      for (; !it.IsAtEnd(); ++it)
      {
        auto pixel = it.Get();
        for (int l = 2; l <= shOrder; l += 2)
          for (int m = -l; m <= l; ++m)
            if (m != 0)
              pixel[l * (l + 1) / 2 + m] = static_cast<float>(pixel[l * (l + 1) / 2 + m] * std::sqrt(2.0));
        it.Set(pixel);
      }
    }

    // The SH image takes over the buffer of the coefficient image.
    mitk::Image::Pointer resultImage = dynamic_cast<mitk::Image*>(shImage.GetPointer());
    mitk::GrabItkImageMemoryChannel(coefficients, nullptr, resultImage.GetPointer(), false);
    return resultImage;
  }

  mitk::Image::Pointer ShImageReader::ReadVectorLayout(const std::string &location)
  {
    auto reader = itk::ImageFileReader<VectorImageType>::New();
    reader->SetFileName(location);
    reader->Update();
    VectorImageType::Pointer img = reader->GetOutput();
    switch (img->GetNumberOfComponentsPerPixel())
    {
    case 6:
      return ConvertVectorShImage<2>(img);
    case 15:
      return ConvertVectorShImage<4>(img);
    case 28:
      return ConvertVectorShImage<6>(img);
    case 45:
      return ConvertVectorShImage<8>(img);
    case 66:
      return ConvertVectorShImage<10>(img);
    case 91:
      return ConvertVectorShImage<12>(img);
    default:
      mitkThrow() << "SH order larger 12 not supported";
    }
  }

  std::vector<itk::SmartPointer<BaseData> > ShImageReader::DoRead()
  {
    mitk::LocaleSwitch localeSwitch("C");
    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    std::string location = GetInputLocation();

    if ( location == "")
    {
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename is empty!");
    }
    else
    {
      try
      {
        if (DiffusionIOMimeTypes::IsShVectorLayout(location))
        {
          result.push_back(ReadVectorLayout(location).GetPointer());
          return result;
        }

        std::string ext = itksys::SystemTools::GetFilenameLastExtension(location);
        typedef itk::ImageFileReader< ShImage::ShOnDiskType > FileReaderType;
        FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetFileName(location);
        if (ext==".shi")
        {
          itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
          reader->SetImageIO(io);
        }
        reader->Update();
        ShImage::ShOnDiskType::Pointer img = reader->GetOutput();

        switch (img->GetLargestPossibleRegion().GetSize()[3])
        {
        case 6:
          result.push_back( ConvertShImage<2>(img).GetPointer() );
          break;
        case 15:
          result.push_back( ConvertShImage<4>(img).GetPointer() );
          break;
        case 28:
          result.push_back( ConvertShImage<6>(img).GetPointer() );
          break;
        case 45:
          result.push_back( ConvertShImage<8>(img).GetPointer() );
          break;
        case 66:
          result.push_back( ConvertShImage<10>(img).GetPointer() );
          break;
        case 91:
          result.push_back( ConvertShImage<12>(img).GetPointer() );
          break;
        default :
          mitkThrow() << "SH order larger 12 not supported";
        }
      }
      catch(std::exception& e)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
      }
      catch(...)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested SH file!");
      }
    }
    return result;
  }

} //namespace MITK

mitk::ShImageReader* mitk::ShImageReader::Clone() const
{
  return new ShImageReader(*this);
}
