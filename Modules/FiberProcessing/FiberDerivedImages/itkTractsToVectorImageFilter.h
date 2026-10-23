#ifndef itkTractsToVectorImageFilter_h
#define itkTractsToVectorImageFilter_h

// MITK
#include <mitkFiberBundle.h>

// ITK
#include <itkImageSource.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>

#include <vnl/vnl_vector_fixed.h>

#include <vector>

namespace itk{

/**
* \brief Voxel-wise main fiber directions of a tractogram (peak image): the segment directions in a
* voxel, weighted by in-voxel length and fiber weight, are clustered by mean shift within the
* angular threshold; each cluster is a peak whose size is the summed weighted length of its members.
* The output has MaxNumDirections x 3 components along its 4th axis.
*/
template< class PixelType >
class TractsToVectorImageFilter : public ImageSource< Image< PixelType, 4 > >
{

public:

  enum class NormalizationMethods {
    /** all peaks divided by the largest peak of the image */
    GLOBAL_MAX,
    /** each peak normalized to length 1 */
    SINGLE_VEC_NORM,
    /** all peaks of a voxel divided by the voxel's largest peak (0-1) */
    MAX_VEC_NORM
  };

  using Self = TractsToVectorImageFilter;
  using Superclass = ImageSource< Image< PixelType, 4 > >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using OutputVectorType = itk::Vector<float,3>;
  using OutputImageType = itk::Image<OutputVectorType, 3>;

  using DirectionType = vnl_vector_fixed< double, 3 >;
  using DirectionContainerType = VectorContainer< unsigned int, DirectionType >;
  using ContainerType = VectorContainer< unsigned int, DirectionContainerType::Pointer >;
  using ItkDirectionImageType = Image< PixelType, 4 >;
  using ItkUcharImgType = itk::Image<unsigned char, 3>;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractsToVectorImageFilter, ImageSource )

  /** Peaks smaller than this fraction of the largest peak of the voxel are discarded. */
  itkSetMacro( SizeThreshold, float)
  itkGetMacro( SizeThreshold, float)
  /** Cosine of the angle within which fiber directions are clustered into one peak */
  itkSetMacro( AngularThreshold, float)
  itkGetMacro( AngularThreshold, float)
  /** Normalization method of the peaks. Written out because itkSetMacro streams the value,
   *  which an enum class does not support. */
  void SetNormalizationMethod(NormalizationMethods method)
  {
    if (m_NormalizationMethod != method)
    {
      m_NormalizationMethod = method;
      this->Modified();
    }
  }
  /** Only the largest peaks of a voxel are kept; at least 1. */
  itkSetMacro( MaxNumDirections, unsigned long)
  itkGetMacro( MaxNumDirections, unsigned long)
  /** Only voxels inside the mask are processed; without it, the grid covers the bundle. */
  itkSetMacro( MaskImage, ItkUcharImgType::Pointer)
  itkSetConstObjectMacro( FiberBundle, mitk::FiberBundle )
  /** Use the mask only for its grid, not its content. */
  itkSetMacro( OnlyUseMaskGeometry, bool)
  /** Number of peaks per voxel */
  itkGetMacro( NumDirectionsImage, ItkUcharImgType::Pointer)
  /** The peaks, MaxNumDirections x 3 components per voxel */
  itkGetMacro( DirectionImage, typename ItkDirectionImageType::Pointer)

  void GenerateData() override;

protected:

  /** Mean-shift clustering of the directions inDirs with their lengths, repeated on the cluster
   *  means until no two clusters merge any more. */
  DirectionContainerType::Pointer FastClustering(DirectionContainerType::Pointer inDirs, std::vector< double > lengths);

  TractsToVectorImageFilter();
  ~TractsToVectorImageFilter() override = default;

  NormalizationMethods                m_NormalizationMethod;
  mitk::FiberBundle::ConstPointer     m_FiberBundle;
  float                               m_AngularThreshold;
  ItkUcharImgType::Pointer            m_MaskImage;
  unsigned long                       m_MaxNumDirections;
  float                               m_SizeThreshold;
  bool                                m_OnlyUseMaskGeometry;

  // outputs
  typename ItkDirectionImageType::Pointer m_DirectionImage;
  ItkUcharImgType::Pointer                m_NumDirectionsImage;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractsToVectorImageFilter.cpp"
#endif

#endif
