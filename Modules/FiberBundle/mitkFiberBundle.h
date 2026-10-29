/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkFiberBundle_h
#define mitkFiberBundle_h

#include <MitkFiberBundleExports.h>
#include <mitkBaseData.h>
#include <mitkLookupTable.h>

#include <itkImage.h>
#include <itkScalableAffineTransform.h>
#include <vnl/vnl_vector_fixed.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

namespace mitk {

class DataNode;
class DataStorage;
class Image;

/**
 * \brief Tractogram: streamlines (fibers) as the lines of a vtkPolyData, with one weight per
 * fiber, one RGBA colour per point and the fiber length statistics.
 */
class MITKFIBERBUNDLE_EXPORT FiberBundle : public BaseData
{
public:

    using ItkUcharImgType = itk::Image<unsigned char, 3>;

    void UpdateOutputInformation() override;
    void SetRequestedRegionToLargestPossibleRegion() override;
    bool RequestedRegionIsOutsideOfTheBufferedRegion() override;
    bool VerifyRequestedRegion() override;
    void SetRequestedRegion(const itk::DataObject*) override;

    mitkClassMacro( FiberBundle, BaseData )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    mitkNewMacro1Param(Self, vtkSmartPointer<vtkPolyData>) // custom constructor

    // colorcoding related methods
    void ColorFibersByFiberWeights(bool opacity, LookupTable::LookupTableType type);
    void ColorFibersByCurvature(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type);
    void ColorFibersByLength(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type);
    void ColorFibersByScalarMap(const mitk::Image* image, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate=true);
    void ColorFibersByOrientation();
    void SetFiberColors(vtkSmartPointer<vtkUnsignedCharArray> fiberColors);
    void SetFiberColors(float r, float g, float b, float alpha=255);
    void SetSingleFiberColor(float r, float g, float b, unsigned int cellId, float alpha=255);
    vtkSmartPointer<vtkUnsignedCharArray> GetFiberColors() const { return m_FiberColors; }

    // fiber compression
    void Compress(float error = 0.0);

    // fiber resampling
    void ResampleSpline(float pointDistance=1);
    void ResampleSpline(float pointDistance, double tension, double continuity, double bias );
    void ResampleLinear(double pointDistance=1);
    /** Resamples every fiber to targetPoints points at equal arc length along it, both end
     *  points included (dipy's set_number_of_points), also fibers that have targetPoints points
     *  already. A fiber without length gets targetPoints copies of its position. */
    void ResampleToNumPoints(unsigned int targetPoints);

    /** The polyline through points resampled to numPoints points at equal arc length along it,
     *  both end points included; a polyline without length or of one point gives numPoints copies
     *  of its first point, an empty one nothing. Used by ResampleToNumPoints and the fiber
     *  clustering filter, so that both resample the same way. */
    static std::vector< vnl_vector_fixed<double, 3> > ResamplePolylineToNumPoints(const std::vector< vnl_vector_fixed<double, 3> >& points, unsigned int numPoints);

    mitk::FiberBundle::Pointer FilterByWeights(float weight_thr, bool invert=false) const;
    bool RemoveShortFibers(float lengthInMM);
    bool RemoveLongFibers(float lengthInMM);
    bool ApplyCurvatureThreshold(float minRadius, bool deleteFibers);
    void MirrorFibers(unsigned int axis);
    void RotateAroundAxis(double x, double y, double z);
    void TranslateFibers(double x, double y, double z);
    void ScaleFibers(double x, double y, double z, bool subtractCenter=true);
    void TransformFibers(double rx, double ry, double rz, double tx, double ty, double tz);
    void TransformFibers(const itk::ScalableAffineTransform< mitk::ScalarType >* transform);
    void RemoveDir(vnl_vector_fixed<double,3> dir, double threshold);

    template< class TType=float >
    void TransformPoint(itk::Point<TType, 3>& point, itk::Matrix< TType, 3, 3>& rot, TType& tx, TType& ty, TType& tz)
    {
      mitk::Point3D center = this->GetGeometry()->GetCenter();

      point[0] -= center[0];
      point[1] -= center[1];
      point[2] -= center[2];
      point = rot*point;
      point[0] += center[0]+tx;
      point[1] += center[1]+ty;
      point[2] += center[2]+tz;
    }

    // add/subtract fibers
    FiberBundle::Pointer AddBundle(const FiberBundle* fib) const;
    mitk::FiberBundle::Pointer AddBundles(const std::vector< mitk::FiberBundle::Pointer >& fibs) const;
    FiberBundle::Pointer SubtractBundle(const FiberBundle* fib) const;

    // fiber subset extraction
    FiberBundle::Pointer           ExtractFiberSubset(const DataNode* roi, const DataStorage* storage) const;
    std::vector<unsigned int>      ExtractFiberIdSubset(const DataNode* roi, const DataStorage* storage) const;
    FiberBundle::Pointer           RemoveFibersOutside(const ItkUcharImgType* mask, bool invert=false) const;
    mitk::FiberBundle::Pointer     SubsampleFibers(float factor, bool random_seed) const;

    // get/set data
    float GetFiberLength(unsigned int index) const { return m_FiberLengths.at(index); }
    vtkSmartPointer<vtkFloatArray> GetFiberWeights() const { return m_FiberWeights; }
    float GetFiberWeight(unsigned int fiber) const;
    void SetFiberWeights(float newWeight);
    void SetFiberWeight(unsigned int fiber, float weight);
    void SetFiberWeights(vtkSmartPointer<vtkFloatArray> weights);
    void SetFiberPolyData(vtkSmartPointer<vtkPolyData>);
    vtkSmartPointer<vtkPolyData> GetFiberPolyData() const;
    itkGetConstMacro( NumFibers, unsigned int)
    itkGetConstMacro( MinFiberLength, float )
    itkGetConstMacro( MaxFiberLength, float )
    itkGetConstMacro( MeanFiberLength, float )
    itkGetConstMacro( MedianFiberLength, float )
    itkGetConstMacro( LengthStDev, float )
    itkGetConstMacro( UpdateTime2D, itk::TimeStamp )
    itkGetConstMacro( UpdateTime3D, itk::TimeStamp )
    void RequestUpdate2D(){ m_UpdateTime2D.Modified(); }
    void RequestUpdate3D(){ m_UpdateTime3D.Modified(); }
    void RequestUpdate(){ m_UpdateTime2D.Modified(); m_UpdateTime3D.Modified(); }

    unsigned int GetNumberOfPoints() const;

    /** Independent copy of fibers, weights, colours, metadata and properties; the same as Clone(). */
    mitk::FiberBundle::Pointer GetDeepCopy() const;

    // compare fiber bundles
    bool Equals(const FiberBundle* fib, double eps=0.01) const;

    /** The fibers with the given ids, and their weights in weights, both in the order of fiberIds. */
    vtkSmartPointer<vtkPolyData>    GeneratePolyDataByIds(const std::vector<unsigned int>& fiberIds, vtkFloatArray* weights) const;

    // Structure to hold metadata of a TrackVis file
    struct TrackVis_header
    {
        char                id_string[6];
        short int           dim[3];
        float               voxel_size[3];
        float               origin[3];
        short int           n_scalars;
        char                scalar_name[10][20];
        short int           n_properties;
        char                property_name[10][20];
        float               vox_to_ras[4][4];
        char                reserved[444];
        char                voxel_order[4];
        char                pad2[4];
        float               image_orientation_patient[6];
        char                pad1[2];
        unsigned char       invert_x;
        unsigned char       invert_y;
        unsigned char       invert_z;
        unsigned char       swap_xy;
        unsigned char       swap_yz;
        unsigned char       swap_zx;
        int                 n_count;
        int                 version;
        int                 hdr_size;
    };

    /** The header set, read from a .trk file or set from a reference geometry; while there is
     *  none, a header describing the bundle's current bounding box. */
    TrackVis_header GetTrackVisHeader() const;
    /** Whether a header was set (read from a .trk file, set from a reference geometry or copied
     *  from a bundle that had one), as opposed to the one derived from the bounding box. */
    bool HasTrackVisHeader() const;
    void SetTrackVisHeader(const TrackVis_header &TrackVisHeader);
    void SetTrackVisHeader(const BaseGeometry *geometry);

    /** Fiber count, length, extent and weight statistics and the TrackVis reference grid, as
     *  PrintSelf prints them before the BaseData part (geometry objects and properties). */
    void PrintStatistics(std::ostream &os, itk::Indent indent) const;
    void PrintSelf(std::ostream &os, itk::Indent indent) const override;

    void setIsRAS(bool newIsRAS);

    bool IsRAS() const;

protected:

    FiberBundle( vtkPolyData* fiberPolyData = nullptr );
    FiberBundle(const FiberBundle& other);
    ~FiberBundle() override;

    mitkCloneMacro(Self)

    void                            UpdateFiberGeometry();

private:

    template <typename TPixel>
    void ColorFibersByScalarMap(const itk::Image<TPixel, 3>* image, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate);

    // actual fiber container
    vtkSmartPointer<vtkPolyData>  m_FiberPolyData;

    unsigned int m_NumFibers;

    vtkSmartPointer<vtkUnsignedCharArray> m_FiberColors;
    vtkSmartPointer<vtkFloatArray> m_FiberWeights;
    std::vector< float > m_FiberLengths;
    float   m_MinFiberLength;
    float   m_MaxFiberLength;
    float   m_MeanFiberLength;
    float   m_MedianFiberLength;
    float   m_LengthStDev;
    itk::TimeStamp m_UpdateTime2D;
    itk::TimeStamp m_UpdateTime3D;

    TrackVis_header     m_TrackVisHeader;

    bool m_IsRAS;
};

} // namespace mitk

#endif
