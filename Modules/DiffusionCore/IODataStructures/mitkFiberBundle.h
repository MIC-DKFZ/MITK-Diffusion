/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef _MITK_FiberBundle_H
#define _MITK_FiberBundle_H

//includes for MITK datastructure
#include <mitkBaseData.h>
#include <MitkDiffusionCoreExports.h>
#include <mitkImage.h>
#include <mitkDataStorage.h>
#include <mitkPlanarFigure.h>
#include <mitkPixelTypeTraits.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkPeakImage.h>

//includes storing fiberdata
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkDataSet.h>
#include <vtkTransform.h>
#include <vtkFloatArray.h>
#include <itkScalableAffineTransform.h>
#include <mitkDiffusionFunctionCollection.h>
#include <mitkLookupTable.h>

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace mitk {

/**
   * \brief Base Class for Fiber Bundles;   */
class MITKDIFFUSIONCORE_EXPORT FiberBundle : public BaseData
{
public:

    typedef itk::Image<unsigned char, 3> ItkUcharImgType;

    // fiber colorcodings
    static const char* FIBER_ID_ARRAY;

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
    void ColorSinglePoint(int f_idx, int p_idx, double rgb[3]);
    void ColorFibersByFiberWeights(bool opacity, LookupTable::LookupTableType type);
    void ColorFibersByCurvature(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type);
    void ColorFibersByLength(bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type);
    void ColorFibersByScalarMap(mitk::Image::Pointer, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate=true);
    template <typename TPixel>
    void ColorFibersByScalarMap(typename itk::Image<TPixel, 3>::Pointer, bool opacity, bool weight_fibers, mitk::LookupTable::LookupTableType type, double max_cap, bool interpolate=true);
    void ColorFibersByOrientation();
    void SetFiberOpacity(vtkDoubleArray *FAValArray);
    void ResetFiberOpacity();
    void SetFiberColors(vtkSmartPointer<vtkUnsignedCharArray> fiberColors);
    void SetFiberColors(float r, float g, float b, float alpha=255);
    vtkSmartPointer<vtkUnsignedCharArray> GetFiberColors() const { return m_FiberColors; }

    // fiber compression
    void Compress(float error = 0.0);

    // fiber resampling
    void ResampleSpline(float pointDistance=1);
    void ResampleSpline(float pointDistance, double tension, double continuity, double bias );
    void ResampleLinear(double pointDistance=1);
    void ResampleToNumPoints(unsigned int targetPoints);

    mitk::FiberBundle::Pointer FilterByWeights(float weight_thr, bool invert=false);
    bool RemoveShortFibers(float lengthInMM);
    bool RemoveLongFibers(float lengthInMM);
    bool ApplyCurvatureThreshold(float minRadius, bool deleteFibers);
    void MirrorFibers(unsigned int axis);
    void RotateAroundAxis(double x, double y, double z);
    void TranslateFibers(double x, double y, double z);
    void ScaleFibers(double x, double y, double z, bool subtractCenter=true);
    void TransformFibers(double rx, double ry, double rz, double tx, double ty, double tz);
    void TransformFibers(itk::ScalableAffineTransform< mitk::ScalarType >::Pointer transform);
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

    template< class TType=float >
    void TransformPoint(itk::Point<TType, 3>& point, TType rx, TType ry, TType rz, TType tx, TType ty, TType tz)
    {
      auto rot = mitk::imv::GetRotationMatrixItk<TType>(rx, ry, rz);
      mitk::Point3D center = this->GetGeometry()->GetCenter();

      point[0] -= center[0];
      point[1] -= center[1];
      point[2] -= center[2];
      point = rot*point;
      point[0] += center[0]+tx;
      point[1] += center[1]+ty;
      point[2] += center[2]+tz;
    }

    itk::Matrix< double, 3, 3 > TransformMatrix(itk::Matrix< double, 3, 3 > m, double rx, double ry, double rz);

    // add/subtract fibers
    FiberBundle::Pointer AddBundle(FiberBundle* fib);
    mitk::FiberBundle::Pointer AddBundles(std::vector< mitk::FiberBundle::Pointer > fibs);
    FiberBundle::Pointer SubtractBundle(FiberBundle* fib);

    // fiber subset extraction
    FiberBundle::Pointer           ExtractFiberSubset(DataNode *roi, DataStorage* storage);
    std::vector<unsigned int>      ExtractFiberIdSubset(DataNode* roi, DataStorage* storage);
    FiberBundle::Pointer           RemoveFibersOutside(ItkUcharImgType* mask, bool invert=false);
    float                          GetOverlap(ItkUcharImgType* mask);
    std::tuple<float, float>       GetDirectionalOverlap(ItkUcharImgType* mask, mitk::PeakImage::ItkPeakImageType* peak_image);
    float                          GetNumEpFractionInMask(ItkUcharImgType* mask, bool different_label);
    mitk::FiberBundle::Pointer     SubsampleFibers(float factor, bool random_seed);

    // get/set data
    float GetFiberLength(unsigned int index) const { return m_FiberLengths.at(index); }
    vtkSmartPointer<vtkFloatArray> GetFiberWeights() const { return m_FiberWeights; }
    float GetFiberWeight(unsigned int fiber) const;
    void SetFiberWeights(float newWeight);
    void SetFiberWeight(unsigned int fiber, float weight);
    void SetFiberWeights(vtkSmartPointer<vtkFloatArray> weights);
    void SetFiberPolyData(vtkSmartPointer<vtkPolyData>, bool updateGeometry = true);
    vtkSmartPointer<vtkPolyData> GetFiberPolyData() const;
    itkGetConstMacro( NumFibers, unsigned int)
    //itkGetMacro( FiberSampling, int)
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

    // copy fiber bundle
    mitk::FiberBundle::Pointer GetDeepCopy();

    // compare fiber bundles
    bool Equals(FiberBundle* fib, double eps=0.01);

    vtkSmartPointer<vtkPolyData>    GeneratePolyDataByIds(std::vector<unsigned int> fiberIds, vtkSmartPointer<vtkFloatArray> weights);

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

    TrackVis_header GetTrackVisHeader();
    void SetTrackVisHeader(const TrackVis_header &TrackVisHeader);
    void SetTrackVisHeader(BaseGeometry *geometry);

    void PrintSelf(std::ostream &os, itk::Indent indent) const override;

protected:

    FiberBundle( vtkPolyData* fiberPolyData = nullptr );
    ~FiberBundle() override;

    void                            GenerateFiberIds();
    void                            UpdateFiberGeometry();

private:

    // actual fiber container
    vtkSmartPointer<vtkPolyData>  m_FiberPolyData;

    // contains fiber ids
    vtkSmartPointer<vtkDataSet>   m_FiberIdDataSet;

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
};

} // namespace mitk

#endif /*  _MITK_FiberBundle_H */
