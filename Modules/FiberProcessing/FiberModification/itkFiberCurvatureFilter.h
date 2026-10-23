/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkFiberCurvatureFilter_h
#define itkFiberCurvatureFilter_h

// MITK
#include <mitkFiberBundle.h>
#include <MitkFiberProcessingExports.h>

// ITK
#include <itkProcessObject.h>

namespace itk{

/**
* \brief Removes the strongly curved parts of fibers. At each point, the mean angle between the
* segment directions within +-Distance/2 mm along the fiber and their mean direction is the
* deviation; points whose deviation reaches AngularDeviation are removed and the fiber is split
* there (or, with RemoveFibers, the whole fiber is removed). Pieces keep their fiber's weight;
* single points are dropped.
*/
class MITKFIBERPROCESSING_EXPORT FiberCurvatureFilter : public ProcessObject
{
public:

    using Self = FiberCurvatureFilter;
    using Superclass = ProcessObject;
    using Pointer = SmartPointer< Self >;
    using ConstPointer = SmartPointer< const Self >;

    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    itkTypeMacro( FiberCurvatureFilter, ProcessObject )

    void Update() override{
        this->GenerateData();
    }

    /** Length of the window along the fiber in mm */
    itkSetMacro( Distance, double )
    /** Threshold of the mean angle in degrees */
    itkSetMacro( AngularDeviation, double )
    /** Remove whole fibers instead of the curved parts */
    itkSetMacro( RemoveFibers, bool )
    /** The reference direction is the mean of the two segments next to the point instead of the
     *  mean of all segments in the window (not a median, despite the name). */
    itkSetMacro( UseMedian, bool )
    itkSetConstObjectMacro( InputFiberBundle, mitk::FiberBundle )
    itkGetMacro( OutputFiberBundle, mitk::FiberBundle::Pointer )

protected:

    void GenerateData() override;

    FiberCurvatureFilter();
    ~FiberCurvatureFilter() override = default;

    mitk::FiberBundle::ConstPointer m_InputFiberBundle;
    mitk::FiberBundle::Pointer m_OutputFiberBundle;
    double  m_AngularDeviation;
    double  m_Distance;
    bool    m_RemoveFibers;
    bool    m_UseMedian;
};
}

#endif
