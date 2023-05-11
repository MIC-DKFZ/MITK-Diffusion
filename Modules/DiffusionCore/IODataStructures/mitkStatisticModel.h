/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef _MITK_StatisticModel_H
#define _MITK_StatisticModel_H

//includes for MITK datastructure
#include <mitkBaseData.h>
#include <MitkDiffusionCoreExports.h>
#include <mitkDataStorage.h>

//includes storing fiberdata
#include <vtkSmartPointer.h>

// OpenCV
#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>


namespace mitk {

/**
   * \brief Base Class for Fiber Bundles;   */
class MITKDIFFUSIONCORE_EXPORT StatisticModel : public BaseData
{
public:


    // mitkClassMacro( StatisticModel, BaseData )
    // itkFactorylessNewMacro(Self)
    // itkCloneMacro(Self)

    StatisticModel();
    virtual ~StatisticModel();

    // Member variable that represents the statistical model
    cv::Ptr<cv::ml::RTrees> m_StatModel;

protected:


private:

};

} // namespace mitk

#endif /*  _MITK_StatisticModel_H */
