/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical Image Computing.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_DiffusionImageNrrdWriter__H_
#define _MITK_DiffusionImageNrrdWriter__H_

#include <mitkAbstractFileWriter.h>
#include <mitkDiffusionPropertyHelper.h>

namespace mitk
{

/**
 * Writes diffusion volumes to a file
 * @ingroup Process
 */
class DiffusionImageNrrdWriter : public mitk::AbstractFileWriter
{
public:

    DiffusionImageNrrdWriter();
    ~DiffusionImageNrrdWriter() override;

    using AbstractFileWriter::Write;
    void Write() override;

    ConfidenceLevel GetConfidenceLevel() const override;

    typedef mitk::DiffusionPropertyHelper::ImageType                        VectorImageType;
    typedef mitk::DiffusionPropertyHelper::GradientDirectionType            GradientDirectionType;
    typedef mitk::DiffusionPropertyHelper::MeasurementFrameType             MeasurementFrameType;
    typedef mitk::DiffusionPropertyHelper::GradientDirectionsContainerType  GradientDirectionContainerType;

  protected:

    DiffusionImageNrrdWriter(const DiffusionImageNrrdWriter& other);
    mitk::DiffusionImageNrrdWriter* Clone() const override;


};


} // end of namespace mitk

#endif
