/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkPlanarFigureCompositeWriter_h
#define mitkPlanarFigureCompositeWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{

/** \brief Writes the operation (AND, OR, NOT) of a planar figure composite (.pfc). */
class PlanarFigureCompositeWriter : public mitk::AbstractFileWriter
{
public:

    PlanarFigureCompositeWriter();
    PlanarFigureCompositeWriter(const PlanarFigureCompositeWriter & other) = default;
    PlanarFigureCompositeWriter * Clone() const override;

    using mitk::AbstractFileWriter::Write;
    void Write() override;
};

}

#endif
