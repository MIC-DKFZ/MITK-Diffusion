/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkPlanarFigureCompositeReader_h
#define mitkPlanarFigureCompositeReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /** \brief Reads the operation (AND, OR, NOT) of a planar figure composite (.pfc). */
  class PlanarFigureCompositeReader : public AbstractFileReader
  {
  public:

    PlanarFigureCompositeReader();
    PlanarFigureCompositeReader(const PlanarFigureCompositeReader& other) = default;
    PlanarFigureCompositeReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
