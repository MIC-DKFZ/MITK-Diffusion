/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkPlanarFigureCompositeWriter.h"
#include "mitkFiberBundleMimeTypes.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <mitkPlanarFigureComposite.h>

mitk::PlanarFigureCompositeWriter::PlanarFigureCompositeWriter()
    : mitk::AbstractFileWriter(mitk::PlanarFigureComposite::GetStaticNameOfClass(), mitk::FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE_NAME(), "Planar Figure Composite Writer")
{
    RegisterService();
}

mitk::PlanarFigureCompositeWriter * mitk::PlanarFigureCompositeWriter::Clone() const
{
    return new mitk::PlanarFigureCompositeWriter(*this);
}

void mitk::PlanarFigureCompositeWriter::Write()
{
    mitk::PlanarFigureComposite::ConstPointer input = dynamic_cast<const mitk::PlanarFigureComposite*>(this->GetInput());
    if (input.IsNull())
      mitkThrow() << "The planar figure composite writer needs a planar figure composite as input.";

    // The XML is written to a file by name, so an output stream gets a local file that is copied
    // to the stream when it goes out of scope. No locale switch: the file holds one integer.
    LocalFile localFile(this);

    boost::property_tree::ptree type;
    type.put("comptype", input->getOperationType());
    boost::property_tree::xml_writer_settings<std::string> writerSettings(' ', 2);
    boost::property_tree::xml_parser::write_xml(localFile.GetFileName(), type, std::locale(), writerSettings);
}
