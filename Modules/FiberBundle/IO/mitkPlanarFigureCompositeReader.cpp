/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkPlanarFigureCompositeReader.h"
#include "mitkFiberBundleMimeTypes.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <mitkPlanarFigureComposite.h>


mitk::PlanarFigureCompositeReader::PlanarFigureCompositeReader()
    : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE_NAME(), "Planar Figure Composite Reader" )
{
    this->RegisterService();
}

mitk::PlanarFigureCompositeReader * mitk::PlanarFigureCompositeReader::Clone() const
{
    return new PlanarFigureCompositeReader(*this);
}


std::vector<itk::SmartPointer<mitk::BaseData> > mitk::PlanarFigureCompositeReader::DoRead()
{
    // No locale switch: the file holds one integer, and a switch by hand stayed in force for the
    // whole process when the XML parser threw.
    const std::string filename = this->GetInputLocation();

    boost::property_tree::ptree tree;
    boost::property_tree::xml_parser::read_xml(filename, tree);
    const int comptype = tree.get<int>("comptype");

    mitk::PlanarFigureComposite::Pointer pfc = mitk::PlanarFigureComposite::New();
    switch(comptype)
    {
    case 0:
        pfc->setOperationType(mitk::PlanarFigureComposite::AND);
        MITK_INFO << "loading AND composition";
        break;
    case 1:
        pfc->setOperationType(mitk::PlanarFigureComposite::OR);
        MITK_INFO << "loading OR composition";
        break;
    case 2:
        pfc->setOperationType(mitk::PlanarFigureComposite::NOT);
        MITK_INFO << "loading NOT composition";
        break;
    default:
        // A composite with a default operation would silently extract other fibers.
        mitkThrow() << filename << " contains no valid composition type (" << comptype << ").";
    }

    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    result.push_back(pfc.GetPointer());
    return result;
}
