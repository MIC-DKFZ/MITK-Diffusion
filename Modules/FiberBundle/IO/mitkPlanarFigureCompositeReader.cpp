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
#include <itkMetaDataObject.h>
#include <itksys/SystemTools.hxx>
#include <vtkCleanPolyData.h>
#include <mitkCustomMimeType.h>
#include "mitkFiberBundleMimeTypes.h"
#include <boost/property_tree/ptree.hpp>
#define RAPIDXML_NO_EXCEPTIONS
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <mitkLexicalCast.h>
#include <mitkPlanarFigureComposite.h>


mitk::PlanarFigureCompositeReader::PlanarFigureCompositeReader()
    : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE(), "Planar Figure Composite Reader" )
{
    m_ServiceReg = this->RegisterService();
}

mitk::PlanarFigureCompositeReader::PlanarFigureCompositeReader(const PlanarFigureCompositeReader &other)
    :mitk::AbstractFileReader(other)
{
}

mitk::PlanarFigureCompositeReader * mitk::PlanarFigureCompositeReader::Clone() const
{
    return new PlanarFigureCompositeReader(*this);
}


std::vector<itk::SmartPointer<mitk::BaseData> > mitk::PlanarFigureCompositeReader::DoRead()
{

    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    try
    {
        const std::string& locale = "C";
        const std::string& currLocale = setlocale( LC_ALL, nullptr );
        setlocale(LC_ALL, locale.c_str());

        std::string filename = this->GetInputLocation();

        std::string ext = itksys::SystemTools::GetFilenameLastExtension(filename);
        ext = itksys::SystemTools::LowerCase(ext);

        boost::property_tree::ptree tree;
        boost::property_tree::xml_parser::read_xml(filename, tree);
        int comptype = tree.get<int>("comptype");

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
            MITK_ERROR << filename << " contains no valid composition type!";
        }

        std::vector<itk::SmartPointer<mitk::BaseData> > result;
        result.push_back(pfc.GetPointer());

        setlocale(LC_ALL, currLocale.c_str());

        return result;
    }
    catch(...)
    {
        throw;
    }
    return result;
}
