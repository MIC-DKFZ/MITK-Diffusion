/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkConnectomicsNetworkReader.h"
#include "mitkConnectomicsNetworkDefinitions.h"
#include <tinyxml2.h>
#include "itksys/SystemTools.hxx"
#include <vtkMatrix4x4.h>
#include "mitkGeometry3D.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"
#include <tinyxml2.h>

namespace mitk
{

  ConnectomicsNetworkReader::ConnectomicsNetworkReader(const ConnectomicsNetworkReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  ConnectomicsNetworkReader::ConnectomicsNetworkReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::CONNECTOMICS_MIMETYPE() ), mitk::DiffusionIOMimeTypes::CONNECTOMICS_MIMETYPE_DESCRIPTION() )
  {
    m_ServiceReg = this->RegisterService();
  }

  ConnectomicsNetworkReader::~ConnectomicsNetworkReader()
  {
  }

  std::vector<itk::SmartPointer<BaseData> > ConnectomicsNetworkReader::DoRead()
  {
    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    std::string location = GetInputLocation();

    std::string ext = itksys::SystemTools::GetFilenameLastExtension(location);
    ext = itksys::SystemTools::LowerCase(ext);

    if ( location == "")
    {
      MITK_ERROR << "No file name specified.";
    }
    else if (ext == ".cnf")
    {
      try
      {
        mitk::ConnectomicsNetwork::Pointer outputNetwork = mitk::ConnectomicsNetwork::New();

        tinyxml2::XMLDocument doc;
        bool loadOkay = doc.LoadFile(location.c_str());
        if(!loadOkay)
        {
          mitkThrow() << "Could not open file " << location << " for reading.";
        }

        tinyxml2::XMLHandle hDoc(&doc);
        tinyxml2::XMLElement* pElem;
        tinyxml2::XMLHandle hRoot = hDoc.FirstChildElement();

        pElem = hRoot.ToElement();

        // save this for later
        hRoot = tinyxml2::XMLHandle(pElem);

        //get file version
        std::string version = std::string(pElem->Attribute(mitk::ConnectomicsNetworkDefinitions::XML_FILE_VERSION));

        // read geometry
        pElem = hRoot.FirstChildElement(mitk::ConnectomicsNetworkDefinitions::XML_GEOMETRY).ToElement();

        mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();

        // read origin
        mitk::Point3D origin;
        origin[0] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_X);
        origin[1] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_Y);
        origin[2] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_Z);
        geometry->SetOrigin(origin);

        // read spacing
        ScalarType spacing[3];
        spacing[0] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_X);
        spacing[1] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_Y);
        spacing[2] = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_Z);
        geometry->SetSpacing(spacing);

        // read transform
		double temp = 0;
        vtkMatrix4x4* m = vtkMatrix4x4::New();
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XX);
        m->SetElement(0,0,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XY);
        m->SetElement(1,0,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XZ);
        m->SetElement(2,0,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YX);
        m->SetElement(0,1,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YY);
        m->SetElement(1,1,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YZ);
        m->SetElement(2,1,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZX);
        m->SetElement(0,2,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZY);
        m->SetElement(1,2,temp);
        temp = pElem->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZZ);
        m->SetElement(2,2,temp);

        m->SetElement(0,3,origin[0]);
        m->SetElement(1,3,origin[1]);
        m->SetElement(2,3,origin[2]);
        m->SetElement(3,3,1);
        geometry->SetIndexToWorldTransformByVtkMatrix(m);

        geometry->SetImageGeometry(true);
        outputNetwork->SetGeometry(geometry);

        // read network
        std::map< int, mitk::ConnectomicsNetwork::VertexDescriptorType > idToVertexMap;
        // read vertices
        pElem = hRoot.FirstChildElement(mitk::ConnectomicsNetworkDefinitions::XML_VERTICES).ToElement();
        {
          // walk through the vertices
          tinyxml2::XMLElement* vertexElement = pElem->FirstChildElement();

          for( ; vertexElement; vertexElement=vertexElement->NextSiblingElement())
          {
            std::vector< float > pos;
            std::string label;
            int vertexID(0);

            temp = vertexElement->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_X);
            pos.push_back(temp);
            temp = vertexElement->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_Y);
            pos.push_back(temp);
            temp = vertexElement->DoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_Z);
            pos.push_back(temp);
            vertexID = vertexElement->IntAttribute(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_ID);
            label = std::string(vertexElement->Attribute(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_LABEL));

            mitk::ConnectomicsNetwork::VertexDescriptorType newVertex = outputNetwork->AddVertex( vertexID );
            outputNetwork->SetLabel( newVertex, label );
            outputNetwork->SetCoordinates( newVertex, pos );

            if ( idToVertexMap.count( vertexID ) > 0 )
            {
              MITK_ERROR << "Aborting network creation, duplicate vertex ID in file.";
              return result;
            }
            idToVertexMap.insert( std::pair< int, mitk::ConnectomicsNetwork::VertexDescriptorType >( vertexID, newVertex) );
          }
        }

        // read edges
        pElem = hRoot.FirstChildElement(mitk::ConnectomicsNetworkDefinitions::XML_EDGES).ToElement();
        {
          // walk through the edges
          auto edgeElement = pElem->FirstChildElement();

          for( ; edgeElement; edgeElement=edgeElement->NextSiblingElement())
          {
            int edgeSourceID(0), edgeTargetID(0), edgeWeight(0);
            double edgeDoubleWeight(0.0);

			edgeSourceID = edgeElement->IntAttribute(mitk::ConnectomicsNetworkDefinitions::XML_EDGE_SOURCE_ID);
			edgeTargetID = edgeElement->IntAttribute(mitk::ConnectomicsNetworkDefinitions::XML_EDGE_TARGET_ID);
			edgeWeight = edgeElement->IntAttribute(mitk::ConnectomicsNetworkDefinitions::XML_EDGE_FIBERCOUNT_ID);			
			
            if(version == "0.1")
            {
              // in version 0.1 only the int weight was saved, assume double weight to be one by default
              edgeDoubleWeight = 1.0;
            }
            else
            {
              edgeElement->QueryDoubleAttribute(mitk::ConnectomicsNetworkDefinitions::XML_EDGE_DOUBLE_WEIGHT_ID, &edgeDoubleWeight);
            }

            mitk::ConnectomicsNetwork::VertexDescriptorType source = idToVertexMap.find( edgeSourceID )->second;
            mitk::ConnectomicsNetwork::VertexDescriptorType target = idToVertexMap.find( edgeTargetID )->second;
            outputNetwork->AddEdge( source, target, edgeSourceID, edgeTargetID, edgeWeight, edgeDoubleWeight);
          }
        }

        outputNetwork->UpdateBounds();
        result.push_back(outputNetwork.GetPointer());
        MITK_INFO << "Network read";
      }
      catch (mitk::Exception& e)
      {
        MITK_ERROR << e.GetDescription();
      }
      catch(...)
      {
        MITK_ERROR << "Unknown error occured while trying to read file.";
      }
    }

    return result;
  }


} //namespace MITK

mitk::ConnectomicsNetworkReader* mitk::ConnectomicsNetworkReader::Clone() const
{
  return new ConnectomicsNetworkReader(*this);
}
