/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkConnectomicsNetworkWriter.h"
#include "mitkConnectomicsNetworkDefinitions.h"
#include <tinyxml2.h>
#include "itksys/SystemTools.hxx"
#include "mitkDiffusionIOMimeTypes.h"

mitk::ConnectomicsNetworkWriter::ConnectomicsNetworkWriter()
  : AbstractFileWriter(mitk::ConnectomicsNetwork::GetStaticNameOfClass(), CustomMimeType( mitk::DiffusionIOMimeTypes::CONNECTOMICS_MIMETYPE() ), mitk::DiffusionIOMimeTypes::CONNECTOMICS_MIMETYPE_DESCRIPTION() )
{
  RegisterService();
}

mitk::ConnectomicsNetworkWriter::ConnectomicsNetworkWriter(const mitk::ConnectomicsNetworkWriter& other)
  : AbstractFileWriter(other)
{
}


mitk::ConnectomicsNetworkWriter::~ConnectomicsNetworkWriter()
{}

mitk::ConnectomicsNetworkWriter* mitk::ConnectomicsNetworkWriter::Clone() const
{
  return new ConnectomicsNetworkWriter(*this);
}

void mitk::ConnectomicsNetworkWriter::Write()
{
  MITK_INFO << "Writing connectomics network";
  InputType::ConstPointer input = dynamic_cast<const InputType*>(this->GetInput());
  if (input.IsNull() )
  {
    MITK_ERROR <<"Sorry, input to ConnectomicsNetworkWriter is nullptr!";
    return;
  }
  if ( this->GetOutputLocation().empty() )
  {
    MITK_ERROR << "Sorry, filename has not been set!" ;
    return ;
  }

  std::string ext = itksys::SystemTools::GetFilenameLastExtension(this->GetOutputLocation());
  ext = itksys::SystemTools::LowerCase(ext);

  // default extension is .cnf
  if(ext == "")
  {
    ext = ".cnf";
    this->SetOutputLocation(this->GetOutputLocation() + ext);
  }

  if (ext == ".cnf")
  {
    // Get geometry of the network
    mitk::BaseGeometry* geometry = input->GetGeometry();

    // Create XML document
    tinyxml2::XMLDocument documentXML;
    { // begin document
      //tinyxml2::XMLDeclaration* declXML = new tinyxml2::XMLDeclaration( "1.0", "", "" ); // TODO what to write here? encoding? etc....
      //documentXML.LinkEndChild( declXML );

      auto mainXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_CONNECTOMICS_FILE);
      mainXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_FILE_VERSION, mitk::ConnectomicsNetworkDefinitions::VERSION_STRING);
      //documentXML.LinkEndChild(mainXML);
	  documentXML.InsertFirstChild(mainXML);

      auto geometryXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_GEOMETRY);
      { // begin geometry
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XX, geometry->GetMatrixColumn(0)[0]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XY, geometry->GetMatrixColumn(0)[1]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_XZ, geometry->GetMatrixColumn(0)[2]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YX, geometry->GetMatrixColumn(1)[0]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YY, geometry->GetMatrixColumn(1)[1]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_YZ, geometry->GetMatrixColumn(1)[2]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZX, geometry->GetMatrixColumn(2)[0]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZY, geometry->GetMatrixColumn(2)[1]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_MATRIX_ZZ, geometry->GetMatrixColumn(2)[2]);

        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_X, geometry->GetOrigin()[0]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_Y, geometry->GetOrigin()[1]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_ORIGIN_Z, geometry->GetOrigin()[2]);

        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_X, geometry->GetSpacing()[0]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_Y, geometry->GetSpacing()[1]);
        geometryXML->SetAttribute(mitk::ConnectomicsNetworkDefinitions::XML_SPACING_Z, geometry->GetSpacing()[2]);

      } // end geometry
      mainXML->InsertEndChild(geometryXML);

      auto verticesXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_VERTICES);
      { // begin vertices section
        VertexVectorType vertexVector = dynamic_cast<const InputType*>(this->GetInput())->GetVectorOfAllNodes();
        for( unsigned int index = 0; index < vertexVector.size(); index++ )
        {
          // not localized as of yet TODO
          auto vertexXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_VERTEX );
          vertexXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_ID , vertexVector[ index ].id );
          vertexXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_LABEL , vertexVector[ index ].label.c_str() );
          vertexXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_X , vertexVector[ index ].coordinates[0] );
          vertexXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_Y , vertexVector[ index ].coordinates[1] );
          vertexXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_VERTEX_Z , vertexVector[ index ].coordinates[2] );
          verticesXML->InsertEndChild(vertexXML);
        }
      } // end vertices section
      mainXML->InsertEndChild(verticesXML);

      auto edgesXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_EDGES);
      { // begin edges section
        EdgeVectorType edgeVector = dynamic_cast<const InputType*>(this->GetInput())->GetVectorOfAllEdges();
        for(unsigned  int index = 0; index < edgeVector.size(); index++ )
        {
          auto edgeXML = documentXML.NewElement(mitk::ConnectomicsNetworkDefinitions::XML_EDGE );
          edgeXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_EDGE_ID , index );
          edgeXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_EDGE_SOURCE_ID , edgeVector[ index ].second.sourceId );
          edgeXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_EDGE_TARGET_ID , edgeVector[ index ].second.targetId );
          edgeXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_EDGE_FIBERCOUNT_ID , edgeVector[ index ].second.fiber_count );
          edgeXML->SetAttribute( mitk::ConnectomicsNetworkDefinitions::XML_EDGE_DOUBLE_WEIGHT_ID , edgeVector[ index ].second.edge_weight );
          edgesXML->InsertEndChild(edgeXML);
        }
      } // end edges section
      mainXML->InsertEndChild(edgesXML);

    } // end document
    documentXML.SaveFile( this->GetOutputLocation().c_str() );
    MITK_INFO << "Connectomics network written";

  }
}
