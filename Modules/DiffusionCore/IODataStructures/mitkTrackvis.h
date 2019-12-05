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
#ifndef _TRACKVIS
#define _TRACKVIS

#include <mitkCommon.h>
#include <mitkFiberBundle.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <itkSize.h>

// Class to handle TrackVis files.
// -------------------------------
class MITKDIFFUSIONCORE_EXPORT TrackVisFiberReader
{
private:
    std::string         m_Filename;
    FILE*               m_FilePointer;
    mitk::FiberBundle::TrackVis_header m_Header;

public:

    short   create(std::string m_Filename, mitk::FiberBundle *fib, bool print_header);
    short   open(std::string m_Filename );
    short   read(mitk::FiberBundle* fib , bool use_matrix, bool print_header);
    short   write(const mitk::FiberBundle* fib );
    void    writeHdr();
    void    updateTotal( int totFibers );
    void    close();
    bool    IsTransformValid();
    void    print_header();

    TrackVisFiberReader();
    ~TrackVisFiberReader();
};

#endif
