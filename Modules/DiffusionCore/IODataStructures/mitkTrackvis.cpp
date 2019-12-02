#include <mitkTrackvis.h>
#include <vtkTransformPolyDataFilter.h>

TrackVisFiberReader::TrackVisFiberReader()  { m_Filename = ""; m_FilePointer = nullptr; }

TrackVisFiberReader::~TrackVisFiberReader() { if (m_FilePointer) fclose( m_FilePointer ); }


// Create a TrackVis file and store standard metadata. The file is ready to append fibers.
// ---------------------------------------------------------------------------------------
short TrackVisFiberReader::create(std::string filename , const mitk::FiberBundle *fib, bool print_header)
{
  if (fib->GetTrackVisHeader().hdr_size > 0)
  {
    m_Header = fib->GetTrackVisHeader();
  }
  else
  {

    vtkSmartPointer< vtkMatrix4x4 > matrix = vtkSmartPointer< vtkMatrix4x4 >::New();
    matrix->Identity();

    for(int i=0; i<3 ;i++)
    {
      if (fib->GetReferenceGeometry().IsNotNull())
      {
        m_Header.dim[i]            = fib->GetReferenceGeometry()->GetExtent(i);
        m_Header.voxel_size[i]     = fib->GetReferenceGeometry()->GetSpacing()[i];
        m_Header.origin[i]         = fib->GetReferenceGeometry()->GetOrigin()[i];
        matrix = fib->GetReferenceGeometry()->GetVtkMatrix();
      }
      else
      {
        m_Header.dim[i]            = fib->GetGeometry()->GetExtent(i);
        m_Header.voxel_size[i]     = 1;
        m_Header.origin[i]         = 0;
      }
    }


    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        m_Header.vox_to_ras[i][j] = matrix->GetElement(i, j);

    m_Header.n_scalars = 0;
    m_Header.n_properties = 0;
    sprintf(m_Header.voxel_order,"LPS");
    m_Header.image_orientation_patient[0] = 1.0;
    m_Header.image_orientation_patient[1] = 0.0;
    m_Header.image_orientation_patient[2] = 0.0;
    m_Header.image_orientation_patient[3] = 0.0;
    m_Header.image_orientation_patient[4] = 1.0;
    m_Header.image_orientation_patient[5] = 0.0;
    m_Header.pad1[0] = 0;
    m_Header.pad1[1] = 0;
    m_Header.pad2[0] = 0;
    m_Header.pad2[1] = 0;
    m_Header.invert_x = 0;
    m_Header.invert_y = 0;
    m_Header.invert_z = 0;
    m_Header.swap_xy = 0;
    m_Header.swap_yz = 0;
    m_Header.swap_zx = 0;
    m_Header.n_count = 0;
    m_Header.version = 2;
    m_Header.hdr_size = 1000;
    std::string id = "TRACK";
    strcpy(m_Header.id_string, id.c_str());
  }

  if (print_header)
    this->print_header();

  // write the header to the file
  m_FilePointer = fopen(filename.c_str(),"w+b");
  if (m_FilePointer == nullptr)
  {
    printf("[ERROR] Unable to create file '%s'\n",filename.c_str());
    return 0;
  }
  if (fwrite((char*)&m_Header, 1, 1000, m_FilePointer) != 1000)
    MITK_ERROR << "TrackVis::create : Error occurding during writing fiber.";

  this->m_Filename = filename;

  return 1;
}


// Open an existing TrackVis file and read metadata information.
// The file pointer is positiond at the beginning of fibers data
// -------------------------------------------------------------
short TrackVisFiberReader::open( std::string filename )
{
  m_FilePointer = std::fopen(filename.c_str(), "rb");
  if (m_FilePointer == nullptr)
  {
    printf("[ERROR] Unable to open file '%s'\n",filename.c_str());
    return 0;
  }
  this->m_Filename = filename;
  return fread((char*)(&m_Header), 1, 1000, m_FilePointer);
}

short TrackVisFiberReader::write(const mitk::FiberBundle *fib)
{
  vtkSmartPointer<vtkPolyData> poly = fib->GetFiberPolyData();
  if (fib->GetReferenceGeometry().IsNotNull())
  {
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(fib->GetReferenceGeometry()->GetVtkTransform()->GetInverse());
    transformFilter->Update();
    poly = transformFilter->GetOutput();
  }

  for (unsigned int i=0; i<fib->GetNumFibers(); i++)
  {
    vtkCell* cell = poly->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    unsigned int     numSaved, pos = 0;
    //float* tmp = new float[3*maxSteps];
    std::vector< float > tmp;
    tmp.reserve(3*numPoints);

    numSaved = numPoints;
    for(unsigned int i=0; i<numSaved ;i++)
    {
      double* p = points->GetPoint(i);

      tmp[pos++] = p[0];
      tmp[pos++] = p[1];
      tmp[pos++] = p[2];
    }

    // write the coordinates to the file
    if ( fwrite((char*)&numSaved, 1, 4, m_FilePointer) != 4 )
    {
      printf( "[ERROR] Problems saving the fiber!\n" );
      return 1;
    }
    if ( fwrite((char*)&(tmp.front()), 1, 4*pos, m_FilePointer) != 4*pos )
    {
      printf( "[ERROR] Problems saving the fiber!\n" );
      return 1;
    }
  }

  return 0;
}

void TrackVisFiberReader::print_header()
{
  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << "see http://trackvis.org/docs/?subsect=fileformat" << std::endl;
  std::cout << "ONLY vox_to_ras AND voxel_order HEADER ENTRIES ARE USED FOR FIBER COORDINATE TRANSFORMATIONS!" << std::endl;

  std::cout << "\nid_string (should be \"TRACK\"): " << m_Header.id_string << std::endl;
  std::cout << "dim: [" << std::defaultfloat << m_Header.dim[0] << " " << m_Header.dim[1] << " " << m_Header.dim[2] << "]" << std::endl;
  std::cout << "voxel_size: [" << m_Header.voxel_size[0] << " " << m_Header.voxel_size[1] << " " << m_Header.voxel_size[2] << "]" << std::endl;
  std::cout << "origin: [" << m_Header.origin[0] << " " << m_Header.origin[1] << " " << m_Header.origin[2] << "]" << std::endl;
  std::cout << "vox_to_world: " << std::scientific << std::endl;
  std::cout << "[[" << m_Header.vox_to_ras[0][0] << ", " << m_Header.vox_to_ras[0][1] << ", " << m_Header.vox_to_ras[0][2] << ", " << m_Header.vox_to_ras[0][3] << "]" << std::endl;
  std::cout << " [" << m_Header.vox_to_ras[1][0] << ", " << m_Header.vox_to_ras[1][1] << ", " << m_Header.vox_to_ras[1][2] << ", " << m_Header.vox_to_ras[1][3] << "]" << std::endl;
  std::cout << " [" << m_Header.vox_to_ras[2][0] << ", " << m_Header.vox_to_ras[2][1] << ", " << m_Header.vox_to_ras[2][2] << ", " << m_Header.vox_to_ras[2][3] << "]" << std::endl;
  std::cout << " [" << m_Header.vox_to_ras[3][0] << ", " << m_Header.vox_to_ras[3][1] << ", " << m_Header.vox_to_ras[3][2] << ", " << m_Header.vox_to_ras[3][3] << "]]" << std::defaultfloat << std::endl;

  std::cout << "voxel_order: " << m_Header.voxel_order[0] << m_Header.voxel_order[1] << m_Header.voxel_order[2] << std::endl;
  std::cout << "pad1: " << m_Header.pad1[0] << m_Header.pad1[1] << std::endl;
  std::cout << "pad2: " << m_Header.pad2[0] << m_Header.pad2[1] << m_Header.pad2[2] << std::endl;
  std::cout << "image_orientation_patient: [" << m_Header.image_orientation_patient[0] << " " << m_Header.image_orientation_patient[1] << " " << m_Header.image_orientation_patient[2] << " " << m_Header.image_orientation_patient[3] << " " << m_Header.image_orientation_patient[4] << " " << m_Header.image_orientation_patient[5] << "]" << std::endl;
  std::cout << "invert_x: " << static_cast<bool>(m_Header.invert_x) << std::endl;
  std::cout << "invert_y: " << static_cast<bool>(m_Header.invert_y) << std::endl;
  std::cout << "invert_z: " << static_cast<bool>(m_Header.invert_z) << std::endl;
  std::cout << "swap_xy: " << static_cast<bool>(m_Header.swap_xy) << std::endl;
  std::cout << "swap_yz: " << static_cast<bool>(m_Header.swap_yz) << std::endl;
  std::cout << "swap_zx: " << static_cast<bool>(m_Header.swap_zx) << std::endl;

  std::cout << "n_count: " << m_Header.n_count << std::endl;
  std::cout << "version: " << m_Header.version << std::endl;
  std::cout << "hdr_size: " << m_Header.hdr_size << std::endl;
  std::cout << "\nNot printed: n_scalars, scalar_name, n_properties, property_name, reserved" << std::endl;

  std::cout << "--------------------------------------------------------" << std::endl;
}

short TrackVisFiberReader::read( mitk::FiberBundle* fib, bool use_matrix, bool print_header)
{
  int numPoints;
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  if (print_header)
    this->print_header();

  while (fread((char*)&numPoints, 1, 4, m_FilePointer)==4)
  {
    if ( numPoints <= 0 )
    {
      printf( "[ERROR] Trying to read a fiber with %d points!\n", numPoints );
      return -1;
    }
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

    float tmp[3];
    for(int i=0; i<numPoints; i++)
    {
      if (fread((char*)tmp, 1, 12, m_FilePointer) == 0)
        MITK_ERROR << "TrackVis::read: Error during read.";

      vtkIdType id = vtkNewPoints->InsertNextPoint(tmp);
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  fiberPolyData->SetPoints(vtkNewPoints);
  fiberPolyData->SetLines(vtkNewCells);

  mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
  vtkSmartPointer< vtkMatrix4x4 > matrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  matrix->Identity();

  if (use_matrix)
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
      {
        if (j<3)
          matrix->SetElement(i, j, m_Header.vox_to_ras[i][j]/m_Header.voxel_size[j]);
        else
          matrix->SetElement(i, j, m_Header.vox_to_ras[i][j]);
      }


  if (m_Header.voxel_order[0]=='R')
  {
    matrix->SetElement(0,0,-matrix->GetElement(0,0));
    matrix->SetElement(0,1,-matrix->GetElement(0,1));
    matrix->SetElement(0,2,-matrix->GetElement(0,2));
    matrix->SetElement(0,3,-matrix->GetElement(0,3));
  }
  if (m_Header.voxel_order[1]=='A')
  {
    matrix->SetElement(1,0,-matrix->GetElement(1,0));
    matrix->SetElement(1,1,-matrix->GetElement(1,1));
    matrix->SetElement(1,2,-matrix->GetElement(1,2));
    matrix->SetElement(1,3,-matrix->GetElement(1,3));
  }
  if (m_Header.voxel_order[2]=='I')
  {
    matrix->SetElement(2,0,-matrix->GetElement(2,0));
    matrix->SetElement(2,1,-matrix->GetElement(2,1));
    matrix->SetElement(2,2,-matrix->GetElement(2,2));
    matrix->SetElement(2,3,-matrix->GetElement(2,3));
  }
  geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputData(fiberPolyData);
  transformFilter->SetTransform(geometry->GetVtkTransform());
  transformFilter->Update();
  fib->SetFiberPolyData(transformFilter->GetOutput());
  fib->SetReferenceGeometry(dynamic_cast<mitk::BaseGeometry*>(geometry.GetPointer()));
  fib->SetTrackVisHeader(m_Header);

  return numPoints;
}



// Update the field in the header to the new FIBER TOTAL.
// ------------------------------------------------------
void TrackVisFiberReader::updateTotal( int totFibers )
{
  fseek(m_FilePointer, 1000-12, SEEK_SET);
  if (fwrite((char*)&totFibers, 1, 4, m_FilePointer) != 4)
    MITK_ERROR << "[ERROR] Problems saving the fiber!";
}


void TrackVisFiberReader::writeHdr()
{
  fseek(m_FilePointer, 0, SEEK_SET);
  if (fwrite((char*)&m_Header, 1, 1000, m_FilePointer) != 1000)
    MITK_ERROR << "[ERROR] Problems saving the fiber!";
}


// Close the TrackVis file, but keep the metadata in the header.
// -------------------------------------------------------------
void TrackVisFiberReader::close()
{
  fclose(m_FilePointer);
  m_FilePointer = nullptr;
}

bool TrackVisFiberReader::IsTransformValid()
{
  if (fabs(m_Header.image_orientation_patient[0])<=0.001 || fabs(m_Header.image_orientation_patient[3])<=0.001 || fabs(m_Header.image_orientation_patient[5])<=0.001)
    return false;
  return true;
}
