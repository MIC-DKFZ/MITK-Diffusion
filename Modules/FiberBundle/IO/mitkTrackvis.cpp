#include <mitkTrackvis.h>
#include <mitkExceptionMacro.h>
#include <vtkTransformFilter.h>

#include <cstdint>
#include <cstring>

namespace
{
  /** Map from the center-based voxmm coordinates of a TrackVis file (voxmm - voxel_size/2) to
   *  MITK's LPS world, shared by reader and writer so that MITK reads back what it wrote.
   *  - Matrix recorded (version 2, vox_to_ras[3][3] != 0) and applied: by the TrackVis definition
   *    vox_to_ras maps voxel indices to RAS whatever the voxel order, so x and y are negated.
   *    Exception: MITK's own files with voxel order "LPS", n_count 0 and the default
   *    image_orientation_patient (written by default from 2019-12 to 2023-01 and since then with
   *    "Use LPS") store MITK's LPS index-to-world matrix; it is applied as it is.
   *  - No matrix (version 1 or vox_to_ras[3][3] == 0) or not applied: identity, with the axes
   *    flipped by voxel order (R, A, I), as MITK always read such files. */
  vtkSmartPointer<vtkMatrix4x4> VoxmmToLps(const mitk::FiberBundle::TrackVis_header& header, bool useMatrix)
  {
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->Identity();

    const bool hasMatrix = header.version >= 2 && header.vox_to_ras[3][3] != 0;
    const bool legacyMitkLps = std::strncmp(header.voxel_order, "LPS", 3) == 0 && header.n_count == 0
        && header.image_orientation_patient[0] == 1 && header.image_orientation_patient[1] == 0
        && header.image_orientation_patient[2] == 0 && header.image_orientation_patient[3] == 0
        && header.image_orientation_patient[4] == 1 && header.image_orientation_patient[5] == 0;

    bool negate[3] = {header.voxel_order[0]=='R', header.voxel_order[1]=='A', header.voxel_order[2]=='I'};
    if (useMatrix && hasMatrix)
    {
      for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
        {
          if (j<3)
            matrix->SetElement(i, j, header.vox_to_ras[i][j]/header.voxel_size[j]);
          else
            matrix->SetElement(i, j, header.vox_to_ras[i][j]);
        }
      negate[0] = !legacyMitkLps;
      negate[1] = !legacyMitkLps;
      negate[2] = false;
    }

    for (int i=0; i<3; ++i)
      if (negate[i])
        for (int j=0; j<4; ++j)
          matrix->SetElement(i, j, -matrix->GetElement(i, j));
    return matrix;
  }
}

TrackVisFiberReader::TrackVisFiberReader()  { m_Filename = ""; m_FilePointer = nullptr; }

TrackVisFiberReader::~TrackVisFiberReader() { if (m_FilePointer) fclose( m_FilePointer ); }

// Create a TrackVis file and store standard metadata. The file is ready to append fibers.
// ---------------------------------------------------------------------------------------
short TrackVisFiberReader::create(std::string filename , mitk::FiberBundle *fib, bool print_header, bool use_lps)
{
  m_Header = fib->GetTrackVisHeader();
  if (use_lps)
    sprintf(m_Header.voxel_order,"LPS");
  // Only coordinates are written, and readers such as nibabel read exactly n_count tracks; a
  // header copied from a loaded file would otherwise keep its counts.
  m_Header.n_scalars = 0;
  m_Header.n_properties = 0;
  m_Header.n_count = static_cast<int>(fib->GetNumFibers());

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
    mitkThrow() << "Unable to open TrackVis file " << filename;
  this->m_Filename = filename;

  if (fread((char*)(&m_Header), 1, 1000, m_FilePointer) != 1000)
    mitkThrow() << filename << " is shorter than a TrackVis header (1000 bytes).";
  if (std::strncmp(m_Header.id_string, "TRACK", 5) != 0)
    mitkThrow() << filename << " is no TrackVis file (id_string is not \"TRACK\").";
  if (m_Header.hdr_size != 1000)
  {
    // The TrackVis format detects the byte order by hdr_size.
    std::uint32_t size = static_cast<std::uint32_t>(m_Header.hdr_size);
    size = ((size >> 24) & 0xffu) | ((size >> 8) & 0xff00u) | ((size << 8) & 0xff0000u) | ((size << 24) & 0xff000000u);
    if (size == 1000)
      mitkThrow() << filename << " is a big-endian TrackVis file, which is not supported.";
    mitkThrow() << filename << " has an invalid TrackVis header size " << m_Header.hdr_size << " (expected 1000).";
  }
  if (m_Header.n_scalars < 0 || m_Header.n_properties < 0)
    mitkThrow() << filename << " has a negative number of scalars or properties in its TrackVis header.";
  return 1;
}

short TrackVisFiberReader::write(const mitk::FiberBundle *fib)
{
  vtkSmartPointer<vtkPolyData> poly = fib->GetFiberPolyData();
  {
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetIndexToWorldTransformByVtkMatrix(VoxmmToLps(m_Header, true));

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(poly);
    transformFilter->SetTransform(geometry->GetVtkTransform()->GetInverse());
    transformFilter->Update();
    poly = transformFilter->GetPolyDataOutput();
  }

  for (unsigned int i=0; i<fib->GetNumFibers(); i++)
  {
    vtkCell* cell = poly->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    unsigned int     numSaved, pos = 0;
    std::vector< float > tmp(3*static_cast<std::size_t>(numPoints));

    numSaved = numPoints;
    for(unsigned int i=0; i<numSaved ;i++)
    {
      double* p = points->GetPoint(i);

      // TRK coordinates are corner based, so we have to shift the center based vtk coordinates by half a voxel
      tmp[pos++] = p[0] + m_Header.voxel_size[0]/2;
      tmp[pos++] = p[1] + m_Header.voxel_size[1]/2;
      tmp[pos++] = p[2] + m_Header.voxel_size[2]/2;
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

short TrackVisFiberReader::read( mitk::FiberBundle* fib, bool use_matrix, bool flip_x, bool flip_y, bool flip_z, bool print_header)
{
  int numPoints = 0;
  int numFibers = 0;
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  if (print_header)
    this->print_header();

  // Each point holds x, y, z and n_scalars values, each track is followed by n_properties
  // values (TrackVis, nibabel, dipy); only the coordinates are imported.
  const std::size_t valuesPerPoint = 3 + static_cast<std::size_t>(m_Header.n_scalars);
  const std::size_t numProperties = static_cast<std::size_t>(m_Header.n_properties);
  if (valuesPerPoint > 3 || numProperties > 0)
    MITK_INFO << "Skipping " << m_Header.n_scalars << " scalars per point and " << m_Header.n_properties << " properties per track.";
  std::vector<float> pointValues(valuesPerPoint);
  std::vector<float> properties(numProperties);

  std::size_t countBytes = 0;
  while ((countBytes = fread((char*)&numPoints, 1, 4, m_FilePointer))==4)
  {
    if ( numPoints <= 0 )
      mitkThrow() << "TrackVis file " << m_Filename << ": track " << numFibers << " has " << numPoints << " points.";
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();

    float tmp[3];
    for(int i=0; i<numPoints; i++)
    {
      if (fread(pointValues.data(), sizeof(float), valuesPerPoint, m_FilePointer) != valuesPerPoint)
        mitkThrow() << "TrackVis file " << m_Filename << " is truncated in track " << numFibers << ".";

      // TRK coordinates are corner based, so we have to shift them back to center based coordinates
      tmp[0] = pointValues[0] - m_Header.voxel_size[0]/2;
      tmp[1] = pointValues[1] - m_Header.voxel_size[1]/2;
      tmp[2] = pointValues[2] - m_Header.voxel_size[2]/2;
      vtkIdType id = vtkNewPoints->InsertNextPoint(tmp);
      container->GetPointIds()->InsertNextId(id);
    }
    if (numProperties > 0 && fread(properties.data(), sizeof(float), numProperties, m_FilePointer) != numProperties)
      mitkThrow() << "TrackVis file " << m_Filename << " is truncated in track " << numFibers << ".";
    vtkNewCells->InsertNextCell(container);
    ++numFibers;
  }
  if (countBytes != 0)
    mitkThrow() << "TrackVis file " << m_Filename << " is truncated after track " << numFibers << ".";
  if (m_Header.n_count > 0 && m_Header.n_count != numFibers)
    MITK_WARN << "TrackVis header of " << m_Filename << " declares " << m_Header.n_count << " tracks, the file holds " << numFibers << ".";

  vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  fiberPolyData->SetPoints(vtkNewPoints);
  fiberPolyData->SetLines(vtkNewCells);

  mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
  vtkSmartPointer< vtkMatrix4x4 > matrix = VoxmmToLps(m_Header, use_matrix);

  if (flip_x)
  {
    matrix->SetElement(0,0,-matrix->GetElement(0,0));
    matrix->SetElement(0,1,-matrix->GetElement(0,1));
    matrix->SetElement(0,2,-matrix->GetElement(0,2));
    matrix->SetElement(0,3,-matrix->GetElement(0,3));
  }
  if (flip_y)
  {
    matrix->SetElement(1,0,-matrix->GetElement(1,0));
    matrix->SetElement(1,1,-matrix->GetElement(1,1));
    matrix->SetElement(1,2,-matrix->GetElement(1,2));
    matrix->SetElement(1,3,-matrix->GetElement(1,3));
  }
  if (flip_z)
  {
    matrix->SetElement(2,0,-matrix->GetElement(2,0));
    matrix->SetElement(2,1,-matrix->GetElement(2,1));
    matrix->SetElement(2,2,-matrix->GetElement(2,2));
    matrix->SetElement(2,3,-matrix->GetElement(2,3));
  }
  geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

  vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetInputData(fiberPolyData);
  transformFilter->SetTransform(geometry->GetVtkTransform());
  transformFilter->Update();
  fib->SetFiberPolyData(transformFilter->GetPolyDataOutput());
  fib->SetTrackVisHeader(geometry.GetPointer());
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
