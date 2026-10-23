/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkDiffusionDicomImportView.h"

// qt includes
#include <QFileDialog>

// itk includes
#include "itkTimeProbesCollectorBase.h"
#include "itkGDCMSeriesFileNames.h"
#include "itksys/SystemTools.hxx"

// mitk includes
#include <mitkProgressTask.h>
#include "mitkStatusBar.h"
#include "mitkProperties.h"
#include "mitkRenderingManager.h"
#include "mitkMemoryUtilities.h"
#include "mitkIOUtil.h"

// diffusion module includes
#include "mitkImage.h"
#include <mitkDiffusionPropertyHelper.h>

#include "mitkDiffusionDICOMFileReader.h"

#include "gdcmDirectory.h"
#include "gdcmScanner.h"
#include "gdcmSorter.h"
#include "gdcmIPPSorter.h"
#include "gdcmAttribute.h"
#include "gdcmVersion.h"

#include <itksys/SystemTools.hxx>
#include <itksys/Directory.hxx>

#include <QMessageBox>
#include <QTreeView>
#include <mitkStringProperty.h>

const std::string QmitkDiffusionDicomImport::VIEW_ID = "org.mitk.views.diffusiondicomimport";

namespace
{
  /** Files of one input folder. Without "Recursive", the files in the folder itself. With
   *  "Recursive", the files in the folder and in all its subfolders at any depth; a prefix
   *  restricts this to the subfolders of the input folder whose name starts with it (and
   *  everything below them), e.g. one subfolder per acquisition. */
  gdcm::Directory::FilenamesType CollectDicomFiles(const std::string& folder, bool recursive, const std::string& prefix)
  {
    if (!recursive || prefix.empty())
    {
      gdcm::Directory directory;
      directory.Load(folder, recursive);
      return directory.GetFilenames();
    }

    gdcm::Directory::FilenamesType files;
    itksys::Directory root;
    if (!root.Load(folder))
      return files;
    for (unsigned long i = 0; i < root.GetNumberOfFiles(); ++i)
    {
      const std::string name = root.GetFile(i);
      const std::string path = folder + "/" + name;
      if (name == "." || name == ".." || !name.starts_with(prefix) || !itksys::SystemTools::FileIsDirectory(path))
        continue;
      gdcm::Directory directory;
      directory.Load(path, true);
      const auto& subfolderFiles = directory.GetFilenames();
      files.insert(files.end(), subfolderFiles.begin(), subfolderFiles.end());
    }
    return files;
  }
}


QmitkDiffusionDicomImport::QmitkDiffusionDicomImport(QObject* /*parent*/, const char* /*name*/)
  : m_Controls(nullptr)
{
}

QmitkDiffusionDicomImport::~QmitkDiffusionDicomImport()
{}

void QmitkDiffusionDicomImport::CreateQtPartControl(QWidget *parent)
{
  m_Parent = parent;
  if (m_Controls == nullptr)
  {
    m_Controls = new Ui::QmitkDiffusionDicomImportControls;
    m_Controls->setupUi(parent);
    this->CreateConnections();

    m_Controls->m_DicomLoadRecursiveCheckbox->setChecked(false);
    m_Controls->m_DicomLoadAverageDuplicatesCheckbox->setChecked(false);

    m_Controls->m_DicomLoadRecursiveCheckbox->setVisible(true);
    m_Controls->m_SubdirPrefixLineEdit->setVisible(false);
    m_Controls->m_SetPrefixButton->setVisible(false);
    m_Controls->m_ResetPrefixButton->setVisible(false);

    AverageClicked();
  }
}



void QmitkDiffusionDicomImport::CreateConnections()
{
  if ( m_Controls )
  {
    connect( m_Controls->m_AddFoldersButton, SIGNAL(clicked()), this, SLOT(DicomLoadAddFolderNames()) );
    connect( m_Controls->m_DeleteFoldersButton, SIGNAL(clicked()), this, SLOT(DicomLoadDeleteFolderNames()) );
    //connect( m_Controls->m_DicomLoadStartLoadButton, SIGNAL(clicked()), this, SLOT(DicomLoadStartLoad()) );
    connect( m_Controls->m_DicomLoadStartLoadButton, SIGNAL(clicked()), this, SLOT(NewDicomLoadStartLoad()) );
    connect( m_Controls->m_DicomLoadAverageDuplicatesCheckbox, SIGNAL(clicked()), this, SLOT(AverageClicked()) );
    connect( m_Controls->m_Remove, SIGNAL(clicked()), this, SLOT(Remove()) );
    connect( m_Controls->m_SetPrefixButton, SIGNAL(clicked()), this, SLOT(SetPrefixButtonPushed()));
    connect( m_Controls->m_ResetPrefixButton, SIGNAL(clicked()), this, SLOT(ResetPrefixButtonPushed()));
    connect( m_Controls->m_DicomLoadRecursiveCheckbox, SIGNAL(clicked()), this, SLOT(RecursiveSettingsChanged()) );
  }
}

void QmitkDiffusionDicomImport::SetFocus()
{
  m_Controls->textBrowser->setFocus();
}

void QmitkDiffusionDicomImport::RecursiveSettingsChanged()
{
  // The prefix is the text of the line edit when the import starts; Set Prefix only locks it.
  m_Controls->m_SubdirPrefixLineEdit->setVisible( m_Controls->m_DicomLoadRecursiveCheckbox->isChecked() );
  m_Controls->m_SetPrefixButton->setVisible( m_Controls->m_DicomLoadRecursiveCheckbox->isChecked() );
  m_Controls->m_SubdirPrefixLineEdit->clear();
  this->m_Controls->m_SubdirPrefixLineEdit->setEnabled(true);
  m_Controls->m_SetPrefixButton->setEnabled(true);
  m_Controls->m_ResetPrefixButton->setVisible(false);
}

void QmitkDiffusionDicomImport::SetPrefixButtonPushed()
{
  if( !this->m_Controls->m_ResetPrefixButton->isVisible() )
    this->m_Controls->m_ResetPrefixButton->setVisible(true);

  this->m_Controls->m_SubdirPrefixLineEdit->setEnabled(false);
  this->m_Controls->m_ResetPrefixButton->setEnabled(true);
  this->m_Controls->m_SetPrefixButton->setEnabled(false);
}

void QmitkDiffusionDicomImport::ResetPrefixButtonPushed()
{
  m_Controls->m_SubdirPrefixLineEdit->clear();

  this->m_Controls->m_SubdirPrefixLineEdit->setEnabled(true);

  this->m_Controls->m_ResetPrefixButton->setEnabled(false);
  this->m_Controls->m_SetPrefixButton->setEnabled(true);
}

void QmitkDiffusionDicomImport::Remove()
{
  int i = m_Controls->listWidget->currentRow();
  delete m_Controls->listWidget->takeItem(i);
}

void QmitkDiffusionDicomImport::AverageClicked()
{
  m_Controls->m_Blur->setEnabled(m_Controls->m_DicomLoadAverageDuplicatesCheckbox->isChecked());
}

void QmitkDiffusionDicomImport::DicomLoadDeleteFolderNames()
{
  m_Controls->listWidget->clear();
}

void QmitkDiffusionDicomImport::DicomLoadAddFolderNames()
{
  // SELECT FOLDER DIALOG
  QFileDialog w( m_Parent, QString("Select folders containing DWI data") );
  w.setFileMode( QFileDialog::Directory );
  w.setOption(QFileDialog::DontUseNativeDialog,true);

  QListView *l = w.findChild<QListView*>("listView");
  if (l) {
    l->setSelectionMode(QAbstractItemView::MultiSelection);
  }
  QTreeView *t = w.findChild<QTreeView*>();
  if (t) {
    t->setSelectionMode(QAbstractItemView::MultiSelection);
  }

  // RETRIEVE SELECTION
  if ( w.exec() != QDialog::Accepted )
    return;

  m_Controls->listWidget->addItems(w.selectedFiles());
}

bool SortBySeriesUID(gdcm::DataSet const & ds1, gdcm::DataSet const & ds2 )
{
  gdcm::Attribute<0x0020,0x000e> at1;
  at1.Set( ds1 );
  gdcm::Attribute<0x0020,0x000e> at2;
  at2.Set( ds2 );
  return at1 < at2;
}

bool SortByAcquisitionNumber(gdcm::DataSet const & ds1, gdcm::DataSet const & ds2 )
{
  gdcm::Attribute<0x0020,0x0012> at1;
  at1.Set( ds1 );
  gdcm::Attribute<0x0020,0x0012> at2;
  at2.Set( ds2 );
  return at1 < at2;
}

bool SortBySeqName(gdcm::DataSet const & ds1, gdcm::DataSet const & ds2 )
{
  gdcm::Attribute<0x0018, 0x0024> at1;
  at1.Set( ds1 );
  gdcm::Attribute<0x0018, 0x0024> at2;
  at2.Set( ds2 );

  std::string str1 = at1.GetValue().Trim();
  std::string str2 = at2.GetValue().Trim();
  return std::lexicographical_compare(str1.begin(), str1.end(),
                                      str2.begin(), str2.end() );
}

void QmitkDiffusionDicomImport::Status(QString status)
{
  mitk::StatusBar::GetInstance()->DisplayText(status.toLatin1());
  MITK_INFO << status.toStdString().c_str();
}

void QmitkDiffusionDicomImport::Status(std::string status)
{
  mitk::StatusBar::GetInstance()->DisplayText(status.c_str());
  MITK_INFO << status.c_str();
}

void QmitkDiffusionDicomImport::Status(const char* status)
{
  mitk::StatusBar::GetInstance()->DisplayText(status);
  MITK_INFO << status;
}

void QmitkDiffusionDicomImport::Error(QString status)
{
  mitk::StatusBar::GetInstance()->DisplayErrorText(status.toLatin1());
  MITK_ERROR << status.toStdString().c_str();
}

void QmitkDiffusionDicomImport::Error(std::string status)
{
  mitk::StatusBar::GetInstance()->DisplayErrorText(status.c_str());
  MITK_ERROR << status.c_str();
}

void QmitkDiffusionDicomImport::Error(const char* status)
{
  mitk::StatusBar::GetInstance()->DisplayErrorText(status);
  MITK_ERROR << status;
}

void QmitkDiffusionDicomImport::PrintMemoryUsage()
{
  size_t processSize = mitk::MemoryUtilities::GetProcessMemoryUsage();
  size_t totalSize =  mitk::MemoryUtilities::GetTotalSizeOfPhysicalRam();
  float percentage = ( (float) processSize / (float) totalSize ) * 100.0;
  MITK_INFO << "Current memory usage: " << GetMemoryDescription( processSize, percentage );
}

std::string QmitkDiffusionDicomImport::FormatMemorySize( size_t size )
{
  double val = size;
  std::string descriptor("B");
  if ( val >= 1000.0 )
  {
    val /= 1024.0;
    descriptor = "KB";
  }
  if ( val >= 1000.0 )
  {
    val /= 1024.0;
    descriptor = "MB";
  }
  if ( val >= 1000.0 )
  {
    val /= 1024.0;
    descriptor = "GB";
  }
  std::ostringstream str;
  str << std::fixed << std::setprecision(2) << val << " " << descriptor;
  return str.str();
}

std::string QmitkDiffusionDicomImport::FormatPercentage( double val )
{
  std::ostringstream str;
  str << std::fixed << std::setprecision(2) << val << " " << "%";
  return str.str();
}

std::string QmitkDiffusionDicomImport::GetMemoryDescription( size_t processSize, float percentage )
{
  std::ostringstream str;
  str << FormatMemorySize(processSize) << " (" << FormatPercentage( percentage ) <<")" ;
  return str.str();
}

void QmitkDiffusionDicomImport::NewDicomLoadStartLoad()
{
  itk::TimeProbesCollectorBase clock;

  // LoadDiffusionDicomSeries switches to the C locale itself, for the time it parses.
  try
  {
    int nrFolders = m_Controls->listWidget->count();

    if(!nrFolders)
    {
      Error(QString("No input folders were selected. ABORTING."));
      return;
    }

    Status(QString("GDCM %1 used for DICOM parsing and sorting!").arg(gdcm::Version::GetVersion()));

    PrintMemoryUsage();
    QString status;
    mitk::DataNode::Pointer node;
    // Indeterminate: the import has no step to report; the task ends when this function returns.
    mitk::ProgressTask progress("Importing diffusion DICOM");

    // The folder list is cleared only after a successful import, so that a failed one can be
    // repeated with other options.
    gdcm::Directory::FilenamesType complete_list;
    for (int folderIndex = 0; folderIndex < nrFolders; ++folderIndex)
    {
      // RETREIVE FOLDERNAME
      QString folderName = m_Controls->listWidget->item(folderIndex)->text();
      const auto files = CollectDicomFiles(folderName.toStdString(),
                                           m_Controls->m_DicomLoadRecursiveCheckbox->isChecked(),
                                           m_Controls->m_SubdirPrefixLineEdit->text().toStdString());
      Status(QString(" ... found %1 different files").arg(files.size()));
      complete_list.insert(complete_list.end(), files.begin(), files.end());
    }

    if (complete_list.empty())
    {
      Error(QString("No files found in the selected folders."));
      QMessageBox::warning(m_Parent, "Diffusion DICOM import", "No files were found in the selected folders.");
      return;
    }

    {
      const auto images = mitk::DiffusionDICOMFileReader::LoadDiffusionDicomSeries( complete_list, this->m_Controls->m_SplitMosaicCheckBox->isChecked(), true );
      if (images.empty())
      {
        Error(QString("No diffusion-weighted image found in the selected folders."));
        QMessageBox::warning(m_Parent, "Diffusion DICOM import", "No diffusion-weighted image was found in the selected folders.");
        return;
      }
      for( const auto& loaded_image : images )
      {
        // "Merge duplicate gradients": volumes whose gradients differ by less than the given
        // precision are averaged, as in the preprocessing view.
        if (m_Controls->m_DicomLoadAverageDuplicatesCheckbox->isChecked())
        {
          mitk::DiffusionPropertyHelper::AverageRedundantGradients(loaded_image, m_Controls->m_Blur->value());
          mitk::DiffusionPropertyHelper::InitializeImage(loaded_image);
        }

        node = mitk::DataNode::New();
        node->SetData( loaded_image );
        node->SetName( loaded_image->GetProperty("name")->GetValueAsString() );

        GetDataStorage()->Add(node);
      }
    }

    Status("Timing information");
    clock.Report();

    if(node.IsNotNull())
    {
      mitk::BaseData::Pointer basedata = node->GetData();
      if (basedata.IsNotNull())
      {
        mitk::RenderingManager::GetInstance()->InitializeViews(
              basedata->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true );
      }
    }

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

    m_Controls->listWidget->clear();
  }
  catch (itk::ExceptionObject &ex)
  {
    Error(QString("%1\n%2\n%3\n%4\n%5\n%6").arg(ex.GetNameOfClass()).arg(ex.GetFile()).arg(ex.GetLine()).arg(ex.GetLocation()).arg(ex.what()).arg(ex.GetDescription()));
    return ;
  }
  catch (const std::exception &ex)
  {
    // e.g. std::bad_alloc for a series too large for the memory; the slot must not let it escape
    Error(QString("DICOM import failed: %1").arg(ex.what()));
    return;
  }

  Status(QString("Finished import with memory:"));

  PrintMemoryUsage();
}

void QmitkDiffusionDicomImport::SetDwiNodeProperties(mitk::DataNode::Pointer node, std::string name)
{

  node->SetProperty( "IsDWIRawVolume", mitk::BoolProperty::New( true ) );

  // set foldername as string property
  mitk::StringProperty::Pointer nameProp = mitk::StringProperty::New( name );
  node->SetProperty( "name", nameProp );
}
