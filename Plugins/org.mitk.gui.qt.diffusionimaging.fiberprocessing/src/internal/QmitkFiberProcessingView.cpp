/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchPart.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkFiberProcessingView.h"

#include <QMessageBox>
#include <QStringList>
#include <functional>

#include <mitkNodePredicateFunction.h>
#include <mitkNodePredicateProperty.h>
#include <mitkImageCast.h>
#include <mitkPointSet.h>
#include <mitkPlanarCircle.h>
#include <mitkPlanarPolygon.h>
#include <mitkPlanarRectangle.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkPlanarFigureRoiImage.h>
#include <mitkDataNodeObject.h>
#include <mitkTensorImage.h>
#include "usModuleRegistry.h"
#include <itkFiberCurvatureFilter.h>

#include "mitkNodePredicateDataType.h"
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>

#include <itkTractsToFiberEndingsImageFilter.h>
#include <itkTractDensityImageFilter.h>
#include <itkImageRegion.h>
#include <itkTractsToRgbaImageFilter.h>
#include <itkFiberExtractionFilter.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

const std::string QmitkFiberProcessingView::VIEW_ID = "org.mitk.views.fiberprocessing";
const std::string id_DataManager = "org.mitk.views.datamanager";
using namespace mitk;

QmitkFiberProcessingView::QmitkFiberProcessingView()
  : QmitkAbstractView()
  , m_Controls( 0 )
  , m_CircleCounter(0)
  , m_PolygonCounter(0)
  , m_UpsamplingFactor(1)
{

}

// Destructor
QmitkFiberProcessingView::~QmitkFiberProcessingView()
{
  RemoveObservers();
}

void QmitkFiberProcessingView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkFiberProcessingViewControls;
    m_Controls->setupUi( parent );

    connect( m_Controls->m_CircleButton, SIGNAL( clicked() ), this, SLOT( OnDrawCircle() ) );
    connect( m_Controls->m_PolygonButton, SIGNAL( clicked() ), this, SLOT( OnDrawPolygon() ) );
    connect(m_Controls->PFCompoANDButton, SIGNAL(clicked()), this, SLOT(GenerateAndComposite()) );
    connect(m_Controls->PFCompoORButton, SIGNAL(clicked()), this, SLOT(GenerateOrComposite()) );
    connect(m_Controls->PFCompoNOTButton, SIGNAL(clicked()), this, SLOT(GenerateNotComposite()) );
    connect(m_Controls->m_GenerateRoiImage, SIGNAL(clicked()), this, SLOT(GenerateRoiImage()) );

    connect(m_Controls->m_JoinBundles, SIGNAL(clicked()), this, SLOT(JoinBundles()) );
    connect(m_Controls->m_SubstractBundles, SIGNAL(clicked()), this, SLOT(SubtractBundles()) );
    connect(m_Controls->m_CopyBundle, SIGNAL(clicked()), this, SLOT(CopyBundles()) );

    connect(m_Controls->m_ExtractFibersButton, SIGNAL(clicked()), this, SLOT(Extract()));
    connect(m_Controls->m_RemoveButton, SIGNAL(clicked()), this, SLOT(Remove()));
    connect(m_Controls->m_ModifyButton, SIGNAL(clicked()), this, SLOT(Modify()));

    connect(m_Controls->m_ExtractionMethodBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
    connect(m_Controls->m_RemovalMethodBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
    connect(m_Controls->m_ModificationMethodBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateGui()));
    connect(m_Controls->m_ExtractionBoxMask, SIGNAL(currentIndexChanged(int)), this, SLOT(OnMaskExtractionChanged()));

    // FiberBundle::ColorFibersByScalarMap reads a 3D scalar image; DWIs, tensor, ODF, peak and
    // RGB images have several components or four dimensions.
    m_Controls->m_ColorMapBox->SetDataStorage(this->GetDataStorage());
    auto isScalarMap = mitk::NodePredicateFunction::New([](const mitk::DataNode* node) {
      auto image = dynamic_cast<const mitk::Image*>(node->GetData());
      return image!=nullptr && image->GetDimension()==3 && image->GetPixelType().GetNumberOfComponents()==1;
    });
    m_Controls->m_ColorMapBox->SetPredicate(isScalarMap);
  }

  UpdateGui();
  OnMaskExtractionChanged();
}

void QmitkFiberProcessingView::OnMaskExtractionChanged()
{
  m_Controls->m_FiberExtractionFractionLabel->setVisible(false);
  m_Controls->m_FiberExtractionFractionBox->setVisible(false);

  m_Controls->m_FiberExtractionThresholdLabel->setVisible(false);
  m_Controls->m_FiberExtractionThresholdBox->setVisible(false);

  m_Controls->m_InterpolateRoiBox->setVisible(false);
  m_Controls->m_BothEnds->setVisible(false);

  m_Controls->m_LabelsBox->setVisible(false);
  m_Controls->m_LabelsLabel->setVisible(false);

  if (m_Controls->m_ExtractionBoxMask->currentIndex() == 2 || m_Controls->m_ExtractionBoxMask->currentIndex() == 3)
  {
    m_Controls->m_FiberExtractionFractionLabel->setVisible(true);
    m_Controls->m_FiberExtractionFractionBox->setVisible(true);

    m_Controls->m_FiberExtractionThresholdLabel->setVisible(true);
    m_Controls->m_FiberExtractionThresholdBox->setVisible(true);

    m_Controls->m_InterpolateRoiBox->setVisible(true);
  }
  else if (m_Controls->m_ExtractionBoxMask->currentIndex() == 0 || m_Controls->m_ExtractionBoxMask->currentIndex() == 1)
  {
    if (m_Controls->m_ExtractionBoxMask->currentIndex() != 3)
      m_Controls->m_BothEnds->setVisible(true);
    m_Controls->m_InterpolateRoiBox->setVisible(true);

    m_Controls->m_FiberExtractionThresholdLabel->setVisible(true);
    m_Controls->m_FiberExtractionThresholdBox->setVisible(true);
  }
  else if (m_Controls->m_ExtractionBoxMask->currentIndex() == 4)
  {
    m_Controls->m_BothEnds->setVisible(true);

    m_Controls->m_LabelsBox->setVisible(true);
    m_Controls->m_LabelsLabel->setVisible(true);
  }
}

void QmitkFiberProcessingView::SetFocus()
{
  m_Controls->toolBoxx->setFocus();
}

void QmitkFiberProcessingView::Modify()
{
  switch (m_Controls->m_ModificationMethodBox->currentIndex())
  {
  case 0:
  {
    ResampleSelectedBundlesSpline();
    break;
  }
  case 1:
  {
    ResampleSelectedBundlesLinear();
    break;
  }
  case 2:
  {
    CompressSelectedBundles();
    break;
  }
  case 3:
  {
    MirrorFibers();
    break;
  }
  case 4:
  {
    WeightFibers();
    break;
  }
  case 5:
  {
    DoWeightColorCoding();
    break;
  }
  case 6:
  {
    DoCurvatureColorCoding();
    break;
  }
  case 7:
  {
    DoLengthColorCoding();
    break;
  }
  case 8:
  {
    DoImageColorCoding();
    break;
  }
  }
}

void QmitkFiberProcessingView::WeightFibers()
{
  float weight = this->m_Controls->m_BundleWeightBox->value();
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->SetFiberWeights(weight);
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::Remove()
{
  switch (m_Controls->m_RemovalMethodBox->currentIndex())
  {
  case 0:
  {
    RemoveDir();
    break;
  }
  case 1:
  {
    PruneBundle();
    break;
  }
  case 2:
  {
    ApplyCurvatureThreshold();
    break;
  }
  case 3:
  {
    RemoveWithMask(false);
    break;
  }
  case 4:
  {
    RemoveWithMask(true);
    break;
  }
  case 5:
  {
    ApplyWeightThreshold();
    break;
  }
  case 6:
  {
    ApplyDensityThreshold();
    break;
  }
  }
}

void QmitkFiberProcessingView::Extract()
{
  switch (m_Controls->m_ExtractionMethodBox->currentIndex())
  {
  case 0:
  {
    ExtractWithPlanarFigure();
    break;
  }
  case 1:
  {
    switch (m_Controls->m_ExtractionBoxMask->currentIndex())
    {
    {
    case 0:
        ExtractWithMask(true, false, false);
        break;
    }
    {
    case 1:
        ExtractWithMask(true, true, false);
        break;
    }
    {
    case 2:
        ExtractWithMask(false, false, false);
        break;
    }
    {
    case 3:
        ExtractWithMask(false, true, false);
        break;
    }
    {
    case 4:
        ExtractWithMask(true, false, true);
        break;
    }
    }
    break;
  }
  }
}

void QmitkFiberProcessingView::PruneBundle()
{
  int minLength = this->m_Controls->m_PruneFibersMinBox->value();
  int maxLength = this->m_Controls->m_PruneFibersMaxBox->value();
  if (minLength>maxLength)
  {
    QMessageBox::information(nullptr, "No output generated:", "The minimum length is larger than the maximum length.");
    return;
  }
  for (auto node : m_SelectedFB)
  {
    // on a copy: when the second step would leave no fibers, the bundle stays as it was
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData())->Clone();
    if (!fib->RemoveShortFibers(minLength) || !fib->RemoveLongFibers(maxLength))
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
    else
      node->SetData(fib);
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::ApplyWeightThreshold()
{
  float thr = this->m_Controls->m_WeightThresholdBox->value();
  std::vector< DataNode::Pointer > nodes = m_SelectedFB;
  for (auto node : nodes)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    mitk::FiberBundle::Pointer newFib = fib->FilterByWeights(thr);
    if (newFib->GetNumFibers()>0)
    {
      switch(m_Controls->m_LookupTableTypeBox->currentIndex())
      {
      case 0:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::JET);
        break;
      case 1:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::HOT_IRON);
        break;
      case 2:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::PLASMA);
        break;
      case 3:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::INFERNO);
        break;
      case 4:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::VIRIDIS);
        break;
      case 5:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::MAGMA);
        break;
      case 6:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::GRAYSCALE);
        break;
      case 7:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::MULTILABEL);
        break;
      default:
        newFib->ColorFibersByFiberWeights(false, mitk::LookupTable::JET);
      }

      node->SetData(newFib);
    }
    else
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::ApplyDensityThreshold()
{
  float thr = this->m_Controls->m_DensityThresholdBox->value();
  float ol = this->m_Controls->m_DensityOverlapBox->value();
  std::vector< DataNode::Pointer > nodes = m_SelectedFB;
  for (auto node : nodes)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    itk::TractDensityImageFilter< ItkFloatImageType >::Pointer generator = itk::TractDensityImageFilter< ItkFloatImageType >::New();
    generator->SetFiberBundle(fib);
    generator->SetMode(itk::TDI_MODE::DENSITY);
    generator->SetOutputAbsoluteValues(false);
    generator->Update();

    itk::FiberExtractionFilter<float>::Pointer extractor = itk::FiberExtractionFilter<float>::New();
    extractor->SetRoiImages({generator->GetOutput()});
    extractor->SetInputFiberBundle(fib);
    extractor->SetOverlapFraction(ol);
    extractor->SetInterpolate(true);
    extractor->SetThreshold(thr);
    extractor->SetNoNegatives(true);
    extractor->Update();

    if (extractor->GetPositives().empty())
    {
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
      continue;
    }

    mitk::FiberBundle::Pointer newFib = extractor->GetPositives().at(0);
    if (newFib->GetNumFibers()>0)
      node->SetData(newFib);
    else
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::ApplyCurvatureThreshold()
{
  double angle = this->m_Controls->m_CurvSpinBox->value();
  double dist = this->m_Controls->m_CurvDistanceSpinBox->value();
  std::vector< DataNode::Pointer > nodes = m_SelectedFB;
  for (auto node : nodes)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    itk::FiberCurvatureFilter::Pointer filter = itk::FiberCurvatureFilter::New();
    filter->SetInputFiberBundle(fib);
    filter->SetAngularDeviation(angle);
    filter->SetDistance(dist);
    filter->SetRemoveFibers(m_Controls->m_RemoveCurvedFibersBox->isChecked());
    filter->Update();
    mitk::FiberBundle::Pointer newFib = filter->GetOutputFiberBundle();
    if (newFib->GetNumFibers()>0)
    {
      newFib->ColorFibersByOrientation();
      node->SetData(newFib);
    }
    else
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::RemoveDir()
{
  vnl_vector_fixed<double,3> dir;
  dir[0] = m_Controls->m_ExtractDirX->value();
  dir[1] = m_Controls->m_ExtractDirY->value();
  dir[2] = m_Controls->m_ExtractDirZ->value();
  // a zero vector would be normalized to zero: no fiber would ever match
  if (dir.magnitude()==0)
  {
    QMessageBox::information(nullptr, "Remove fibers in direction", "Enter a direction other than (0, 0, 0).");
    return;
  }

  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->RemoveDir(dir,cos((float)m_Controls->m_ExtractAngle->value()*itk::Math::pi/180));
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::RemoveWithMask(bool removeInside)
{
  if (m_RoiImageNode.IsNull())
    return;

  mitk::Image::Pointer mitkMask = dynamic_cast<mitk::Image*>(m_RoiImageNode->GetData());
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    ItkUCharImageType::Pointer mask = ItkUCharImageType::New();
    mitk::CastToItkImage(mitkMask, mask);
    mitk::FiberBundle::Pointer newFib = fib->RemoveFibersOutside(mask, removeInside);
    if (newFib->GetNumFibers()<=0)
    {
      QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
      continue;
    }
    node->SetData(newFib);
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::ExtractWithMask(bool onlyEnds, bool invert, bool labelmap)
{
  // Adding the results to the data storage may change the selection (NodeAdded), so the loop
  // works on copies.
  const mitk::DataNode::Pointer roiNode = m_RoiImageNode;
  const auto bundleNodes = m_SelectedFB;
  if (roiNode.IsNull())
    return;

  mitk::Image::Pointer mitkMask = dynamic_cast<mitk::Image*>(roiNode->GetData());
  const std::string roi_name = roiNode->GetName();
  for (auto node : bundleNodes)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    ItkFloatImageType::Pointer mask = ItkFloatImageType::New();
    mitk::CastToItkImage(mitkMask, mask);

    itk::FiberExtractionFilter<float>::Pointer extractor = itk::FiberExtractionFilter<float>::New();
    extractor->SetInputFiberBundle(fib);
    extractor->SetRoiImages({mask});
    extractor->SetRoiImageNames({roi_name});
    extractor->SetThreshold(m_Controls->m_FiberExtractionThresholdBox->value());
    extractor->SetOverlapFraction(m_Controls->m_FiberExtractionFractionBox->value());
    extractor->SetBothEnds(m_Controls->m_BothEnds->isChecked());
    extractor->SetInterpolate(m_Controls->m_InterpolateRoiBox->isChecked());
    extractor->SetMinFibersPerTract(m_Controls->m_MinExtractedFibersBox->value());
    extractor->SetSplitByRoi(true);
    if (invert)
      extractor->SetNoPositives(true);
    else
      extractor->SetNoNegatives(true);

    if (labelmap)
    {
      std::string labels_string = m_Controls->m_LabelsBox->text().toStdString();
      if (labels_string!="ALL")
      {
        std::vector<std::string> strs;
        boost::split(strs,labels_string,boost::is_any_of(" ,;\t"));

        std::vector< unsigned short > labels_vector;
        QStringList invalid;
        for (auto v : strs)
        {
          if (v.empty())
            continue;
          try{
            unsigned short l = boost::lexical_cast<unsigned short>(v);
            labels_vector.push_back(l);
          }
          catch(...)
          {
            invalid.push_back(QString::fromStdString(v));
          }
        }
        // an empty label list means "all labels" to the extractor
        if (!invalid.empty() || labels_vector.empty())
        {
          QMessageBox::information(nullptr, "No output generated:", QString("Invalid labels: \"%1\". Enter label values separated by spaces, commas or semicolons, or ALL.").arg(invalid.join(" ")));
          return;
        }
        extractor->SetLabels(labels_vector);
      }

      extractor->SetInterpolate(false);
      extractor->SetInputType(itk::FiberExtractionFilter<float>::INPUT::LABEL_MAP);
      extractor->SetSplitLabels(true);
      onlyEnds = true;
    }

    if (onlyEnds)
      extractor->SetMode(itk::FiberExtractionFilter<float>::MODE::ENDPOINTS);
    extractor->Update();

    std::vector< mitk::FiberBundle::Pointer > newFibs;
    if (invert)
      newFibs = extractor->GetNegatives();
    else
      newFibs = extractor->GetPositives();

    if (newFibs.empty())
    {
      QMessageBox::information(nullptr, "No output generated:", "No fibers could be extracted.");
      continue;
    }

    auto labels = extractor->GetPositiveLabels();
    for (unsigned int i=0; i<newFibs.size(); ++i)
    {
      DataNode::Pointer newNode = DataNode::New();
      auto fib = newFibs.at(i);

      newNode->SetData(fib);
      std::string name = roi_name;
      if (i<labels.size())
        name = labels.at(i);

      if (invert)
      {
        if (onlyEnds)
          name = "not-ending-in-" + name;
        else
          name = "not-passing-" + name;
      }
      else if (!labelmap)
      {
        if (onlyEnds)
          name = "ending-in-" + name;
        else
          name = "passing-" + name;
      }

      float currentThickness = 0;
      node->GetFloatProperty("Fiber2DSliceThickness", currentThickness);
      newNode->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(currentThickness));

      GetDataStorage()->Add(newNode, node);
    }
    node->SetVisibility(false);
  }
}

void QmitkFiberProcessingView::GenerateRoiImage()
{
  if (m_SelectedPF.empty())
    return;

  // WritePfToImage paints every figure of a composite, so the ROI image is the union of the
  // figures; an AND or NOT composite would silently differ from the fiber extraction with it.
  std::function<bool(const mitk::DataNode*)> hasAndOrNot = [this, &hasAndOrNot](const mitk::DataNode* node) {
    auto composite = dynamic_cast<mitk::PlanarFigureComposite*>(node->GetData());
    if (composite==nullptr)
      return false;
    if (composite->getOperationType()!=mitk::PlanarFigureComposite::OR)
      return true;
    auto children = this->GetDataStorage()->GetDerivations(node);
    for (unsigned int i=0; i<children->Size(); i++)
      if (hasAndOrNot(children->at(i)))
        return true;
    return false;
  };
  for (const auto& node : m_SelectedPF)
  {
    if (hasAndOrNot(node))
    {
      QMessageBox::information(nullptr, "Generate ROI Image", "The ROI image is the union of the selected planar figures. AND and NOT composites cannot be converted; select the figures or an OR composite.");
      return;
    }
  }

  mitk::BaseGeometry::Pointer geometry;
  if (!m_SelectedFB.empty())
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.front()->GetData());
    geometry = fib->GetGeometry();
  }
  else if (m_SelectedImage)
    geometry = m_SelectedImage->GetGeometry();
  else
    return;

  itk::Vector<double,3> spacing = geometry->GetSpacing();
  spacing /= m_UpsamplingFactor;

  mitk::Point3D newOrigin = geometry->GetOrigin();
  mitk::Geometry3D::BoundsArrayType bounds = geometry->GetBounds();
  newOrigin[0] += bounds.GetElement(0);
  newOrigin[1] += bounds.GetElement(2);
  newOrigin[2] += bounds.GetElement(4);

  // Column i of the geometry matrix is direction i scaled by the original spacing along i.
  itk::Matrix<double, 3, 3> direction;
  itk::ImageRegion<3> imageRegion;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      direction[j][i] = geometry->GetMatrixColumn(i)[j]/geometry->GetSpacing()[i];
  imageRegion.SetSize(0, geometry->GetExtent(0)*m_UpsamplingFactor);
  imageRegion.SetSize(1, geometry->GetExtent(1)*m_UpsamplingFactor);
  imageRegion.SetSize(2, geometry->GetExtent(2)*m_UpsamplingFactor);

  auto roiImage = ItkUCharImageType::New();
  roiImage->SetSpacing( spacing );   // Set the image spacing
  roiImage->SetOrigin( newOrigin );     // Set the image origin
  roiImage->SetDirection( direction );  // Set the image direction
  roiImage->SetRegions( imageRegion );
  roiImage->Allocate();
  roiImage->FillBuffer( 0 );

  std::string name = m_SelectedPF.at(0)->GetName();
  WritePfToImage(m_SelectedPF.at(0), *roiImage);
  for (unsigned int i=1; i<m_SelectedPF.size(); i++)
  {
    name += "+";
    name += m_SelectedPF.at(i)->GetName();
    WritePfToImage(m_SelectedPF.at(i), *roiImage);
  }

  DataNode::Pointer node = DataNode::New();
  Image::Pointer tmpImage = Image::New();
  tmpImage->InitializeByItk(roiImage.GetPointer());
  tmpImage->SetVolume(roiImage->GetBufferPointer());
  node->SetData(tmpImage);
  node->SetName(name);
  this->GetDataStorage()->Add(node);
}

void QmitkFiberProcessingView::WritePfToImage(mitk::DataNode::Pointer node, ItkUCharImageType& image)
{
  if (auto figure = dynamic_cast<mitk::PlanarFigure*>(node->GetData()))
  {
    mitk::AddPlanarFigureToRoiImage(*figure, image);
  }
  else if (dynamic_cast<mitk::PlanarFigureComposite*>(node->GetData()))
  {
    DataStorage::SetOfObjects::ConstPointer children = GetDataStorage()->GetDerivations(node);
    for (unsigned int i=0; i<children->Size(); i++)
    {
      WritePfToImage(children->at(i), image);
    }
  }
}

void QmitkFiberProcessingView::UpdateGui()
{
  m_Controls->m_FibLabel->setText("<font color='red'>mandatory</font>");
  m_Controls->m_PfLabel->setText("<font color='grey'>needed for extraction</font>");
  m_Controls->m_InputData->setTitle("Please Select Input Data");

  m_Controls->m_RemoveButton->setEnabled(false);

  m_Controls->m_PlanarFigureButtonsFrame->setEnabled(false);
  m_Controls->PFCompoANDButton->setEnabled(false);
  m_Controls->PFCompoORButton->setEnabled(false);
  m_Controls->PFCompoNOTButton->setEnabled(false);

  m_Controls->m_GenerateRoiImage->setEnabled(false);
  m_Controls->m_ExtractFibersButton->setEnabled(false);
  m_Controls->m_ModifyButton->setEnabled(false);

  m_Controls->m_CopyBundle->setEnabled(false);
  m_Controls->m_JoinBundles->setEnabled(false);
  m_Controls->m_SubstractBundles->setEnabled(false);

  // disable alle frames
  m_Controls->m_BundleWeightFrame->setVisible(false);
  m_Controls->m_ExtactionFramePF->setVisible(false);
  m_Controls->m_RemoveDirectionFrame->setVisible(false);
  m_Controls->m_RemoveLengthFrame->setVisible(false);
  m_Controls->m_RemoveCurvatureFrame->setVisible(false);
  m_Controls->m_RemoveByWeightFrame->setVisible(false);
  m_Controls->m_RemoveByDensityFrame->setVisible(false);
  m_Controls->m_SmoothFibersFrame->setVisible(false);
  m_Controls->m_CompressFibersFrame->setVisible(false);
  m_Controls->m_ColorFibersFrame->setVisible(false);
  m_Controls->m_MirrorFibersFrame->setVisible(false);
  m_Controls->m_MaskExtractionFrame->setVisible(false);
  m_Controls->m_ColorMapBox->setVisible(false);
  m_Controls->m_ValueAsWeightBox->setVisible(false);

  bool pfSelected = !m_SelectedPF.empty();
  bool fibSelected = !m_SelectedFB.empty();
  bool multipleFibsSelected = (m_SelectedFB.size()>1);
  bool maskSelected = m_RoiImageNode.IsNotNull();
  bool imageSelected = m_SelectedImage.IsNotNull();

  // toggle visibility of elements according to selected method
  switch ( m_Controls->m_ExtractionMethodBox->currentIndex() )
  {
  case 0:
    m_Controls->m_ExtactionFramePF->setVisible(true);
    break;
  case 1:
    m_Controls->m_MaskExtractionFrame->setVisible(true);
    break;
  }

  switch ( m_Controls->m_RemovalMethodBox->currentIndex() )
  {
  case 0:
    m_Controls->m_RemoveDirectionFrame->setVisible(true);
    if ( fibSelected )
      m_Controls->m_RemoveButton->setEnabled(true);
    break;
  case 1:
    m_Controls->m_RemoveLengthFrame->setVisible(true);
    if ( fibSelected )
      m_Controls->m_RemoveButton->setEnabled(true);
    break;
  case 2:
    m_Controls->m_RemoveCurvatureFrame->setVisible(true);
    if ( fibSelected )
      m_Controls->m_RemoveButton->setEnabled(true);
    break;
  case 3:
    break;
  case 4:
    break;
  case 5:
    m_Controls->m_RemoveByWeightFrame->setVisible(true);
    if ( fibSelected )
      m_Controls->m_RemoveButton->setEnabled(true);
    break;
  case 6:
    m_Controls->m_RemoveByDensityFrame->setVisible(true);
    if ( fibSelected )
      m_Controls->m_RemoveButton->setEnabled(true);
    break;
  }

  switch ( m_Controls->m_ModificationMethodBox->currentIndex() )
  {
  case 0:
    m_Controls->m_SmoothFibersFrame->setVisible(true);
    break;
  case 1:
    m_Controls->m_SmoothFibersFrame->setVisible(true);
    break;
  case 2:
    m_Controls->m_CompressFibersFrame->setVisible(true);
    break;
  case 3:
    m_Controls->m_MirrorFibersFrame->setVisible(true);
    if (m_SelectedSurfaces.size()>0)
      m_Controls->m_ModifyButton->setEnabled(true);
    break;
  case 4:
    m_Controls->m_BundleWeightFrame->setVisible(true);
    break;
  case 5:
    m_Controls->m_ColorFibersFrame->setVisible(true);
    break;
  case 6:
    m_Controls->m_ValueAsWeightBox->setVisible(true);
    m_Controls->m_ColorFibersFrame->setVisible(true);
    break;
  case 7:
    m_Controls->m_ValueAsWeightBox->setVisible(true);
    m_Controls->m_ColorFibersFrame->setVisible(true);
    break;
  case 8:
    m_Controls->m_ValueAsWeightBox->setVisible(true);
    m_Controls->m_ColorFibersFrame->setVisible(true);
    m_Controls->m_ColorMapBox->setVisible(true);
    break;
  }

  // are fiber bundles selected?
  if ( fibSelected )
  {
    m_Controls->m_CopyBundle->setEnabled(true);
    m_Controls->m_ModifyButton->setEnabled(true);
    m_Controls->m_PlanarFigureButtonsFrame->setEnabled(true);
    m_Controls->m_FibLabel->setText(QString(m_SelectedFB.at(0)->GetName().c_str()));

    // one bundle and one planar figure needed to extract fibers
    if (pfSelected && m_Controls->m_ExtractionMethodBox->currentIndex()==0)
    {
      m_Controls->m_InputData->setTitle("Input Data");
      m_Controls->m_PfLabel->setText(QString(m_SelectedPF.at(0)->GetName().c_str()));
      m_Controls->m_ExtractFibersButton->setEnabled(true);
    }

    // more than two bundles needed to join/subtract
    if (multipleFibsSelected)
    {
      m_Controls->m_FibLabel->setText("multiple bundles selected");

      m_Controls->m_JoinBundles->setEnabled(true);
      m_Controls->m_SubstractBundles->setEnabled(true);
    }

    if (maskSelected && m_Controls->m_ExtractionMethodBox->currentIndex()==1)
    {
      m_Controls->m_InputData->setTitle("Input Data");
      m_Controls->m_PfLabel->setText(QString(m_RoiImageNode->GetName().c_str()));
      m_Controls->m_ExtractFibersButton->setEnabled(true);
    }

    if (maskSelected && (m_Controls->m_RemovalMethodBox->currentIndex()==3 || m_Controls->m_RemovalMethodBox->currentIndex()==4) )
    {
      m_Controls->m_InputData->setTitle("Input Data");
      m_Controls->m_PfLabel->setText(QString(m_RoiImageNode->GetName().c_str()));
      m_Controls->m_RemoveButton->setEnabled(true);
    }
  }

  // are planar figures selected?
  if (pfSelected)
  {
    if ( fibSelected || m_SelectedImage.IsNotNull())
      m_Controls->m_GenerateRoiImage->setEnabled(true);

    if (m_SelectedPF.size() > 1)
    {
      m_Controls->PFCompoANDButton->setEnabled(true);
      m_Controls->PFCompoORButton->setEnabled(true);
    }
    else
      m_Controls->PFCompoNOTButton->setEnabled(true);
  }

  // is image selected
  if (imageSelected || maskSelected)
  {
    m_Controls->m_PlanarFigureButtonsFrame->setEnabled(true);
  }
}

void QmitkFiberProcessingView::NodeRemoved(const mitk::DataNode* node )
{
  // Only the removed node leaves the selection. The data storage calls this before the node
  // leaves it, so the topmost fallback must skip the node.
  RemoveObservers();
  auto isRemoved = [node](const mitk::DataNode::Pointer& selected) { return selected.GetPointer() == node; };
  std::erase_if(m_SelectedFB, isRemoved);
  std::erase_if(m_SelectedPF, isRemoved);
  std::erase_if(m_SelectedSurfaces, [node](const mitk::Surface::Pointer& surface) { return surface.GetPointer() == node->GetData(); });
  if (m_RoiImageNode.GetPointer() == node)
    m_RoiImageNode = nullptr;
  if (m_SelectedImage.IsNotNull() && m_SelectedImage.GetPointer() == node->GetData())
    m_SelectedImage = nullptr;
  if (m_InteractiveNode.GetPointer() == node)
    m_InteractiveNode = nullptr;

  SelectTopmostIfEmpty(node);
  AddObservers();
  UpdateGui();
}

void QmitkFiberProcessingView::NodeAdded(const mitk::DataNode* )
{
  // A node added anywhere (a result of this view, a loaded file) keeps the selection, also the
  // bundles a slot of this view is looping over; it only fills an empty selection.
  if (!m_Controls->m_InteractiveBox->isChecked())
  {
    RemoveObservers();
    SelectTopmostIfEmpty(nullptr);
    AddObservers();
    UpdateGui();
  }
}

void QmitkFiberProcessingView::OnEndInteraction()
{
  if (m_Controls->m_InteractiveBox->isChecked())
    ExtractWithPlanarFigure(true);
}

void QmitkFiberProcessingView::AddObservers()
{
  typedef itk::SimpleMemberCommand< QmitkFiberProcessingView > SimpleCommandType;
  for (auto node : m_SelectedPF)
  {
    mitk::PlanarFigure* figure = dynamic_cast<mitk::PlanarFigure*>(node->GetData());
    if (figure!=nullptr)
    {
      // add observer for event when interaction with figure ends
      SimpleCommandType::Pointer endInteractionCommand = SimpleCommandType::New();
      endInteractionCommand->SetCallbackFunction( this, &QmitkFiberProcessingView::OnEndInteraction);
      m_FigureObservers.emplace_back(figure, figure->AddObserver( mitk::EndInteractionPlanarFigureEvent(), endInteractionCommand ));
    }
  }
}

void QmitkFiberProcessingView::RemoveObservers()
{
  // only the observers of this view: other views (e.g. Measurement) observe the same figures
  for (const auto& [figure, tag] : m_FigureObservers)
    figure->RemoveObserver(tag);
  m_FigureObservers.clear();
}

void QmitkFiberProcessingView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer>& nodes)
{
  RemoveObservers();

  //reset existing Vectors containing FiberBundles and PlanarFigures from a previous selection
  std::vector<mitk::DataNode::Pointer>  lastSelectedFB = m_SelectedFB;
  m_SelectedFB.clear();
  m_SelectedPF.clear();
  m_SelectedSurfaces.clear();
  m_SelectedImage = nullptr;
  m_RoiImageNode = nullptr;

  for (auto node: nodes)
  {
    if ( dynamic_cast<mitk::FiberBundle*>(node->GetData()) )
      m_SelectedFB.push_back(node);
    else if (dynamic_cast<mitk::PlanarPolygon*>(node->GetData()) || dynamic_cast<mitk::PlanarFigureComposite*>(node->GetData()) || dynamic_cast<mitk::PlanarCircle*>(node->GetData()))
      m_SelectedPF.push_back(node);
    else if (dynamic_cast<mitk::Image*>(node->GetData()))
    {
      m_SelectedImage = dynamic_cast<mitk::Image*>(node->GetData());
      if (m_SelectedImage->GetDimension()==3)
        m_RoiImageNode = node;
    }
    else if (dynamic_cast<mitk::Surface*>(node->GetData()))
      m_SelectedSurfaces.push_back(dynamic_cast<mitk::Surface*>(node->GetData()));
  }

  // if we perform interactive fiber extraction, we want to avoid auto-selection of the extracted bundle
  if (m_SelectedFB.empty() && m_Controls->m_InteractiveBox->isChecked())
    m_SelectedFB = lastSelectedFB;

  SelectTopmostIfEmpty(nullptr);

  AddObservers();
  UpdateGui();
}

void QmitkFiberProcessingView::SelectTopmostIfEmpty(const mitk::DataNode* excluded)
{
  auto isPlanarFigure = [](const mitk::DataNode* node) {
    return dynamic_cast<mitk::PlanarPolygon*>(node->GetData()) || dynamic_cast<mitk::PlanarFigureComposite*>(node->GetData()) || dynamic_cast<mitk::PlanarCircle*>(node->GetData());
  };

  // the top-level node of the given kind with the highest layer
  auto topmost = [this, excluded](auto isKind) {
    mitk::DataNode::Pointer result;
    int maxLayer = 0;
    auto nodes = this->GetDataStorage()->GetAll();
    for (unsigned int i=0; i<nodes->Size(); i++)
    {
      mitk::DataNode::Pointer node = nodes->at(i);
      if (node.GetPointer()==excluded || !isKind(node.GetPointer()) || this->GetDataStorage()->GetSources(node)->Size()>0)
        continue;

      int layer = 0;
      node->GetPropertyValue("layer", layer);
      if (layer>=maxLayer)
      {
        maxLayer = layer;
        result = node;
      }
    }
    return result;
  };

  // if no fibers or surfaces are selected, select topmost
  if (m_SelectedFB.empty() && m_SelectedSurfaces.empty())
  {
    auto node = topmost([](const mitk::DataNode* n) { return dynamic_cast<mitk::FiberBundle*>(n->GetData())!=nullptr; });
    if (node.IsNotNull())
      m_SelectedFB.push_back(node);
  }

  // if no planar figure is selected, select topmost
  if (m_SelectedPF.empty())
  {
    auto node = topmost(isPlanarFigure);
    if (node.IsNotNull())
      m_SelectedPF.push_back(node);
  }
}

void QmitkFiberProcessingView::OnDrawPolygon()
{
  mitk::PlanarPolygon::Pointer figure = mitk::PlanarPolygon::New();
  figure->ClosedOn();
  this->AddFigureToDataStorage(figure, QString("Polygon%1").arg(++m_PolygonCounter));
}

void QmitkFiberProcessingView::OnDrawCircle()
{
  mitk::PlanarCircle::Pointer figure = mitk::PlanarCircle::New();
  this->AddFigureToDataStorage(figure, QString("Circle%1").arg(++m_CircleCounter));
}

void QmitkFiberProcessingView::AddFigureToDataStorage(mitk::PlanarFigure* figure, const QString& name, const char *, mitk::BaseProperty* )
{
  // initialize figure's geometry with empty geometry
  mitk::PlaneGeometry::Pointer emptygeometry = mitk::PlaneGeometry::New();
  figure->SetPlaneGeometry( emptygeometry );

  //set desired data to DataNode where Planarfigure is stored
  mitk::DataNode::Pointer newNode = mitk::DataNode::New();
  newNode->SetName(name.toStdString());
  newNode->SetData(figure);
  newNode->SetBoolProperty("planarfigure.3drendering", true);
  newNode->SetBoolProperty("planarfigure.3drendering.fill", true);

  mitk::PlanarFigureInteractor::Pointer figureInteractor = dynamic_cast<mitk::PlanarFigureInteractor*>(newNode->GetDataInteractor().GetPointer());
  if(figureInteractor.IsNull())
  {
    figureInteractor = mitk::PlanarFigureInteractor::New();
    us::Module* planarFigureModule = us::ModuleRegistry::GetModule( "MitkPlanarFigure" );
    figureInteractor->LoadStateMachine("PlanarFigureInteraction.xml", planarFigureModule );
    figureInteractor->SetEventConfig( "PlanarFigureConfig.xml", planarFigureModule );
    figureInteractor->SetDataNode(newNode);
  }

  // figure drawn on the topmost layer / image
  GetDataStorage()->Add(newNode );

  RemoveObservers();
  for(unsigned int i = 0; i < m_SelectedPF.size(); i++)
    m_SelectedPF[i]->SetSelected(false);

  newNode->SetSelected(true);
  m_SelectedPF.clear();
  m_SelectedPF.push_back(newNode);
  AddObservers();
  UpdateGui();
}

void QmitkFiberProcessingView::ExtractWithPlanarFigure(bool interactive)
{
  if ( m_SelectedFB.empty() ){
    QMessageBox::information( nullptr, "Warning", "No fiber bundle selected!");
    return;
  }
  if ( m_SelectedPF.empty() ){
    QMessageBox::information( nullptr, "Warning", "No planar figure selected!");
    return;
  }

  try
  {
    std::vector<mitk::DataNode::Pointer> fiberBundles = m_SelectedFB;
    mitk::DataNode::Pointer planarFigure = m_SelectedPF.at(0);
    for (unsigned int i=0; i<fiberBundles.size(); i++)
    {
      mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(fiberBundles.at(i)->GetData());
      mitk::FiberBundle::Pointer extFB = fib->ExtractFiberSubset(planarFigure, GetDataStorage());

      float currentThickness = 0;
      fiberBundles.at(i)->GetFloatProperty("Fiber2DSliceThickness", currentThickness);

      if (interactive && m_Controls->m_InteractiveBox->isChecked())
      {
        if (m_InteractiveNode.IsNull())
        {

          m_InteractiveNode = mitk::DataNode::New();
          QString name("Interactive");
          m_InteractiveNode->SetName(name.toStdString());
          m_InteractiveNode->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(currentThickness));
          GetDataStorage()->Add(m_InteractiveNode);
        }
        // fainter for larger bundles, opaque below 25 fibers (and for an empty bundle)
        float op = fib->GetNumFibers()>25 ? 5.0/sqrt(fib->GetNumFibers()) : 1.0;
        float currentOp = 0;
        fiberBundles.at(i)->GetFloatProperty("opacity", currentOp);

        if (currentOp!=op)
        {
          fib->SetFiberColors(255, 255, 255);
          fiberBundles.at(i)->SetFloatProperty("opacity", op);
        }
        m_InteractiveNode->SetData(extFB);
      }
      else
      {
        if (extFB->GetNumFibers()<=0)
        {
          QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers.");
          continue;
        }

        mitk::DataNode::Pointer node;
        node = mitk::DataNode::New();
        node->SetData(extFB);
        QString name(fiberBundles.at(i)->GetName().c_str());
        name += "*";
        node->SetName(name.toStdString());
        node->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(currentThickness));
        fiberBundles.at(i)->SetVisibility(false);
        GetDataStorage()->Add(node);
      }
    }
  }
  catch(const std::exception& e)
  {
    // an unplaced figure has no control points: PlanarFigure throws itk::ExceptionObject,
    // the extraction mitk::Exception or std::out_of_range
    MITK_ERROR << e.what();
    QMessageBox::warning( nullptr, "Fiber extraction failed", "Did you only create the planar figure, using the circle or polygon button, but forgot to actually place it in the image afterwards? \nAfter creating a planar figure, simply left-click at the desired position in the image or on the tractogram to place it.");
  }
}

void QmitkFiberProcessingView::GenerateAndComposite()
{
  mitk::PlanarFigureComposite::Pointer PFCAnd = mitk::PlanarFigureComposite::New();
  PFCAnd->setOperationType(mitk::PlanarFigureComposite::AND);

  mitk::DataNode::Pointer newPFCNode;
  newPFCNode = mitk::DataNode::New();
  newPFCNode->SetName("AND");
  newPFCNode->SetData(PFCAnd);

  AddCompositeToDatastorage(newPFCNode, m_SelectedPF);
  RemoveObservers();
  m_SelectedPF.clear();
  m_SelectedPF.push_back(newPFCNode);
  AddObservers();
  UpdateGui();
}

void QmitkFiberProcessingView::GenerateOrComposite()
{
  mitk::PlanarFigureComposite::Pointer PFCOr = mitk::PlanarFigureComposite::New();
  PFCOr->setOperationType(mitk::PlanarFigureComposite::OR);

  mitk::DataNode::Pointer newPFCNode;
  newPFCNode = mitk::DataNode::New();
  newPFCNode->SetName("OR");
  newPFCNode->SetData(PFCOr);

  RemoveObservers();
  AddCompositeToDatastorage(newPFCNode, m_SelectedPF);
  m_SelectedPF.clear();
  m_SelectedPF.push_back(newPFCNode);
  UpdateGui();
}

void QmitkFiberProcessingView::GenerateNotComposite()
{
  mitk::PlanarFigureComposite::Pointer PFCNot = mitk::PlanarFigureComposite::New();
  PFCNot->setOperationType(mitk::PlanarFigureComposite::NOT);

  mitk::DataNode::Pointer newPFCNode;
  newPFCNode = mitk::DataNode::New();
  newPFCNode->SetName("NOT");
  newPFCNode->SetData(PFCNot);

  RemoveObservers();
  AddCompositeToDatastorage(newPFCNode, m_SelectedPF);
  m_SelectedPF.clear();
  m_SelectedPF.push_back(newPFCNode);
  AddObservers();
  UpdateGui();
}

void QmitkFiberProcessingView::AddCompositeToDatastorage(mitk::DataNode::Pointer pfc, std::vector<mitk::DataNode::Pointer> children, mitk::DataNode::Pointer parentNode )
{
  pfc->SetSelected(true);
  if (parentNode.IsNotNull())
    GetDataStorage()->Add(pfc, parentNode);
  else
    GetDataStorage()->Add(pfc);

  for (auto child : children)
  {
    if (dynamic_cast<PlanarFigure*>(child->GetData()))
    {
      mitk::DataNode::Pointer newChild;
      newChild = mitk::DataNode::New();
      newChild->SetData(dynamic_cast<PlanarFigure*>(child->GetData()));
      newChild->SetName( child->GetName() );
      newChild->SetBoolProperty("planarfigure.3drendering", true);
      newChild->SetBoolProperty("planarfigure.3drendering.fill", true);

      GetDataStorage()->Add(newChild, pfc);
      GetDataStorage()->Remove(child);
    }
    else if (dynamic_cast<PlanarFigureComposite*>(child->GetData()))
    {
      mitk::DataNode::Pointer newChild;
      newChild = mitk::DataNode::New();
      newChild->SetData(dynamic_cast<PlanarFigureComposite*>(child->GetData()));
      newChild->SetName( child->GetName() );

      std::vector< mitk::DataNode::Pointer > grandChildVector;
      mitk::DataStorage::SetOfObjects::ConstPointer grandchildren = GetDataStorage()->GetDerivations(child);
      for( mitk::DataStorage::SetOfObjects::const_iterator it = grandchildren->begin(); it != grandchildren->end(); ++it )
        grandChildVector.push_back(*it);

      AddCompositeToDatastorage(newChild, grandChildVector, pfc);
      GetDataStorage()->Remove(child);
    }
  }
  UpdateGui();
}

void QmitkFiberProcessingView::CopyBundles()
{
  if ( m_SelectedFB.empty() ){
    QMessageBox::information( nullptr, "Warning", "Select at least one fiber bundle!");
    MITK_WARN("QmitkFiberProcessingView") << "Select at least one fiber bundle!";
    return;
  }

  // adding the copies may change the selection (NodeAdded), so the loop works on a copy
  const auto bundleNodes = m_SelectedFB;
  for (auto node : bundleNodes)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    mitk::FiberBundle::Pointer newFib = fib->Clone();

    node->SetVisibility(false);
    QString name("");
    name += QString(node->GetName().c_str());
    name += "_copy";

    mitk::DataNode::Pointer fbNode = mitk::DataNode::New();
    fbNode->SetData(newFib);
    fbNode->SetName(name.toStdString());
    fbNode->SetVisibility(true);
    GetDataStorage()->Add(fbNode);
  }
  UpdateGui();
}

void QmitkFiberProcessingView::JoinBundles()
{
  if ( m_SelectedFB.size()<2 ){
    QMessageBox::information( nullptr, "Warning", "Select at least two fiber bundles!");
    MITK_WARN("QmitkFiberProcessingView") << "Select at least two fiber bundles!";
    return;
  }

  m_SelectedFB.at(0)->SetVisibility(false);
  mitk::FiberBundle::Pointer newBundle = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());

  std::vector< mitk::FiberBundle::Pointer > tractograms;
  for (unsigned int i=1; i<m_SelectedFB.size(); i++)
  {
    m_SelectedFB.at(i)->SetVisibility(false);
    tractograms.push_back(dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(i)->GetData()));
  }
  newBundle = newBundle->AddBundles(tractograms);

  mitk::DataNode::Pointer fbNode = mitk::DataNode::New();
  fbNode->SetData(newBundle);
  fbNode->SetName("Joined_Tractograms");
  fbNode->SetVisibility(true);
  GetDataStorage()->Add(fbNode);
  UpdateGui();
}

void QmitkFiberProcessingView::SubtractBundles()
{
  if ( m_SelectedFB.size()<2 ){
    QMessageBox::information( nullptr, "Warning", "Select at least two fiber bundles!");
    MITK_WARN("QmitkFiberProcessingView") << "Select at least two fiber bundles!";
    return;
  }

  mitk::FiberBundle::Pointer newBundle = dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(0)->GetData());
  m_SelectedFB.at(0)->SetVisibility(false);
  QString name("");
  name += QString(m_SelectedFB.at(0)->GetName().c_str());
  for (unsigned int i=1; i<m_SelectedFB.size(); i++)
  {
    newBundle = newBundle->SubtractBundle(dynamic_cast<mitk::FiberBundle*>(m_SelectedFB.at(i)->GetData()));
    if (newBundle.IsNull())
      break;
    name += "-"+QString(m_SelectedFB.at(i)->GetName().c_str());
    m_SelectedFB.at(i)->SetVisibility(false);
  }
  if (newBundle.IsNull())
  {
    QMessageBox::information(nullptr, "No output generated:", "The resulting fiber bundle contains no fibers. Did you select the fiber bundles in the correct order? X-Y is not equal to Y-X!");
    return;
  }

  mitk::DataNode::Pointer fbNode = mitk::DataNode::New();
  fbNode->SetData(newBundle);
  fbNode->SetName(name.toStdString());
  fbNode->SetVisibility(true);
  GetDataStorage()->Add(fbNode);
  UpdateGui();
}

void QmitkFiberProcessingView::ResampleSelectedBundlesSpline()
{
  double factor = this->m_Controls->m_SmoothFibersBox->value();
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->ResampleSpline(factor);
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::ResampleSelectedBundlesLinear()
{
  double factor = this->m_Controls->m_SmoothFibersBox->value();
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->ResampleLinear(factor);
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::CompressSelectedBundles()
{
  double factor = this->m_Controls->m_ErrorThresholdBox->value();
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->Compress(factor);
    fib->ColorFibersByOrientation();
  }
  RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberProcessingView::DoImageColorCoding()
{
  if (m_Controls->m_ColorMapBox->GetSelectedNode().IsNull())
  {
    QMessageBox::information(nullptr, "Bundle coloring aborted:", "No image providing the scalar values for coloring the selected bundle available.");
    return;
  }

  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    switch(m_Controls->m_LookupTableTypeBox->currentIndex())
    {
    case 0:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET, m_Controls->m_PercentileBox->value());
      break;
    case 1:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::HOT_IRON, m_Controls->m_PercentileBox->value());
      break;
    case 2:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::PLASMA, m_Controls->m_PercentileBox->value());
      break;
    case 3:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::INFERNO, m_Controls->m_PercentileBox->value());
      break;
    case 4:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::VIRIDIS, m_Controls->m_PercentileBox->value());
      break;
    case 5:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MAGMA, m_Controls->m_PercentileBox->value());
      break;
    case 6:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::GRAYSCALE, m_Controls->m_PercentileBox->value());
      break;
    case 7:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MULTILABEL, 1.0);
      break;
    default:
      fib->ColorFibersByScalarMap(dynamic_cast<mitk::Image*>(m_Controls->m_ColorMapBox->GetSelectedNode()->GetData()), m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET, m_Controls->m_PercentileBox->value());
    }

  }

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}


void QmitkFiberProcessingView::DoCurvatureColorCoding()
{
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    switch(m_Controls->m_LookupTableTypeBox->currentIndex())
    {
    case 0:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET);
      break;
    case 1:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::HOT_IRON);
      break;
    case 2:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::PLASMA);
      break;
    case 3:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::INFERNO);
      break;
    case 4:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::VIRIDIS);
      break;
    case 5:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MAGMA);
      break;
    case 6:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::GRAYSCALE);
      break;
    case 7:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MULTILABEL);
      break;
    default:
      fib->ColorFibersByCurvature(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET);
    }
  }

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}

void QmitkFiberProcessingView::DoLengthColorCoding()
{
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    switch(m_Controls->m_LookupTableTypeBox->currentIndex())
    {
    case 0:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET);
      break;
    case 1:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::HOT_IRON);
      break;
    case 2:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::PLASMA);
      break;
    case 3:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::INFERNO);
      break;
    case 4:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::VIRIDIS);
      break;
    case 5:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MAGMA);
      break;
    case 6:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::GRAYSCALE);
      break;
    case 7:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::MULTILABEL);
      break;
    default:
      fib->ColorFibersByLength(m_Controls->m_FiberOpacityBox->isChecked(), m_Controls->m_ValueAsWeightBox->isChecked(), mitk::LookupTable::JET);
    }
  }

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}

void QmitkFiberProcessingView::DoWeightColorCoding()
{
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());

    switch(m_Controls->m_LookupTableTypeBox->currentIndex())
    {
    case 0:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::JET);
      break;
    case 1:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::HOT_IRON);
      break;
    case 2:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::PLASMA);
      break;
    case 3:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::INFERNO);
      break;
    case 4:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::VIRIDIS);
      break;
    case 5:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::MAGMA);
      break;
    case 6:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::GRAYSCALE);
      break;
    case 7:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::MULTILABEL);
      break;
    default:
      fib->ColorFibersByFiberWeights(m_Controls->m_FiberOpacityBox->isChecked(), mitk::LookupTable::JET);
    }
  }

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}

void QmitkFiberProcessingView::MirrorFibers()
{
  unsigned int axis = this->m_Controls->m_MirrorFibersBox->currentIndex();
  for (auto node : m_SelectedFB)
  {
    mitk::FiberBundle::Pointer fib = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    fib->MirrorFibers(axis);
  }

  for (auto surf : m_SelectedSurfaces)
  {
    vtkSmartPointer<vtkPolyData> poly = surf->GetVtkPolyData();
    vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();

    for (int i=0; i<poly->GetNumberOfPoints(); i++)
    {
      double* point = poly->GetPoint(i);
      point[axis] *= -1;
      vtkNewPoints->InsertNextPoint(point);
    }
    poly->SetPoints(vtkNewPoints);
    surf->CalculateBoundingBox();
    surf->Modified();
  }

  if (auto renderWindowPart = this->GetRenderWindowPart())
  {
    renderWindowPart->RequestUpdate();
  }
}
