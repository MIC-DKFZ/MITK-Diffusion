/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkPreprocessingView.h"

#include <QFontDatabase>
#include <QMessageBox>
#include <QTableWidget>
#include <QTableWidgetItem>

#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkNodePredicateDataType.h>
#include <mitkPointSet.h>
#include <mitkProperties.h>

#include <vtkMatrix4x4.h>

#include <array>
#include <optional>

const std::string QmitkPreprocessingView::VIEW_ID = "org.mitk.views.diffusionpreprocessing";

namespace
{
  using PropHelper = mitk::DiffusionPropertyHelper;

  /**
   * The values of three spin boxes with 4 decimals; components the user did not edit keep their
   * exact current value. Empty if nothing was edited.
   */
  std::optional<mitk::Vector3D> ReadEditedVector(const mitk::Vector3D& current, const std::array<double, 3>& boxValues)
  {
    mitk::Vector3D result = current;
    bool edited = false;
    for (unsigned int i=0; i<3; ++i)
    {
      if (std::abs(boxValues[i] - current[i]) > 0.5e-4)
      {
        result[i] = boxValues[i];
        edited = true;
      }
    }
    return edited ? std::optional<mitk::Vector3D>(result) : std::nullopt;
  }

  /**
   * The values of a 3x3 table; entries whose text is still the displayed current value keep the
   * exact current value. Empty if nothing was edited; throws for text that is not a number.
   */
  std::optional<mitk::DiffusionImageEditing::MatrixType> ReadEditedMatrix(const QTableWidget* table, const mitk::DiffusionImageEditing::MatrixType& current)
  {
    mitk::DiffusionImageEditing::MatrixType result = current;
    bool edited = false;
    for (int r=0; r<3; ++r)
    {
      for (int c=0; c<3; ++c)
      {
        const QTableWidgetItem* item = table->item(r, c);
        if (item == nullptr || item->text() == QString::number(current[r][c]))
          continue;
        bool ok = false;
        const double value = item->text().toDouble(&ok);
        if (!ok)
          mitkThrow() << "\"" << item->text().toStdString() << "\" is not a number.";
        result[r][c] = value;
        edited = true;
      }
    }
    return edited ? std::optional<mitk::DiffusionImageEditing::MatrixType>(result) : std::nullopt;
  }

  /** Direction matrix of the image: columns are the world directions of the image axes. */
  mitk::DiffusionImageEditing::MatrixType ImageDirection(const mitk::BaseGeometry* geometry)
  {
    mitk::DiffusionImageEditing::MatrixType direction;
    for (int r=0; r<3; r++)
      for (int c=0; c<3; c++)
        direction[r][c] = geometry->GetVtkMatrix()->GetElement(r,c)/geometry->GetSpacing()[c];
    return direction;
  }

  /** The normalized gradients, one row per component, the values of a gradient in one column. */
  QString GradientText(const PropHelper::GradientDirectionsContainerType* gradients)
  {
    QString text;
    for (unsigned int component=0; component<3; ++component)
    {
      for (unsigned int j=0; j<gradients->Size(); j++)
      {
        auto g = gradients->at(j);
        g = g.normalize();
        if (g[component]>=0)
          text += " ";
        text += QString::number(g[component]);
        text += j<gradients->Size()-1 ? " " : "\n";
      }
    }
    return text;
  }

  void FillTable(QTableWidget* table, const mitk::DiffusionImageEditing::MatrixType& matrix)
  {
    for (int r=0; r<3; r++)
    {
      for (int c=0; c<3; c++)
      {
        auto item = new QTableWidgetItem(QString::number(matrix[r][c]));
        item->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
        table->setItem(r, c, item);
      }
    }
  }
}

QmitkPreprocessingView::QmitkPreprocessingView()
  : QmitkAbstractView()
{
}

QmitkPreprocessingView::~QmitkPreprocessingView() = default;

void QmitkPreprocessingView::CreateQtPartControl(QWidget *parent)
{
  if (m_Controls)
    return;

  m_Controls = std::make_unique<Ui::QmitkPreprocessingViewControls>();
  m_Controls->setupUi(parent);
  this->CreateConnections();

  m_Controls->m_MeasurementFrameTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_Controls->m_MeasurementFrameTable->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_Controls->m_DirectionMatrixTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_Controls->m_DirectionMatrixTable->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);

  this->OnImageSelectionChanged();
}

void QmitkPreprocessingView::SetFocus()
{
  m_Controls->m_MirrorGradientToHalfSphereButton->setFocus();
}

void QmitkPreprocessingView::CreateConnections()
{
  m_Controls->m_SelctedImageComboBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_MergeDwiBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_AlignImageBox->SetDataStorage(this->GetDataStorage());

  auto isMitkImage = mitk::TNodePredicateDataType<mitk::Image>::New();
  m_Controls->m_SelctedImageComboBox->SetPredicate(isMitkImage);
  m_Controls->m_MergeDwiBox->SetPredicate(isMitkImage);
  m_Controls->m_AlignImageBox->SetPredicate(isMitkImage);

  m_Controls->m_ResampleIntFrame->setVisible(false);

  connect(m_Controls->m_ButtonAverageGradients, &QAbstractButton::clicked, this, &QmitkPreprocessingView::AverageGradients);
  connect(m_Controls->m_ButtonExtractB0, &QAbstractButton::clicked, this, &QmitkPreprocessingView::ExtractB0);
  connect(m_Controls->m_ReduceGradientsButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoReduceGradientDirections);
  connect(m_Controls->m_ShowGradientsButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoShowGradientDirections);
  connect(m_Controls->m_MirrorGradientToHalfSphereButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoHalfSphereGradientDirections);
  connect(m_Controls->m_SwapGradientsButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoSwapGradientDimensions);
  connect(m_Controls->m_MergeDwisButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::MergeDwis);
  connect(m_Controls->m_CreateLengthCorrectedDwi, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoLengthCorrection);
  connect(m_Controls->m_ResampleImageButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoResampleImage);
  connect(m_Controls->m_ResampleTypeBox, &QComboBox::currentIndexChanged, this, &QmitkPreprocessingView::DoUpdateInterpolationGui);
  connect(m_Controls->m_CropImageButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoCropImage);
  connect(m_Controls->m_RemoveGradientButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoRemoveGradient);
  connect(m_Controls->m_ExtractGradientButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoExtractGradient);
  connect(m_Controls->m_FlipAxis, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoFlipAxis);
  connect(m_Controls->m_FlipGradientsButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoFlipGradientDirections);
  connect(m_Controls->m_ClearRotationButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoClearRotationOfGradients);
  connect(m_Controls->m_ModifyHeader, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoApplyHeader);
  connect(m_Controls->m_AlignImageButton, &QAbstractButton::clicked, this, &QmitkPreprocessingView::DoAlignImages);
  connect(m_Controls->m_SelctedImageComboBox, &QmitkDataStorageComboBox::OnSelectionChanged, this, &QmitkPreprocessingView::OnImageSelectionChanged);

  const QFont fixedFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
  m_Controls->m_OriginalGradientsText->setFont(fixedFont);
  m_Controls->m_WorkingGradientsText->setFont(fixedFont);
}

mitk::DataNode::Pointer QmitkPreprocessingView::SelectedNode() const
{
  return m_Controls->m_SelctedImageComboBox->GetSelectedNode();
}

mitk::Image::Pointer QmitkPreprocessingView::SelectedImage() const
{
  auto node = this->SelectedNode();
  return node.IsNotNull() ? dynamic_cast<mitk::Image*>(node->GetData()) : nullptr;
}

void QmitkPreprocessingView::AddResult(mitk::DataNode* parent, mitk::BaseData* data, const QString& name, bool reinitializeViews)
{
  auto node = mitk::DataNode::New();
  node->SetData(data);
  node->SetName(name.toStdString());
  this->GetDataStorage()->Add(node, parent);

  if (reinitializeViews)
    mitk::RenderingManager::GetInstance()->InitializeViews(data->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkPreprocessingView::ApplyEdit(const QString& description, const QString& suffix, const std::function<mitk::Image::Pointer(const mitk::Image*)>& edit, bool reinitializeViews)
{
  auto node = this->SelectedNode();
  auto image = this->SelectedImage();
  if (image.IsNull())
    return;

  mitk::Image::Pointer result;
  try
  {
    result = edit(image);
  }
  catch (const std::exception& e)
  {
    QMessageBox::warning(nullptr, "Warning", description + " failed: " + e.what());
    return;
  }
  this->AddResult(node, result, QString::fromStdString(node->GetName()) + suffix, reinitializeViews);
}

void QmitkPreprocessingView::DoAlignImages()
{
  auto targetNode = m_Controls->m_AlignImageBox->GetSelectedNode();
  if (targetNode.IsNull())
    return;
  auto targetImage = dynamic_cast<mitk::Image*>(targetNode->GetData());
  if (targetImage == nullptr)
    return;

  // Like every other edit, the result is a new image; only the origin is taken over.
  mitk::DiffusionImageEditing::HeaderChange change;
  change.origin = targetImage->GetGeometry()->GetOrigin();
  this->ApplyEdit("Aligning", "_aligned", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::ApplyHeader(image, change); }, true);
}

void QmitkPreprocessingView::DoFlipAxis()
{
  const std::array<bool, 3> flipAxes = { m_Controls->m_FlipX->isChecked(), m_Controls->m_FlipY->isChecked(), m_Controls->m_FlipZ->isChecked() };
  this->ApplyEdit("Flipping", "_flipped", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::FlipAxes(image, flipAxes); }, true);
}

void QmitkPreprocessingView::DoRemoveGradient()
{
  const unsigned int index = m_Controls->m_RemoveGradientBox->value();
  this->ApplyEdit("Removing the volume", "_removedgradients", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::RemoveVolumes(image, {index}); }, true);
}

void QmitkPreprocessingView::DoExtractGradient()
{
  const unsigned int index = m_Controls->m_ExtractGradientBox->value();
  this->ApplyEdit("Extracting the volume", "_direction-" + QString::number(index), [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::ExtractVolume(image, index); }, true);
}

void QmitkPreprocessingView::DoCropImage()
{
  const std::array<unsigned int, 3> lower = { static_cast<unsigned int>(m_Controls->m_XstartBox->value()), static_cast<unsigned int>(m_Controls->m_YstartBox->value()), static_cast<unsigned int>(m_Controls->m_ZstartBox->value()) };
  const std::array<unsigned int, 3> upper = { static_cast<unsigned int>(m_Controls->m_XendBox->value()), static_cast<unsigned int>(m_Controls->m_YendBox->value()), static_cast<unsigned int>(m_Controls->m_ZendBox->value()) };
  this->ApplyEdit("Cropping", "_cropped", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::Crop(image, lower, upper); }, true);
}

void QmitkPreprocessingView::DoUpdateInterpolationGui(int i)
{
  auto image = this->SelectedImage();
  if (image.IsNull())
    return;

  // Factor and spacing are doubles, the size integers; the boxes start at the current values.
  const bool size = i == 2;
  m_Controls->m_ResampleIntFrame->setVisible(size);
  m_Controls->m_ResampleDoubleFrame->setVisible(!size);

  const mitk::BaseGeometry* geometry = image->GetGeometry();
  if (i == 1)
  {
    m_Controls->m_ResampleDoubleX->setValue(geometry->GetSpacing()[0]);
    m_Controls->m_ResampleDoubleY->setValue(geometry->GetSpacing()[1]);
    m_Controls->m_ResampleDoubleZ->setValue(geometry->GetSpacing()[2]);
  }
  else if (size)
  {
    m_Controls->m_ResampleIntX->setValue(geometry->GetExtent(0));
    m_Controls->m_ResampleIntY->setValue(geometry->GetExtent(1));
    m_Controls->m_ResampleIntZ->setValue(geometry->GetExtent(2));
  }
}

void QmitkPreprocessingView::DoResampleImage()
{
  using mitk::DiffusionImageEditing::Interpolation;
  using mitk::DiffusionImageEditing::ResampleGrid;

  ResampleGrid grid = ResampleGrid::SamplingFactor;
  mitk::Vector3D value;
  switch (m_Controls->m_ResampleTypeBox->currentIndex())
  {
  case 1:
    grid = ResampleGrid::Spacing;
    [[fallthrough]];
  case 0:
    value[0] = m_Controls->m_ResampleDoubleX->value();
    value[1] = m_Controls->m_ResampleDoubleY->value();
    value[2] = m_Controls->m_ResampleDoubleZ->value();
    break;
  case 2:
    grid = ResampleGrid::Size;
    value[0] = m_Controls->m_ResampleIntX->value();
    value[1] = m_Controls->m_ResampleIntY->value();
    value[2] = m_Controls->m_ResampleIntZ->value();
    break;
  default:
    return;
  }

  const std::array<std::pair<Interpolation, QString>, 4> interpolations = {{
    {Interpolation::NearestNeighbour, "NearestNeighbour"},
    {Interpolation::Linear, "Linear"},
    {Interpolation::BSpline, "BSpline"},
    {Interpolation::WindowedSinc, "WindowedSinc"}}};
  const int index = m_Controls->m_InterpolatorBox->currentIndex();
  const auto& [interpolation, name] = interpolations[index >= 0 && index < 4 ? index : 0];

  this->ApplyEdit("Resampling", "_resampled_" + name, [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::Resample(image, grid, value, interpolation); }, false);
}

void QmitkPreprocessingView::DoApplyHeader()
{
  auto image = this->SelectedImage();
  if (image.IsNull())
    return;

  // Only edited fields are passed on: the spin boxes and tables show rounded values, and
  // re-applying them would also round what the user did not touch.
  mitk::DiffusionImageEditing::HeaderChange change;
  change.rotateGradientsWithImage = m_Controls->m_RotateGradientsWithImageBox->isChecked();

  const mitk::BaseGeometry* geometry = image->GetGeometry();
  change.spacing = ReadEditedVector(geometry->GetSpacing(), { m_Controls->m_HeaderSpacingX->value(), m_Controls->m_HeaderSpacingY->value(), m_Controls->m_HeaderSpacingZ->value() });
  const auto origin = ReadEditedVector(geometry->GetOrigin().GetVectorFromOrigin(), { m_Controls->m_HeaderOriginX->value(), m_Controls->m_HeaderOriginY->value(), m_Controls->m_HeaderOriginZ->value() });
  if (origin)
  {
    mitk::Point3D newOrigin;
    for (unsigned int i=0; i<3; ++i)
      newOrigin[i] = (*origin)[i];
    change.origin = newOrigin;
  }

  this->ApplyEdit("Applying the header", "_newheader", [&](const mitk::Image* input)
  {
    auto headerChange = change;
    headerChange.direction = ReadEditedMatrix(m_Controls->m_DirectionMatrixTable, ImageDirection(input->GetGeometry()));
    if (PropHelper::IsDiffusionWeightedImage(input))
      headerChange.measurementFrame = ReadEditedMatrix(m_Controls->m_MeasurementFrameTable, PropHelper::GetMeasurementFrame(input));
    return mitk::DiffusionImageEditing::ApplyHeader(input, headerChange);
  }, true);
}

void QmitkPreprocessingView::DoLengthCorrection()
{
  // The rounding value is only read here. The b-value table shows the image's shells unrounded,
  // because "Reduce" looks its rows up in the image's b-value map.
  const int roundingValue = m_Controls->m_B_ValueMap_Rounder_SpinBox->value();
  this->ApplyEdit("Rounding the b-values", "_rounded", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::RoundBValues(image, roundingValue); }, false);
}

void QmitkPreprocessingView::CleanBValueTableWidget()
{
  m_Controls->m_B_ValueMap_TableWidget->clear();
  m_Controls->m_B_ValueMap_TableWidget->setRowCount(1);
  m_Controls->m_B_ValueMap_TableWidget->setHorizontalHeaderLabels({"b-Value", "Number of gradients"});
  m_Controls->m_B_ValueMap_TableWidget->setItem(0, 0, new QTableWidgetItem("-"));
  m_Controls->m_B_ValueMap_TableWidget->setItem(0, 1, new QTableWidgetItem("-"));
}

void QmitkPreprocessingView::UpdateGradientDetails()
{
  auto image = this->SelectedImage();
  if (image.IsNull() || !PropHelper::IsDiffusionWeightedImage(image))
  {
    m_Controls->m_WorkingGradientsText->clear();
    m_Controls->m_OriginalGradientsText->clear();
    return;
  }

  m_Controls->m_WorkingGradientsText->setPlainText(GradientText(PropHelper::GetGradientContainer(image)));
  m_Controls->m_OriginalGradientsText->setPlainText(GradientText(PropHelper::GetOriginalGradientContainer(image)));
}

void QmitkPreprocessingView::UpdateBValueTableWidget()
{
  auto image = this->SelectedImage();
  if (image.IsNull() || !PropHelper::IsDiffusionWeightedImage(image))
  {
    this->CleanBValueTableWidget();
    return;
  }

  const auto& bValueMap = PropHelper::GetBValueMap(image);
  m_Controls->m_B_ValueMap_TableWidget->clear();
  m_Controls->m_B_ValueMap_TableWidget->setRowCount(bValueMap.size());
  m_Controls->m_B_ValueMap_TableWidget->setHorizontalHeaderLabels({"b-Value", "Number of gradients"});

  int row = 0;
  for (const auto& [bValue, indices] : bValueMap)
  {
    auto item = new QTableWidgetItem(QString::number(bValue));
    item->setFlags(item->flags() & ~Qt::ItemIsEditable);
    m_Controls->m_B_ValueMap_TableWidget->setItem(row, 0, item);
    m_Controls->m_B_ValueMap_TableWidget->setItem(row, 1, new QTableWidgetItem(QString::number(indices.size())));
    ++row;
  }
}

void QmitkPreprocessingView::ClearHeaderTables()
{
  for (int r=0; r<3; r++)
  {
    for (int c=0; c<3; c++)
    {
      m_Controls->m_MeasurementFrameTable->setItem(r, c, new QTableWidgetItem());
      m_Controls->m_DirectionMatrixTable->setItem(r, c, new QTableWidgetItem());
    }
  }
}

void QmitkPreprocessingView::OnImageSelectionChanged()
{
  auto node = this->SelectedNode();
  auto image = this->SelectedImage();
  m_ShownNode = node;

  // Without an image every action is off and nothing of an earlier image stays on display.
  const bool foundImage = image.IsNotNull();
  const bool foundDwi = foundImage && PropHelper::IsDiffusionWeightedImage(image);

  m_Controls->m_ButtonAverageGradients->setEnabled(foundDwi);
  m_Controls->m_ButtonExtractB0->setEnabled(foundDwi);
  m_Controls->m_CheckExtractAll->setEnabled(foundDwi);
  m_Controls->m_MeasurementFrameTable->setEnabled(foundDwi);
  m_Controls->m_ReduceGradientsButton->setEnabled(foundDwi);
  m_Controls->m_ShowGradientsButton->setEnabled(foundDwi);
  m_Controls->m_MirrorGradientToHalfSphereButton->setEnabled(foundDwi);
  m_Controls->m_SwapGradientsButton->setEnabled(foundDwi);
  m_Controls->m_MergeDwisButton->setEnabled(foundDwi);
  m_Controls->m_B_ValueMap_Rounder_SpinBox->setEnabled(foundDwi);
  m_Controls->m_CreateLengthCorrectedDwi->setEnabled(foundDwi);
  m_Controls->m_RemoveGradientButton->setEnabled(foundDwi);
  m_Controls->m_ExtractGradientButton->setEnabled(foundDwi);
  m_Controls->m_FlipGradientsButton->setEnabled(foundDwi);
  // Only images whose gradients include the image rotation (legacy NRRD, DICOM) have one to clear.
  m_Controls->m_ClearRotationButton->setEnabled(foundDwi && mitk::DiffusionImageEditing::AppliesImageRotationToGradients(image));

  m_Controls->m_DirectionMatrixTable->setEnabled(foundImage);
  m_Controls->m_ModifyHeader->setEnabled(foundImage);
  m_Controls->m_AlignImageBox->setEnabled(foundImage);
  m_Controls->m_AlignImageButton->setEnabled(foundImage);

  // Multi-component and 3D+t images other than DWIs are not supported by these operations.
  const bool singleVolume = foundDwi || (foundImage && image->GetPixelType().GetNumberOfComponents() == 1 && image->GetTimeSteps() == 1);
  m_Controls->m_FlipAxis->setEnabled(singleVolume);
  m_Controls->m_CropImageButton->setEnabled(singleVolume);
  m_Controls->m_ResampleImageButton->setEnabled(singleVolume);

  this->ClearHeaderTables();
  this->UpdateBValueTableWidget();
  this->UpdateGradientDetails();
  m_Controls->m_InputData->setTitle(foundImage ? "Input Data" : "Please Select Input Data");
  if (!foundImage)
    return;

  m_Controls->m_B_ValueMap_Rounder_SpinBox->setValue(1);
  this->DoUpdateInterpolationGui(m_Controls->m_ResampleTypeBox->currentIndex());

  const mitk::BaseGeometry* geometry = image->GetGeometry();
  m_Controls->m_HeaderSpacingX->setValue(geometry->GetSpacing()[0]);
  m_Controls->m_HeaderSpacingY->setValue(geometry->GetSpacing()[1]);
  m_Controls->m_HeaderSpacingZ->setValue(geometry->GetSpacing()[2]);
  m_Controls->m_HeaderOriginX->setValue(geometry->GetOrigin()[0]);
  m_Controls->m_HeaderOriginY->setValue(geometry->GetOrigin()[1]);
  m_Controls->m_HeaderOriginZ->setValue(geometry->GetOrigin()[2]);
  m_Controls->m_XstartBox->setMaximum(geometry->GetExtent(0)-1);
  m_Controls->m_YstartBox->setMaximum(geometry->GetExtent(1)-1);
  m_Controls->m_ZstartBox->setMaximum(geometry->GetExtent(2)-1);
  m_Controls->m_XendBox->setMaximum(geometry->GetExtent(0)-1);
  m_Controls->m_YendBox->setMaximum(geometry->GetExtent(1)-1);
  m_Controls->m_ZendBox->setMaximum(geometry->GetExtent(2)-1);

  FillTable(m_Controls->m_DirectionMatrixTable, ImageDirection(geometry));

  if (foundDwi)
  {
    FillTable(m_Controls->m_MeasurementFrameTable, PropHelper::GetMeasurementFrame(image));
    m_Controls->m_RemoveGradientBox->setMaximum(PropHelper::GetGradientContainer(image)->Size()-1);
    m_Controls->m_ExtractGradientBox->setMaximum(PropHelper::GetGradientContainer(image)->Size()-1);
  }
}

void QmitkPreprocessingView::Activated()
{
}

void QmitkPreprocessingView::Deactivated()
{
}

void QmitkPreprocessingView::Visible()
{
  // Edited header values and settings stay unless another image was chosen meanwhile.
  if (m_Controls && this->SelectedNode() != m_ShownNode)
    this->OnImageSelectionChanged();
}

void QmitkPreprocessingView::Hidden()
{
}

void QmitkPreprocessingView::DoClearRotationOfGradients()
{
  this->ApplyEdit("Clearing the rotation", "_OriginalGradients", [](const mitk::Image* image) { return mitk::DiffusionImageEditing::ClearGradientRotation(image); }, false);
}

void QmitkPreprocessingView::DoFlipGradientDirections()
{
  const std::array<bool, 3> components = { m_Controls->m_FlipGradBoxX->isChecked(), m_Controls->m_FlipGradBoxY->isChecked(), m_Controls->m_FlipGradBoxZ->isChecked() };
  this->ApplyEdit("Flipping the gradients", "_GradientFlip", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::FlipGradients(image, components); }, false);
}

void QmitkPreprocessingView::DoSwapGradientDimensions()
{
  std::array<unsigned int, 3> order = {0, 1, 2};
  const QStringList axes = m_Controls->m_SwapGradientsEdit->text().split(' ', Qt::SkipEmptyParts);
  bool valid = axes.size() == 3;
  for (int i=0; valid && i<3; ++i)
  {
    if (axes[i] == "x")
      order[i] = 0;
    else if (axes[i] == "y")
      order[i] = 1;
    else if (axes[i] == "z")
      order[i] = 2;
    else
      valid = false;
  }
  if (!valid)
  {
    QMessageBox::warning(nullptr, "Warning", "String should contain x, y and z in arbitrary order separated by a space.");
    return;
  }
  this->ApplyEdit("Reordering the gradient components", "_GradientDimSwapped", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::PermuteGradientComponents(image, order); }, false);
}

void QmitkPreprocessingView::DoHalfSphereGradientDirections()
{
  this->ApplyEdit("Mirroring the gradients", "_halfsphere", [](const mitk::Image* image) { return mitk::DiffusionImageEditing::MirrorGradientsToHalfSphere(image); }, false);
}

void QmitkPreprocessingView::DoShowGradientDirections()
{
  auto node = this->SelectedNode();
  auto image = this->SelectedImage();
  if (image.IsNull() || !PropHelper::IsDiffusionWeightedImage(image))
    return;

  int maxIndex = 0;
  unsigned int maxSize = image->GetDimension(0);
  if (maxSize<image->GetDimension(1))
  {
    maxSize = image->GetDimension(1);
    maxIndex = 1;
  }
  if (maxSize<image->GetDimension(2))
  {
    maxSize = image->GetDimension(2);
    maxIndex = 2;
  }

  auto originSet = mitk::PointSet::New();
  auto gradientContainer = PropHelper::GetGradientContainer(image);
  const mitk::BaseGeometry* geometry = image->GetGeometry();
  const mitk::Point3D center = geometry->GetCenter();
  int shellCount = 1;

  for (const auto& [bValue, indices] : PropHelper::GetBValueMap(image))
  {
    auto pointset = mitk::PointSet::New();
    for (unsigned int j=0; j<indices.size(); j++)
    {
      // The gradients are world vectors drawn around the image centre; the centre computed per
      // world axis from origin and spacing is wrong for a non-identity direction matrix.
      vnl_vector_fixed< double, 3 > v = gradientContainer->at(indices[j]);
      mitk::Point3D ip = center;
      for (unsigned int d=0; d<3; ++d)
        ip[d] += v[d]*maxSize*geometry->GetSpacing()[maxIndex]/2;

      if (v.magnitude()>mitk::eps)
        pointset->InsertPoint(j, ip);
      else if (originSet->IsEmpty())
        originSet->InsertPoint(j, ip);
    }
    if (bValue < mitk::eps)
      continue;

    auto newNode = mitk::DataNode::New();
    newNode->SetData(pointset);
    newNode->SetName(node->GetName() + "_Shell_" + QString::number(bValue).toStdString());
    newNode->SetProperty("pointsize", mitk::FloatProperty::New((float)maxSize / 50));

    const int b0 = shellCount % 2;
    const int b1 = shellCount%4 >= 2 ? 1 : 0;
    const int b2 = shellCount>4 ? 1 : 0;
    newNode->SetProperty("color", mitk::ColorProperty::New(b2, b1, b0));
    this->GetDataStorage()->Add(newNode, node);
    shellCount++;
  }

  auto newNode = mitk::DataNode::New();
  newNode->SetData(originSet);
  newNode->SetName(node->GetName() + "_Origin");
  newNode->SetProperty("pointsize", mitk::FloatProperty::New((float)maxSize/50));
  newNode->SetProperty("color", mitk::ColorProperty::New(1,1,1));
  this->GetDataStorage()->Add(newNode, node);
}

void QmitkPreprocessingView::DoReduceGradientDirections()
{
  // One row per shell of the image's b-value map: b-value and the number of directions to keep.
  std::map<unsigned int, unsigned int> directionsPerShell;
  QString suffix;
  for (int i=0; i<m_Controls->m_B_ValueMap_TableWidget->rowCount(); i++)
  {
    bool bOk = false;
    bool numberOk = false;
    const unsigned int bValue = m_Controls->m_B_ValueMap_TableWidget->item(i,0)->text().toUInt(&bOk);
    const unsigned int number = m_Controls->m_B_ValueMap_TableWidget->item(i,1)->text().toUInt(&numberOk);
    if (!bOk || !numberOk)
    {
      QMessageBox::warning(nullptr, "Warning", "The number of gradients to keep must be a whole number for every shell.");
      return;
    }
    directionsPerShell[bValue] = number;
    suffix += "_" + QString::number(number);
  }
  if (directionsPerShell.empty())
    return;

  const bool useFirstN = m_Controls->m_KeepFirstNBox->isChecked();
  this->ApplyEdit("Reducing the gradients", suffix, [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::ReduceGradients(image, directionsPerShell, useFirstN); }, false);

  // The table shows the shells of the selected image again, not the numbers just applied.
  this->UpdateBValueTableWidget();
}

void QmitkPreprocessingView::MergeDwis()
{
  auto node = this->SelectedNode();
  auto node2 = m_Controls->m_MergeDwiBox->GetSelectedNode();
  if (node2.IsNull())
    return;
  const mitk::Image* image2 = dynamic_cast<mitk::Image*>(node2->GetData());
  if (image2 == nullptr)
    return;

  // The result is a child of the selected image, named after both inputs.
  this->ApplyEdit("Merging", "+" + QString::fromStdString(node2->GetName()), [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::MergeDwis({ image, image2 }); }, false);
}

void QmitkPreprocessingView::ExtractB0()
{
  if (m_Controls->m_CheckExtractAll->isChecked())
  {
    // A reinitialization makes the time steps of the new image reachable in the render windows.
    this->ApplyEdit("Extracting the b0 volumes", "_B0_ALL", [](const mitk::Image* image) { return mitk::DiffusionImageEditing::ExtractB0Volumes(image); }, true);
    return;
  }
  this->ApplyEdit("Extracting the b0 volume", "_B0", [](const mitk::Image* image) { return mitk::DiffusionImageEditing::ExtractB0(image); }, false);
}

void QmitkPreprocessingView::AverageGradients()
{
  const double precision = m_Controls->m_Blur->value();
  this->ApplyEdit("Averaging the repetitions", "_averaged", [&](const mitk::Image* image) { return mitk::DiffusionImageEditing::AverageRepeatedGradients(image, precision); }, false);
}
