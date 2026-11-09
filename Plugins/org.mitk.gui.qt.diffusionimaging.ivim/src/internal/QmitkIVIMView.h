/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _QMITKIVIMVIEW_H_INCLUDED
#define _QMITKIVIMVIEW_H_INCLUDED

#include <QmitkAbstractView.h>
#include <mitkILifecycleAwarePart.h>

#include "ui_QmitkIVIMViewControls.h"

#include <itkVectorImage.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIvimFit.h>
#include <mitkKurtosisFit.h>
#include <QmitkSliceNavigationListener.h>
#include <QmitkPlotWidget.h>

#include <memory>
#include <optional>

/**
  \brief Fits the IVIM or the kurtosis model to the signal at the crosshair, or to the mean signal
  of a mask, plots it, and generates parameter maps of a DWI.
*/
class QmitkIVIMView : public QmitkAbstractView, public mitk::ILifecycleAwarePart
{
  Q_OBJECT

public:

  static const std::string VIEW_ID;

  QmitkIVIMView();
  ~QmitkIVIMView() override;

  using GradientDirectionType = mitk::DiffusionPropertyHelper::GradientDirectionType;
  using DirContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using VecImgType = itk::VectorImage<short, 3>;

  void CreateQtPartControl(QWidget *parent) override;

  void SetFocus() override;

  void Activated() override;
  void Deactivated() override;
  void Visible() override;
  void Hidden() override;

protected slots:

  void OnSliceChanged();

  /** Generates the parameter maps of the selected DWI. */
  void FittIVIMStart();
  void AutoThreshold();

  void OnModelTabChanged(int tab);

  void OnIvimFitChanged(int val);

  void DStarSlider();
  void BThreshSlider();
  void S0ThreshSlider();
  void ClipboardStatisticsButtonClicked();
  void ClipboardCurveButtonClicked();
  void SavePlotButtonClicked();

  void UpdateGui();

protected:

  void InitChart(const char* yAxisTitle);
  void AddSecondFitPlot();
  void RemoveSecondFitPlot();
  void ActivateSliceListener();

  /** The values of the sliders. */
  double GetDStar() const;
  double GetBThreshold() const;
  double GetSignalThreshold() const;

  mitk::ivim::FitParameters GetIvimFitParameters() const;
  /** Fits and plots the signal of the selected voxel or the ROI average. */
  void FitIvimVoxel(const vnl_vector<double>& signal, DirContainerType::ConstPointer dirs, double bval);

  mitk::kurtosis::FitConfiguration GetKurtosisFitConfiguration() const;
  /** Fits and plots the signal of the selected voxel or the ROI average. */
  void FitKurtosisVoxel(const vnl_vector<double>& signal, DirContainerType::ConstPointer dirs, double bval);

  /** Adds a map below the DWI with a jet lookup table and, if given, the level window [0, max]. */
  void AddMapNode(mitk::Image* map, const QString& name, std::optional<double> levelWindowMax, mitk::DataNode* parent);

  std::unique_ptr<Ui::QmitkIVIMViewControls> m_Controls;

  mitk::ivim::Protocol m_IvimProtocol;
  mitk::ivim::VoxelFit m_IvimFit;
  vnl_vector<double> m_KurtosisBValues;
  vnl_vector<double> m_KurtosisMeasurements;
  mitk::kurtosis::VoxelFit m_KurtosisFit;

  bool m_Visible = false;
  bool m_ListenerActive = false;

  /** Keeps the preview while the controls are set up and while a mask is selected, whose mean
      signal does not change with the crosshair. */
  bool m_HoldUpdate = false;
  QmitkSliceNavigationListener m_SliceChangeListener;

  unsigned int m_DPartCurveId = 0;
  unsigned int m_FitCurveId = 0;
  unsigned int m_SignalCurveId = 0;
  unsigned int m_SecondFitSignalCurveId = 0;
  QwtPlotCurve* m_SecondFitSignalCurve = nullptr;
};

#endif // _QMITKIVIMVIEW_H_INCLUDED

