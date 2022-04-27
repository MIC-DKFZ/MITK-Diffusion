/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/



#ifndef QmitkInteractiveFiberDissectionView_h
#define QmitkInteractiveFiberDissectionView_h
#include "ui_QmitkInteractiveFiberDissectionViewControls.h"

#include <QInputDialog> //Pointset
#include <QmitkPointListWidget.h>//Pointset

#include <QmitkAbstractView.h>//Pointset
#include <QmitkSingleNodeSelectionWidget.h>//Pointset


#include <mitkFiberBundle.h>
#include <mitkDataInteractor.h>
#include <mitkIRenderWindowPartListener.h> //Pointset

#include <itkCastImageFilter.h>
#include <itkVTKImageImport.h>
#include <itkVTKImageExport.h>
#include <itkRegionOfInterestImageFilter.h>

#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkImageImport.h>
#include <vtkImageExport.h>
#include <vtkImageStencil.h>
#include <vtkSmartPointer.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelectedThresholds.h>
#include <vtkFloatArray.h>
#include <vtkCellPicker.h>

/*!
\brief View to process fiber bundles. Supplies methods to extract fibers from the bundle, fiber resampling, mirroring, join and subtract bundles and much more.
*/
class QmitkInteractiveFiberDissectionView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  typedef itk::Image< unsigned char, 3 >    ItkUCharImageType;
  typedef itk::Image< float, 3 >            ItkFloatImageType;

  static const std::string VIEW_ID;

  QmitkInteractiveFiberDissectionView();
  virtual ~QmitkInteractiveFiberDissectionView();

  virtual void CreateQtPartControl(QWidget *parent) override;

  ///
  /// Sets the focus to an internal widget.
  ///
  virtual void SetFocus() override;

protected slots:

  void RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart) override; //Pointset
  void RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart) override; //Pointset
  void OnAddPointSetClicked();//Pointset
  void CreateStreamline();
//  void AddtoBundle();
  void RemovefromBundle();
  void ExtractRandomFibersFromTractogram();




  void UpdateGui();     ///< update button activity etc. dpending on current datamanager selection




protected:
  void OnCurrentSelectionChanged(QmitkSingleNodeSelectionWidget::NodeList nodes);//Pointset


  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  void OnEndInteraction();
  void OnAddBundle();


  Ui::QmitkInteractiveFiberDissectionViewControls* m_Controls;



  int m_IterationCounter;                                      ///< used for data node naming


  std::vector<mitk::DataNode::Pointer>  m_SelectedFB;       ///< selected fiber bundle nodes

  mitk::Image::Pointer                  m_SelectedImage;
  mitk::DataNode::Pointer               m_SelectedPS;
  mitk::DataNode::Pointer               m_SelectedImageNode;
  mitk::FiberBundle::Pointer            m_positiveBundle;
  mitk::FiberBundle::Pointer            m_negativeBundle;


  std::vector< mitk::DataNode::Pointer> m_positivSelectedBundles;
  std::vector< mitk::DataNode::Pointer> m_negativeSelectedBundles;
  vtkSmartPointer<vtkPolyData>          m_positiveFibersData;
  vtkSmartPointer<vtkPolyData>          m_negativeFibersData;

  vtkSmartPointer<vtkCellPicker>        m_picker1;


};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

