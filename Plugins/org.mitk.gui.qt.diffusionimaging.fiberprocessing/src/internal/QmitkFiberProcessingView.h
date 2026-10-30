/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef QmitkFiberProcessingView_h
#define QmitkFiberProcessingView_h

#include <QmitkAbstractView.h>
#include "ui_QmitkFiberProcessingViewControls.h"

#include <mitkPlanarFigureComposite.h>
#include <mitkFiberBundle.h>
#include <mitkImage.h>
#include <mitkSurface.h>

#include <vtkSmartPointer.h>

/*!
\brief View to process fiber bundles. Supplies methods to extract fibers from the bundle, fiber resampling, mirroring, join and subtract bundles and much more.
*/
class QmitkFiberProcessingView : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  typedef itk::Image< unsigned char, 3 >    ItkUCharImageType;
  typedef itk::Image< float, 3 >            ItkFloatImageType;

  static const std::string VIEW_ID;

  QmitkFiberProcessingView();
  virtual ~QmitkFiberProcessingView();

  virtual void CreateQtPartControl(QWidget *parent) override;

  ///
  /// Sets the focus to an internal widget.
  ///
  virtual void SetFocus() override;

protected slots:

  void OnDrawCircle();          ///< add circle interactors etc.
  void OnDrawPolygon();         ///< add circle interactors etc.
  void GenerateAndComposite();
  void GenerateOrComposite();
  void GenerateNotComposite();

  void CopyBundles();               ///< add copies of selected bundles to data storage
  void JoinBundles();               ///< merge selected fiber bundles
  void SubtractBundles();          ///< subtract bundle A from bundle B. Not commutative! Defined by order of selection.
  void GenerateRoiImage();          ///< generate binary image of selected planar figures.

  void Remove();
  void Extract();
  void Modify();
  void UpdateGui();     ///< update button activity etc. dpending on current datamanager selection
  void OnMaskExtractionChanged();

  virtual void AddFigureToDataStorage(mitk::PlanarFigure* figure, const QString& name, const char *propertyKey = nullptr, mitk::BaseProperty *property = nullptr );

protected:

  void MirrorFibers();              ///< mirror bundle on the specified plane
  void ResampleSelectedBundlesSpline();   ///<
  void ResampleSelectedBundlesLinear();   ///<
  void DoImageColorCoding();        ///< color fibers by selected scalar image
  void DoWeightColorCoding();       ///< color fibers by their respective weights
  void DoLengthColorCoding();       ///< color fibers by length
  void DoCurvatureColorCoding();    ///< color fibers by curvature
  void CompressSelectedBundles();   ///< remove points below certain error threshold
  void WeightFibers();
  void ApplyWeightThreshold();
  void ApplyDensityThreshold();

  void RemoveWithMask(bool removeInside);
  void RemoveDir();
  void ApplyCurvatureThreshold();   ///< remove/split fibers with a too high curvature threshold
  void PruneBundle();               ///< remove too short/too long fibers

  void ExtractWithMask(bool onlyEnds, bool invert, bool labelmap);
  void ExtractWithPlanarFigure(bool interactive=false);

  void OnEndInteraction();

  /// \brief called by QmitkAbstractView when DataManager's selection has changed
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes) override;

  Ui::QmitkFiberProcessingViewControls* m_Controls;

  int m_CircleCounter;                                      ///< used for data node naming
  int m_PolygonCounter;                                     ///< used for data node naming
  std::vector<mitk::DataNode::Pointer>  m_SelectedFB;       ///< selected fiber bundle nodes
  std::vector<mitk::DataNode::Pointer>  m_SelectedPF;       ///< selected planar figure nodes
  std::vector<mitk::Surface::Pointer>   m_SelectedSurfaces;
  mitk::Image::Pointer                  m_SelectedImage;
  float                                 m_UpsamplingFactor; ///< upsampling factor for all image generations
  mitk::DataNode::Pointer               m_RoiImageNode;

  /** The planar figures this view observes, with the tag of its observer */
  std::vector<std::pair<mitk::PlanarFigure::Pointer, unsigned long>> m_FigureObservers;
  mitk::DataNode::Pointer               m_InteractiveNode;

  void AddCompositeToDatastorage(mitk::DataNode::Pointer pfc, std::vector<mitk::DataNode::Pointer> children, mitk::DataNode::Pointer parentNode=nullptr);
  void debugPFComposition(mitk::PlanarFigureComposite::Pointer , int );
  /** Adds the planar figures of the node (a figure or an OR composite) to the ROI image */
  void WritePfToImage(mitk::DataNode::Pointer node, ItkUCharImageType& image);
  mitk::DataNode::Pointer GenerateTractDensityImage(mitk::FiberBundle::Pointer fib, bool binary, bool absolute);
  mitk::DataNode::Pointer GenerateColorHeatmap(mitk::FiberBundle::Pointer fib);
  mitk::DataNode::Pointer GenerateFiberEndingsImage(mitk::FiberBundle::Pointer fib);
  mitk::DataNode::Pointer GenerateFiberEndingsPointSet(mitk::FiberBundle::Pointer fib);

  /** Fills an empty bundle or planar figure selection with the topmost top-level node of that
   *  kind, other than excluded. */
  void SelectTopmostIfEmpty(const mitk::DataNode* excluded);

  void NodeAdded( const mitk::DataNode* node ) override;
  void NodeRemoved(const mitk::DataNode* node) override;
  void RemoveObservers();
  void AddObservers();
};



#endif // _QMITKFIBERTRACKINGVIEW_H_INCLUDED

