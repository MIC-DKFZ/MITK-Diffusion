/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/



#include <berryISelectionListener.h>
#include <berryIStructuredSelection.h>

#include <QmitkAbstractView.h>
#include "ui_QmitkFiberfoxViewControls.h"
#include <itkVectorImage.h>
#include <itkVectorContainer.h>
#include <itkOrientationDistributionFunction.h>

#ifndef Q_MOC_RUN
#include <mitkFiberBundle.h>
#include <mitkPlanarEllipse.h>
#include <mitkDiffusionNoiseModel.h>
#include <mitkDiffusionSignalModel.h>
#include <mitkRicianNoiseModel.h>
#include <itkTractsToDWIImageFilter.h>
#include <mitkTensorModel.h>
#include <mitkBallModel.h>
#include <mitkStickModel.h>
#include <mitkAstroStickModel.h>
#include <mitkDotModel.h>
#include <mitkFiberfoxParameters.h>
#include <itkStatisticsImageFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#endif

#include <QThread>
#include <QObject>
#include <QTimer>
#include <QTime>

/*!
\brief View for fiber based diffusion software phantoms (Fiberfox). See "Fiberfox: Facilitating the creation of realistic white matter software phantoms" (DOI: 10.1002/mrm.25045) for details.
*/

// Forward Qt class declarations

class QmitkFiberfoxView;

class QmitkFiberfoxWorker : public QObject
{
    Q_OBJECT

public:

    QmitkFiberfoxWorker(QmitkFiberfoxView* view);

public slots:

    void run();

private:

    QmitkFiberfoxView*                  m_View;
};

class QmitkFiberfoxView : public QmitkAbstractView
{

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

    QmitkFiberfoxView();
    virtual ~QmitkFiberfoxView();

    virtual void CreateQtPartControl(QWidget *parent) override;
    void SetFocus() override;

    typedef mitk::DiffusionPropertyHelper::GradientDirectionType            GradientDirectionType;
    typedef mitk::DiffusionPropertyHelper::GradientDirectionsContainerType  GradientDirectionContainerType;
    typedef itk::Vector<double,3>           GradientType;
    typedef std::vector<GradientType>       GradientListType;
    typedef itk::VectorImage< short, 3 >    ItkDwiType;
    typedef itk::Image<double, 3>           ItkDoubleImgType;
    typedef itk::Image<float, 3>            ItkFloatImgType;
    typedef itk::Image<unsigned char, 3>    ItkUcharImgType;

protected slots:

    void SetOutputPath();           ///< path where image is automatically saved to after the simulation is finished
    void LoadParameters();          ///< load fiberfox parameters
    void SaveParameters();          ///< save fiberfox parameters
    void SetBvalsEdit();
    void SetBvecsEdit();

    void BeforeThread();
    void AfterThread();
    void KillThread();              ///< abort simulation
    void UpdateSimulationStatus();  ///< print simulation progress and satus messages

    void GenerateImage();           ///< start image simulation
    void Comp1ModelFrameVisibility(int index);  ///< only show parameters of selected signal model for compartment 1
    void Comp2ModelFrameVisibility(int index);  ///< only show parameters of selected signal model for compartment 2
    void Comp3ModelFrameVisibility(int index);  ///< only show parameters of selected signal model for compartment 3
    void Comp4ModelFrameVisibility(int index);  ///< only show parameters of selected signal model for compartment 4
    void ShowAdvancedOptions(int state);

    /** update GUI elements */
    void OnAddNoise(int value);
    void OnAddGhosts(int value);
    void OnAddDistortions(int value);
    void OnAddEddy(int value);
    void OnAddSpikes(int value);
    void OnAddAliasing(int value);
    void OnAddMotion(int value);
    void OnAddDrift(int value);
    void OnAddRinging(int value);
    void OnMaskSelected(int value);
    void OnFibSelected(int value);
    void OnTemplateSelected(int value);
    void OnBvalsBvecsCheck(int value);
    void OnTlineChanged();
    void OnGradientSettingsChanged();   ///< loaded gradients and prototypes give way to the changed settings

protected:

    Ui::QmitkFiberfoxViewControls* m_Controls;

    void StartSimulation();     ///< GenerateImage without its exception handling
    void SimulateForExistingDwi(mitk::DataNode* imageNode);     ///< add artifacts to existing diffusion weighted image
    void SimulateImageFromFibers(mitk::DataNode* fiberNode);    ///< simulate new diffusion weighted image
    bool ConfirmMemory();        ///< asks before a run whose estimated memory exceeds the available memory; false if cancelled
    void UpdateParametersFromGui();  ///< update fiberfox paramater object: the Read* functions in this order, then SetResultProperties
    void ReadOutputSettings();       ///< output path and additional images
    void ReadImageGeometry();        ///< mask, image geometry (template, mask or boxes) and gradients
    void ReadArtifacts();            ///< relaxation, k-space artifacts and motion
    void ReadSequence();             ///< acquisition type, coils and timing
    void ReadNoise();                ///< noise distribution and variance
    void ReadCompartments();         ///< signal model and volume fraction map of each compartment
    void SetResultProperties();      ///< "Fiberfox.*" properties of the result node, from the parameters read
    void SetGradientsFromGui();      ///< bvals/bvecs files, loaded gradients or a generated half shell
    /** Shows SignalGenerationParameters::m_KspaceTimeSegments: exact (0), automatic (negative) or a fixed count. */
    void ShowKspaceTimeSegments(int segments);
    void EnableKspacePhase();        ///< the k-space phase option applies only with distortions or eddy currents

    /** Signal model types the compartment boxes offer. */
    enum class CompartmentModel { None, Stick, Zeppelin, Tensor, Ball, Astrosticks, Dot, Prototype };
    /** The controls of one compartment: the model type and parameter widget of each entry of its box. */
    struct CompartmentControls
    {
      int id;   ///< m_CompartmentId of its signal model, 1 to 4
      std::string description;
      QComboBox* box;
      QmitkDataStorageComboBoxWithSelectNone* volumeFraction;
      QWidget* fractionFrame;   ///< shown only with a model selected; nullptr where always shown
      std::vector< std::pair<CompartmentModel, QWidget*> > entries;   ///< by box index; no widget for None
    };
    /**
     * The one mapping between box entries, parameter widgets and signal models, used to create the
     * models, to show loaded ones and to show the selected widget. Compartment id at index id-1.
     */
    std::vector<CompartmentControls> Compartments() const;
    static std::string ModelName(CompartmentModel type);   ///< as in result names and properties
    static CompartmentModel SelectedModel(const CompartmentControls& compartment);
    static QWidget* ModelWidget(const CompartmentControls& compartment, CompartmentModel type);   ///< nullptr if not offered
    /** Signal model of the selected entry with the widget's parameters, nullptr for none. */
    mitk::DiffusionSignalModel<>* NewCompartmentModel(const CompartmentControls& compartment);
    /**
     * Selects the entry of the model's type and shows the model's parameters; nothing if not offered.
     * A tensor model is shown as a zeppelin if zeppelin is set and the compartment offers one.
     */
    void ShowCompartmentModel(const CompartmentControls& compartment, mitk::DiffusionSignalModel<>* model, bool zeppelin);
    /** Shows the parameter widget of the entry with the given index only. */
    void ShowCompartmentWidget(int compartmentId, int index);
    void UpdateGui();                                           ///< enable/disbale buttons etc. according to current datamanager selection
    void PlanarFigureSelected( itk::Object* object, const itk::EventObject& );
    void EnableCrosshairNavigation();               ///< enable crosshair navigation if planar figure interaction ends
    void DisableCrosshairNavigation();              ///< disable crosshair navigation if planar figure interaction starts

    void SaveParameters(QString filename);

    /** The widget of a compartment's prototype model, nullptr for compartment 2 (none). */
    QmitkPrototypeSignalParametersWidget* GetPrototypeWidget(int compartmentId);
    /** A prototype model of the widget's settings, with the loaded kernels if they still apply. */
    mitk::RawShModel<>* NewPrototypeModel(QmitkPrototypeSignalParametersWidget* widget, int compartmentId);
    /** Prototype signal models of m_Parameters that have no kernels yet. */
    std::vector< mitk::RawShModel<>* > GetPrototypeModelsWithoutKernels();
    /** Samples the kernels of the given models from the DWI; false if a model got no valid kernel. */
    bool SamplePrototypeKernels(mitk::Image::Pointer diffImg, const std::vector< mitk::RawShModel<>* >& models);


    mitk::DataNode::Pointer                             m_SelectedImageNode;

    QString m_ParameterFile;    ///< parameter file name

    // GUI thread
    QmitkFiberfoxWorker     m_Worker;   ///< runs filter
    QThread                 m_Thread;   ///< worker thread
    bool                    m_ThreadIsRunning;
    QTimer*                 m_SimulationTimer;
    QTime                   m_SimulationTime;
    QString                 m_SimulationStatusText;
    std::string             m_SimulationError;  ///< message of the exception that stopped the last simulation, empty if none

    /** Image filters that do all the simulations. */
    itk::TractsToDWIImageFilter< short >::Pointer           m_TractsToDwiFilter;

    friend class QmitkFiberfoxWorker;
    mitk::FiberfoxParameters  m_Parameters;
    mitk::SignalGenerationParameters::GradientListType  m_LoadedGradients;  ///< of the loaded .ffp, until a gradient setting changes
    double                    m_LoadedBvalue = 0;
    /** Prototype model with kernels from a loaded .ffp and the widget settings it was loaded with. */
    struct LoadedPrototype
    {
      std::unique_ptr< mitk::DiffusionSignalModel<> > m_Model;
      std::vector<double>       m_WidgetSettings;
      std::vector<int>          m_Shells;   ///< non-zero b-values of the file's gradients
    };
    std::map< int, LoadedPrototype > m_LoadedPrototypes;   ///< by compartment id
    /** m_Parameters of the last started simulation before ApplyDirectionMatrix (gradients as FSL
     * bvecs), saved as .ffp next to the simulated image */
    mitk::FiberfoxParameters  m_SimulatedParameters;
};
