/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_FiberfoxParameters_H
#define _MITK_FiberfoxParameters_H

#include <itkImageRegion.h>
#include <itkMatrix.h>
#include <mitkDiffusionNoiseModel.h>
#include <mitkDiffusionSignalModel.h>
#include <mitkDataNode.h>
#include <mitkRicianNoiseModel.h>
#include <mitkNoncentralChiNoiseModel.h>
#include <mitkImage.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkPlanarEllipse.h>
#include <mitkStickModel.h>
#include <mitkTensorModel.h>
#include <mitkAstroStickModel.h>
#include <mitkBallModel.h>
#include <mitkDotModel.h>
#include <mitkRawShModel.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <limits>
#include <MitkMriSimulationExports.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace mitk
{
  class MITKMRISIMULATION_EXPORT FiberfoxParameters;

  /** Signal generation */
  class MITKMRISIMULATION_EXPORT SignalGenerationParameters
  {
    friend FiberfoxParameters;
  public:
    typedef itk::Image<float, 3>                    ItkFloatImgType;
    typedef itk::Image<unsigned char, 3>            ItkUcharImgType;
    typedef itk::Vector<double,3>                   GradientType;
    typedef std::vector<GradientType>               GradientListType;

    enum CoilSensitivityProfile : int {
      COIL_CONSTANT,
      COIL_LINEAR,
      COIL_EXPONENTIAL
    };

    enum AcquisitionType : int
    {
      SingleShotEpi,
      ConventionalSpinEcho,
      FastSpinEcho,
    };

    SignalGenerationParameters()
      : m_AcquisitionType(SignalGenerationParameters::SingleShotEpi)
      , m_SignalScale(100)
      , m_tEcho(100)
      , m_tRep(4000)
      , m_tInv(0)
      , m_tLine(1)
      , m_tInhom(50)
      , m_EchoTrainLength(8)
      , m_ReversePhase(false)
      , m_PartialFourier(1.0)
      , m_NoiseVariance(0.001f)
      , m_NumberOfCoils(1)
      , m_CoilSensitivityProfile(SignalGenerationParameters::COIL_CONSTANT)
      , m_CoilSensitivity(0.3f)
      , m_SimulateKspaceAcquisition(false)
      , m_AxonRadius(0)
      , m_DoDisablePartialVolume(false)
      , m_Spikes(0)
      , m_SpikeAmplitude(1)
      , m_KspaceLineOffset(0)
      , m_EddyStrength(0.002f)
      , m_Tau(70)
      , m_CroppingFactor(1)
      , m_Drift(0.06)
      , m_DoAddGibbsRinging(false)
      , m_ZeroRinging(0)
      , m_KspaceTimeSegments(0)
      , m_DoSimulateRelaxation(true)
      , m_DoAddMotion(false)
      , m_DoRandomizeMotion(true)
      , m_DoAddDrift(false)
      , m_FrequencyMap(nullptr)
      , m_MaskImage(nullptr)
      , m_Bvalue(1000)
    {
      m_ImageRegion.SetSize(0, 12);
      m_ImageRegion.SetSize(1, 12);
      m_ImageRegion.SetSize(2, 3);
      m_ImageSpacing.Fill(2.0);
      m_ImageOrigin.Fill(0.0);
      m_ImageDirection.SetIdentity();
      m_Translation.Fill(0.0);
      m_Rotation.Fill(0.0);
      SetNumWeightedVolumes(6);
    }

    /** input/output image specifications */
    itk::ImageRegion<3>                 m_CroppedRegion;            ///< Image size with reduced FOV.
    itk::ImageRegion<3>                 m_ImageRegion;              ///< Image size.
    itk::Vector<float,3>                m_ImageSpacing;             ///< Image voxel size.
    itk::Point<float,3>                 m_ImageOrigin;              ///< Image origin.
    itk::Matrix<double, 3, 3>           m_ImageDirection;           ///< Image rotation matrix.

    /** Other acquisitions parameters */
    AcquisitionType                     m_AcquisitionType;          ///< determines k-space trajectory and maximum echo position(s)
    float                               m_SignalScale;              ///< Scaling factor for output signal (before noise is added).
    float                               m_tEcho;                    ///< Echo time TE.
    float                               m_tRep;                     ///< Echo time TR.
    float                               m_tInv;                       ///< Inversion time
    float                               m_tLine;                    ///< k-space line readout time (dwell time).
    float                               m_tInhom;                   ///< T2'
    unsigned int                        m_EchoTrainLength;          ///< Only relevant for Fast Spin Echo sequence (number of k-space lines acquired with one RF pulse)
    bool                                m_ReversePhase;             ///< If true, the phase readout direction will be inverted (-y instead of y)
    float                               m_PartialFourier;           ///< Partial fourier factor (0.5-1)
    float                               m_NoiseVariance;            ///< Variance of the complex Gaussian k-space noise; not added when FiberfoxParameters::m_NoiseModel is set
    unsigned int                        m_NumberOfCoils;            ///< Number of coils in multi-coil acquisition
    CoilSensitivityProfile              m_CoilSensitivityProfile;   ///< Choose between constant, linear or exponential sensitivity profile of the used coils
    float                               m_CoilSensitivity;          ///< signal remaining in slice center
    bool                                m_SimulateKspaceAcquisition;///< Flag to enable/disable k-space acquisition simulation
    double                              m_AxonRadius;               ///< Determines compartment volume fractions (0 == automatic axon radius estimation)
    bool                                m_DoDisablePartialVolume;   ///< Disable partial volume effects. Each voxel is either all fiber or all non-fiber.

    /** Artifacts and other effects */
    unsigned int                        m_Spikes;                   ///< Number of spikes randomly appearing in the image
    float                               m_SpikeAmplitude;           ///< amplitude of spikes relative to the largest signal intensity (magnitude of complex)
    float                               m_KspaceLineOffset;         ///< Causes N/2 ghosts. Larger offset means stronger ghost.
    float                               m_EddyStrength;             ///< Strength of eddy current induced gradients in mT/m.
    float                               m_Tau;                      ///< Eddy current decay constant (in ms)
    float                               m_CroppingFactor;           ///< FOV size in y-direction is multiplied by this factor. Causes aliasing artifacts.
    float                               m_Drift;                    ///< Global signal decrease by the end of the acquisition.
    bool                                m_DoAddGibbsRinging;        ///< Add Gibbs ringing artifact
    int                                 m_ZeroRinging;              ///< If > 0, ringing is simulated by by setting the defined percentage of higher frequencies to 0 in k-space. Otherwise, the input to the k-space simulation is generated with twice the resolution and cropped during k-space simulation (much slower).
    /**
     * How the k-space simulation sums the phase of distortions (field map) and eddy currents: 0, the
     * default and the value for parameter files without the setting, sums the exact phase of every
     * pixel for every k-space sample; a number L > 0 approximates the phase by L time segments, a few
     * transforms per slice instead of the full sum (itk::KspaceImageFilter::ComputeTimeSegmentation);
     * AUTOMATIC_TIME_SEGMENTS chooses L per slice from the phase spread, with the error bound stated there.
     */
    int                                 m_KspaceTimeSegments;
    static constexpr int                AUTOMATIC_TIME_SEGMENTS = -1;
    bool                                m_DoSimulateRelaxation;     ///< Add T2 relaxation effects
    bool                                m_DoAddMotion;              ///< Enable motion artifacts.
    bool                                m_DoRandomizeMotion;        ///< Toggles between random and linear motion.
    bool                                m_DoAddDrift;               ///< Add quadratic signal drift.
    std::vector< bool >                 m_MotionVolumes;            ///< Indicates the image volumes that are affected by motion; empty: all, missing entries: not moved
    ///< with positive numbers, inverted logic with negative numbers.
    itk::Vector<float,3>                m_Translation;              ///< Maximum translational motion.
    itk::Vector<float,3>                m_Rotation;                 ///< Maximum rotational motion.
    ItkFloatImgType::Pointer            m_FrequencyMap;             ///< If != nullptr, distortions are added to the image using this frequency map.
    ItkUcharImgType::Pointer            m_MaskImage;                ///< Signal is only genrated inside of the mask image.

    std::vector< int > GetBaselineIndices();                 ///< Returns list of nun-diffusion-weighted image volume indices
    unsigned int GetFirstBaselineIndex();                    ///< Returns index of first non-diffusion-weighted image volume
    bool IsBaselineIndex(unsigned int idx);                  ///< Checks if image volume with given index is non-diffusion-weighted volume or not.
    unsigned int GetNumWeightedVolumes();                    ///< Get number of diffusion-weighted image volumes
    unsigned int GetNumBaselineVolumes();                    ///< Get number of non-diffusion-weighted image volumes
    unsigned int GetNumVolumes();                            ///< Get number of baseline and diffusion-weighted image volumes
    GradientListType GetGradientDirections();                ///< Return gradient direction container
    mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer GetItkGradientContainer();
    GradientType GetGradientDirection(unsigned int i);
    std::vector< int > GetBvalues();                         ///< Returns a vector with all unique b-values (determined by the gradient magnitudes)
    double GetBvalue();
    /**
     * The gradients of .ffp and bvecs files are FSL bvecs of the simulated geometry (voxel frame, x
     * negated for a direction matrix with positive determinant), like every bvecs file MITK reads.
     * This turns them into the world (LPS) vectors the simulation uses (and stores as the output's
     * gradients). Without it, the gradients are taken as world vectors.
     */
    void ApplyDirectionMatrix();
    /**
     * Sets m_MotionVolumes, one flag per volume, from the motion-volume text of the Fiberfox view
     * or of a .ffp file (leading and trailing spaces ignored): "random" moves every volume but the
     * first with probability 1/2, a blank text moves every volume, positive numbers name the
     * moved volumes and negative numbers the volumes that stay. Any other text (no numbers, both
     * signs) moves no volume and returns false.
     */
    bool SetMotionVolumes(const std::string& text, itk::Statistics::MersenneTwisterRandomVariateGenerator* randgen);
    /**
     * Throws a mitk::Exception if relaxation is simulated, an inversion pulse is set (TI > 0) and
     * TI is not shorter than TR: the inversion comes before the excitation within one repetition,
     * and the inversion-recovery signal 1 - 2 exp(-TI/T1) + exp(-TR/T1) holds only for TI < TR.
     * Without relaxation TI is not used, so any TI passes.
     */
    void CheckInversionTime() const;

  protected:

    unsigned int                        m_NumGradients;             ///< Number of diffusion-weighted image volumes.
    unsigned int                        m_NumBaseline;              ///< Number of non-diffusion-weighted image volumes.
    GradientListType                    m_GradientDirections;       ///< Total number of image volumes.
    double                              m_Bvalue;                   ///< Acquisition b-value

    void SetNumWeightedVolumes(int numGradients);            ///< Automaticall calls GenerateGradientHalfShell() afterwards.
    void SetGradienDirections(GradientListType gradientList);
    void SetGradienDirections(mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer gradientList);
    void GenerateGradientHalfShell();                        ///< Generates half shell of gradient directions (with m_NumGradients non-zero directions)
  };

  /** Fiber generation */
  class MITKMRISIMULATION_EXPORT FiberGenerationParameters
  {
  public:

    enum FiberDistribution
    {
      DISTRIBUTE_UNIFORM, // distribute fibers uniformly in the ROIs
      DISTRIBUTE_GAUSSIAN // distribute fibers using a 2D gaussian
    };

    typedef std::vector< std::vector< mitk::PlanarEllipse::Pointer > >    FiducialListType;
    typedef std::vector< std::vector< unsigned int > >                    FlipListType;

    FiberGenerationParameters()
      : m_Distribution(DISTRIBUTE_UNIFORM)
      , m_Density(100)
      , m_Variance(100)
      , m_Sampling(1)
      , m_Tension(0)
      , m_Continuity(0)
      , m_Bias(0)
    {
      m_Rotation.Fill(0.0);
      m_Translation.Fill(0.0);
      m_Scale.Fill(1.0);
    }

    FiberDistribution   m_Distribution;
    unsigned int        m_Density;
    double              m_Variance;
    double              m_Sampling;
    double              m_Tension;
    double              m_Continuity;
    double              m_Bias;
    mitk::Vector3D      m_Rotation;
    mitk::Vector3D      m_Translation;
    mitk::Vector3D      m_Scale;
    FlipListType        m_FlipList;        ///< contains flags indicating a flip of the 2D fiber x-coordinates (needed to resolve some unwanted fiber twisting)
    FiducialListType    m_Fiducials;       ///< container of the planar ellipses used as fiducials for the fiber generation process
  };

  /** GUI persistence, input, output, ... */
  class MITKMRISIMULATION_EXPORT MiscFiberfoxParameters
  {
  public:
    MiscFiberfoxParameters()
      : m_ResultNode(DataNode::New())
      , m_ParentNode(nullptr)
      , m_SignalModelString("")
      , m_ArtifactModelString("")
      , m_OutputPath("/tmp/")
      , m_OutputPrefix("fiberfox")
      , m_AfterSimulationMessage("")
      , m_BvalsFile("")
      , m_BvecsFile("")
      , m_OutputAdditionalImages(false)
      , m_CheckAdvancedSignalOptionsBox(false)
      , m_DoAddNoise(false)
      , m_DoAddGhosts(false)
      , m_DoAddAliasing(false)
      , m_DoAddSpikes(false)
      , m_DoAddEddyCurrents(false)
      , m_DoAddDistortions(false)
      , m_MotionVolumesBox("random")
      , m_CheckRealTimeFibersBox(true)
      , m_CheckAdvancedFiberOptionsBox(false)
      , m_CheckConstantRadiusBox(false)
      , m_CheckIncludeFiducialsBox(true)
    {}

    DataNode::Pointer   m_ResultNode;                       ///< Stores resulting image.
    DataNode::Pointer   m_ParentNode;                       ///< Parent node of result node.
    std::string         m_SignalModelString;                ///< Appendet to the name of the result node
    std::string         m_ArtifactModelString;              ///< Appendet to the name of the result node
    std::string         m_OutputPath;                       ///< Image is automatically saved to the specified folder after simulation is finished.
    std::string         m_OutputPrefix;  /** Prefix for filename of output files and logfile. */
    std::string         m_AfterSimulationMessage;           ///< Store messages that are displayed after the simulation has finished (e.g. warnings, automatic parameter adjustments etc.)
    std::string         m_BvalsFile;
    std::string         m_BvecsFile;

    /** member variables that store the check-state of GUI checkboxes */
    // image generation
    bool                m_OutputAdditionalImages;
    bool                m_CheckAdvancedSignalOptionsBox;
    bool                m_DoAddNoise;
    bool                m_DoAddGhosts;
    bool                m_DoAddAliasing;
    bool                m_DoAddSpikes;
    bool                m_DoAddEddyCurrents;
    bool                m_DoAddDistortions;
    std::string         m_MotionVolumesBox;
    // fiber generation
    bool                m_CheckRealTimeFibersBox;
    bool                m_CheckAdvancedFiberOptionsBox;
    bool                m_CheckConstantRadiusBox;
    bool                m_CheckIncludeFiducialsBox;
  };

  /**
  * \brief Datastructure to manage the Fiberfox signal generation parameters.
  *
  */
  class MITKMRISIMULATION_EXPORT FiberfoxParameters
  {

  public:


    typedef itk::Image<float, 3>                            ItkFloatImgType;
    typedef itk::Image<double, 3>                           ItkDoubleImgType;
    typedef itk::Image<unsigned char, 3>                    ItkUcharImgType;
    typedef DiffusionSignalModel<double>                    DiffusionModelType;
    typedef std::vector< DiffusionModelType* >              DiffusionModelListType;
    typedef DiffusionNoiseModel<double>                     NoiseModelType;

    FiberfoxParameters();
    /** Copy and assignment clone the signal models and the noise model. */
    FiberfoxParameters(const FiberfoxParameters &params);
    FiberfoxParameters& operator=(const FiberfoxParameters &params);
    ~FiberfoxParameters();

    /** Not templated parameters */
    FiberGenerationParameters           m_FiberGen;             ///< Fiber generation parameters
    SignalGenerationParameters          m_SignalGen;            ///< Signal generation parameters
    MiscFiberfoxParameters              m_Misc;                 ///< GUI realted and I/O parameters

    /** Templated parameters. The parameters own the models in both lists (created with new) and delete them. */
    DiffusionModelListType              m_FiberModelList;       ///< Intra- and inter-axonal compartments.
    DiffusionModelListType              m_NonFiberModelList;    ///< Extra-axonal compartments.
    std::shared_ptr< NoiseModelType >   m_NoiseModel;           ///< If != nullptr, noise is added to the image.

    void GenerateGradientHalfShell();
    void SetNumWeightedVolumes(int numGradients);            ///< Automaticall calls GenerateGradientHalfShell() afterwards.
    void SetGradienDirections(mitk::SignalGenerationParameters::GradientListType gradientList);
    void SetGradienDirections(mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::ConstPointer gradientList);
    void SetBvalue(double Bvalue);
    void UpdateSignalModels();
    void ClearFiberParameters();
    void ClearSignalParameters();
    void ApplyDirectionMatrix();
    /**
     * Gradients and reference b-value of a diffusion-weighted image, for a simulation that reproduces
     * its world gradients; m_SignalGen.m_ImageDirection has to be set to the simulated geometry first.
     * Stored as FSL bvecs of that geometry (the form of .ffp files, see ApplyDirectionMatrix) or,
     * with asFslBvecs false, as world vectors for a simulation without ApplyDirectionMatrix.
     */
    void SetGradientsFromImage(const mitk::Image* dwi, bool asFslBvecs=true);
    /**
     * For the .ffp saved next to a simulated image: with motion, the motion-volume text becomes the
     * list of the volumes that m_SignalGen.m_MotionVolumes moves, so that loading the file moves the
     * same volumes ("random" would draw anew). An empty list stays blank (every volume moves, in
     * SetMotionVolumes and TractsToDWIImageFilter alike); a list without a moved volume switches
     * motion off, because a blank text would move every volume.
     */
    void ListMovedVolumes();

    void PrintSelf();                           ///< Print parameters to stdout.
    void SaveParameters(std::string filename);  ///< Save image generation parameters to .ffp file.
    void LoadParameters(std::string filename, bool fix_seed=false);  ///< Load image generation parameters from .ffp file.
    template< class ParameterType >
    ParameterType ReadVal(boost::property_tree::ptree::value_type const& v, std::string tag, ParameterType defaultValue, bool essential=false);
    std::string                         m_MissingTags;

  private:
    void DeleteSignalModels();
    void CopyFrom(const FiberfoxParameters &params);
  };
}

#endif

