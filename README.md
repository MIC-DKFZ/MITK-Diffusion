MITK Diffusion
==============

Copyright © German Cancer Research Center (DKFZ), [Division of Medical Image Computing (MIC)](https://www.dkfz.de/en/mic/index.php). Please make sure that your usage of this code is in compliance with the code [license](https://github.com/MIC-DKFZ/MITK-Diffusion/blob/master/LICENSE.txt).

[![DOI](https://zenodo.org/badge/195387051.svg)](https://zenodo.org/badge/latestdoi/195387051)

The MITK Diffusion application [[1,2]](#References) offers a selection of image analysis algorithms for the processing of diffusion-weighted MR images. It encompasses the research of the Division of Medical Image Computing at the German Cancer Research Center (DKFZ).

* [Downloads](#Downloads)
* [Requirements](#Requirements)
* [Features](#Features)
* [Related Links](#Related-Links)
* [Image Gallery](#Image-Gallery)
* [Building MITK Diffusion from source](#Building-MITK-Diffusion-from-source)
* [User Manual](http://docs.mitk.org/diffusion/nightly/)
* [Report a Bug](https://phabricator.mitk.org/maniphest/task/edit/form/29/)
* [References](#References)
* [Contact](#Contact)

## Downloads

Please have a look at the [requirements](#Requirements) for running MITK Diffusion with all its features successfully!

The nightly installers come as executable setup wizards that install MITK Diffusion on your system or alternatively as simple .tar.gz or .zip archive where you can execute MITK Diffusion and the command line apps "manually". Should there be no new installer for a while, please [contact](#Contact) us and report the issue.

**Windows 64-bit**:
* [Windows 10 executable (with Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_Windows-10_Python.exe.html)
* [Windows 10 zip (with Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_Windows-10_Python.zip.html)
* [Windows 10 executable (without Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_Windows-10_NoPython.exe.html)
* [Windows 10 zip (without Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_Windows-10_NoPython.zip.html)

**Linux 64-bit**:
* [Ubuntu 18.04 executable (with Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_ubuntu-18.04_Python.run.html)
* [Ubuntu 18.04 tarball (with Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_ubuntu-18.04_Python.tar.gz.html)
* [Ubuntu 18.04 executable (without Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_ubuntu-18.04_NoPython.run.html)
* [Ubuntu 18.04 tarball (without Python support)](http://www.mitk.org/download/diffusion/nightly/MITK-Diffusion_ubuntu-18.04_NoPython.tar.gz.html)


If you encounter any bugs, please report them in our [bugtracking](https://phabricator.mitk.org/maniphest/task/edit/form/29/) system or use the [MITK-users mailing list](http://mitk.org/wiki/MITK_Mailinglist). We are grateful for any feedback!


### Requirements

* Install Python 3.X: `sudo apt install python3 python3-pip` (Ubuntu) or from https://www.python.org/downloads/windows/ (Windows)
* Download our Python requirements file: [PythonRequirements.txt](https://github.com/MIC-DKFZ/MITK-Diffusion/tree/master/PythonRequirements.txt)
* Install the Python requirements: `pip3 install -r PythonRequirements.txt`
* If your are behind a proxy use `pip3 --proxy <proxy> install -r PythonRequirements.txt`

**For Windows users**:
MITK Diffusion requires the Microsoft Visual C++ 2017 Redistributable to be installed on the system. The MITK Diffusion installer automatically installs this redistributable for you if not already present on the system, but it needs administrative privileges to do so. So to install the redistributable, **run the MITK Diffusion installer as administrator**.

## Features

**Support for most established image formats**
* Images: DICOM, NIFTI, NRRD (peak and SH images compatible with MRtrix)
* Tractograms: fib/vtk, tck and trk.

**Image preprocessing**
* Registration
* Head-motion correction
* Denoising
* Skull stripping and brain mask segmentation (Linux only)
* Resampling, cropping, flipping and merging
* Header modifications
* Single volume extraction

**Diffusion gradient/b-value processing**
* b-value rounding
* Gradient direction flipping
* Gradient direction subsampling
* Averaging of gradient directions/volumes
* Gradient direction and b-value visualization

**ODF reconstruction and signal modelling**
* Tensor and Q-ball reconstruction
* Other reconstructions via Dipy wrapping (CSD, 3D SHORE, SFM) (Linux only)
* ODF peak calculation
* MRtrix or camino results can be imported

**Quantification of diffusion-weighted/tensor/ODF images**
* Intravoxel Incoherent Motion (IVIM) and diffusion kurtosis analysis
* Calculation of many other derived indices such as ADC, MD, GFA, FA, RA, AD, RD
* Image statistics

**Segmentation**
* Automatic white matter bundle segmentation (TractSeg) [[3]](#References) (Linux only)
* Automatic brain mask segmentation (Linux only)
* Manual image segmentation and operations on segmentations
* SOON: automatic brain tissue segmentation

**Fiber tractography**
* Global tractography [[4]](#References)
* Streamline tractography
    * Interactive (similar to [[5]](#References)) or seed image based
    * Deterministic or probabilistic
    * Peak, ODF, tensor and raw dMRI based. The latter one in conjunction with machine learning based tractography [[6]](#References)
    * Various possibilities for anatomical constraints.
    * Tractography priors in form of additional peak images, e.g. obtained using TractSeg

**Fiber processing**
* Tract dissection (parcellation or ROI based)
* Tract filtering by
    * length
    * curvature
    * direction
    * weight
    * density
* Tract resampling and compression
* Tract transformation
    * Mirroring
    * Rotating and translating
    * Registration (apply transform of previously performed image registration)
* Tract coloring
    * Curvature
    * Length
    * Weight
    * Scalar map (e.g. FA)
* Other operations
    * Join
    * Subtract
    * Copy
* Fiber clustering [[7]](#References)
* Fiber fitting and weighting similar to SIFT2 and LiFE [[8,9]](#References)
* Principal direction extraction (fibers --> peaks)
* Tract derived images:
    * Tract density images
    * Tract endpoint images
    * Tract envelopes

**Fiberfox dMRI simulations** [[10]](#References)
* Multi-compartment signal modeling
* Simulation of the k-space acquisition including
    * Compartment specific relaxation effects
    * Artifacts such as noise, spikes, ghosts, aliasing, distortions, signal drift, head motion, eddy currents and Gibbs ringing
    * Definition of important acquisition parameters such as bvalues and gradient directions, TE, TR, dwell time, partial Fourier, ...
* Manual definition of fiber configurations, e.g. for evaluation purposes
* Automatic generation of random fiber configurations

**Other features**
* Brain network statistics and visualization (connectomics)
* Interactive Python console (Linux only)
* Integrated screenshot maker
* Command line tools for most functionalities

## Related Links
* Great python package for logging your (MITK) command line experiments:
    * https://github.com/MIC-DKFZ/cmdint
    * `pip3 install cmdint`
* TractSeg reference data of 72 semiautomatically defined bundles in 105 HCP subjects: https://zenodo.org/record/1285152
* TractSeg python package: https://github.com/MIC-DKFZ/TractSeg
* Simulated dMRI images and ground truth of random fiber phantoms in various configurations: https://doi.org/10.5281/zenodo.2533250
* ISMRM 2015 Tractography Challenge Data: https://doi.org/10.5281/zenodo.572345 & https://doi.org/10.5281/zenodo.1007149


## Image Gallery

![](http://mitk.org/images/8/8f/MitkDiffusion.png) Screenshot of the MITK Diffusion Welcome Screen
<br/><br/><br/>

![](http://mitk.org/images/0/09/ScalarMaps.png) Scalar map visualization
<br/><br/><br/>

![](http://mitk.org/images/3/3b/Data_Tensors.png) Tensor Visualization
<br/><br/><br/>

![](http://mitk.org/images/5/5c/Data_ODF.png) ODF visualization
<br/><br/><br/>

![](http://mitk.org/images/7/73/Data_Peaks.png) Peak visualization (uniform white coloring)
<br/><br/><br/>

![](http://mitk.org/images/6/68/StreamlineTractography.png) Interactive tractography in MITK Diffusion. The tractogram updates automatically on parameter change and movement of the spherical seed region.
<br/><br/><br/>

![](http://mitk.org/images/3/3a/Extraction_1.png) Tract dissection using manually drawn ROIs.
<br/><br/><br/>

![](http://mitk.org/images/f/f5/FiberFit.png) Automatic streamline weighting (similar to SIFT2 or LiFE)
<br/><br/><br/>

![](http://mitk.org/images/3/33/Fiberfox.png) Illustration of the dMRI phantom simulation process using Fiberfox.
<br/><br/><br/>

<img src="http://mitk.org/images/e/ec/Eddy.gif" width="250" height="200" /><img src="http://mitk.org/images/a/ad/Spike_and_motion.gif" width="250" height="200" /><img src="http://mitk.org/images/8/8e/Drift.gif" width="250" height="200" />
<img src="http://mitk.org/images/0/08/645551_axial_motion_ghost_eddy.gif" width="250" height="200" /><img src="http://mitk.org/images/5/55/Ringing.gif" width="250" height="200" /><img src="http://mitk.org/images/3/37/683256_Distortions_Sagittal.gif" width="250" height="200" />\
Illustration of simulated dMRI images with various artifacts (a bit excessive for illustration purposes): eddy current distortions (1), motion and spike (2), intensity drift (3), motion, eddy and noise (4), ringing (5), B0 inhomogeneity distortions (6), from left to right.
<br/><br/><br/>

![](http://mitk.org/images/0/08/RandomFibers_Example.png)
Automatically generated random fiber configuration for Fiberfox simulations.
<br/><br/><br/>


## Building MITK Diffusion from source
* Install [Qt](https://www.qt.io/) on your system (>= 5.12.6).
* Clone MITK from [github](https://github.com/MIC-DKFZ/MITK-Diffusion.git) using [Git version control](https://git-scm.com/).
* Clone MITK Diffusion from [github](https://github.com/MITK/MITK.git).
* Configure the MITK Superbuild using [CMake](https://cmake.org/) (>= 3.14.5).
    * Choose the MITK source code directory and an empty binary directory.
    * Click "Configure".
    * Set the option MITK_EXTENSION_DIRS to "/path/to/my/mitk-diffusion-repository".
    * Click "Configure".
    * Set the option MITK_BUILD_CONFIGURATION to "DiffusionRelease".
    * Click "Generate".
    * macOS specifics (unclear if this still applies):
        * Use python 3.**6**, since python 3.**7** leads to build errors on macOS.
        * The cmake variables for python 3 might need to be set manually. It is probably enough to specify PYTHON_EXECUTABLE.
        * Openmp needs to be installed manually since it is not included in apple clang anymore: "brew install libomp" should do the trick. It might be necessary to set the corresponding make variables later in the MITK build manually:
            * OpenMP_CXX_FLAGS: -Xpreprocessor -fopenmp -I"/path/to/python3/includes/"
            * OpenMP_C_FLAGS: -Xpreprocessor -fopenmp -I"/path/to/python3/includes/"
            * OpenMPCXX_LIB_NAMES: libomp
            * OpenMPC_LIB_NAMES: libomp
            * OpenMP_libomp_LIBRARY: /path/to/libomp.dylib
* Start the Superbuild:
    * Linux/maxOS: Open a console window, navigate to the build folder and type "make -j8" (optionally supply the number threads to be used for a parallel build with -j).
    * Windows (requires visual studio): Open the MITK Superbuild solution file and build all projects.
    * The Superbuild may take some time.
    * After the Superbuild has finished, change the cmake binary directory from "/path/to/my/build" to "/path/to/my/build/MITK-build" and configure+generate again.
    * Build again in "/path/to/my/build/MITK-build" using make (linux/mac) or the VS solution file.
    * The build may again take some time and should yield the binaries in "/path/to/my/build/MITK-build/bin"

More detailed build instructions can be found in the [documentation](http://docs.mitk.org/nightly/BuildInstructionsPage.html).


Continuous integration: http://cdash.mitk.org/index.php?project=MITK&display=project


## References

All publications of the Division of Medical Image Computing can be found [https://www.dkfz.de/en/mic/publications/ here].

[1] Fritzsche, Klaus H., Peter F. Neher, Ignaz Reicht, Thomas van Bruggen, Caspar Goch, Marco Reisert, Marco Nolden, et al. “MITK Diffusion Imaging.” Methods of Information in Medicine 51, no. 5 (2012): 441.

[2] Fritzsche, K., and H.-P. Meinzer. “MITK-DI A New Diffusion Imaging Component for MITK.” In Bildverarbeitung Für Die Medizin, n.d.

[3] Wasserthal, Jakob, Peter Neher, and Klaus H. Maier-Hein. “TractSeg - Fast and Accurate White Matter Tract Segmentation.” NeuroImage 183 (August 4, 2018): 239–53.

[4] Neher, P. F., B. Stieltjes, M. Reisert, I. Reicht, H.P. Meinzer, and K. Maier-Hein. “MITK Global Tractography.” In SPIE Medical Imaging: Image Processing, 2012.

[5] Chamberland, M., K. Whittingstall, D. Fortin, D. Mathieu, und M. Descoteaux. „Real-time multi-peak tractography for instantaneous connectivity display“. Front Neuroinform 8 (2014): 59. doi:10.3389/fninf.2014.00059.

[6] Neher, Peter F., Marc-Alexandre Côté, Jean-Christophe Houde, Maxime Descoteaux, and Klaus H. Maier-Hein. “Fiber Tractography Using Machine Learning.” NeuroImage. Accessed July 17, 2017. doi:10.1016/j.neuroimage.2017.07.028.

[7] Garyfallidis, Eleftherios, Matthew Brett, Marta Morgado Correia, Guy B. Williams, and Ian Nimmo-Smith. “QuickBundles, a Method for Tractography Simplification.” Frontiers in Neuroscience 6 (2012).

[8] Smith, Robert E., Jacques-Donald Tournier, Fernando Calamante, and Alan Connelly. “SIFT2: Enabling Dense Quantitative Assessment of Brain White Matter Connectivity Using Streamlines Tractography.” NeuroImage 119, no. Supplement C (October 1, 2015): 338–51.

[9] Pestilli, Franco, Jason D. Yeatman, Ariel Rokem, Kendrick N. Kay, and Brian A. Wandell. “Evaluation and Statistical Inference for Human Connectomes.” Nature Methods 11, no. 10 (October 2014): 1058–63.

[10] Neher, Peter F., Frederik B. Laun, Bram Stieltjes, and Klaus H. Maier-Hein. “Fiberfox: Facilitating the Creation of Realistic White Matter Software Phantoms.” Magnetic Resonance in Medicine 72, no. 5 (November 2014): 1460–70. doi:10.1002/mrm.25045.


## Contact
If you have questions about the application or if you would like to give us feedback, don't hesitate to contact us using [our mailing list](http://mitk.org/wiki/MITK_Mailinglist) or, for questions that are of no interest for the community, [directly](https://www.dkfz.de/en/mic/team/people/Peter_Neher.html).
