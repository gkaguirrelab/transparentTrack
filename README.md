# transparentTrack
Code to analyze pupil size and gaze location in IR video of the eye.

These MATLAB routines are designed to operate upon infra-red videos of the human eye and extract the elliptical boundary of the pupil and the location of the IR "glint" (first Purkinje image). Additional routines support calibration of absolute pupil size and gaze position, resulting in extracted time-series data that provide eye gaze in degrees of visual angle relative to a viewed screen and pupil size in mm^2.

Notably, this software is computationally intensive and is designed to be run post-hoc upon videos collected during an experimental session. A particular design goal is to provide an accurate fit to the pupil boundary when it is partially obscured by the eyelid. This circumstance is encountered when the pupil is large (as is seen in data collected under low-light conditions) or in people with retinal disease.

The central computation is to fit an ellipse to the pupil boundary in each image frame. This fitting operation is performed iteratively, aided by ever more informed constraints on the parameters that define the ellipse. The ellipse function is recast in "transparent" form, with parameters that define the ellipse by X center, Y center, area, eccentricity (e.g., aspect ratio), and theta (angle). Expressing the parameters of the ellipse in this way allows us to place linear and non-linear constraints upon different parameters. At a high level of description, the fitting approach involves:

- **Intensity segmentation to extract the boundary of the pupil**. A preliminary circle fit via the Hough transform and an adaptive size window is used.
- **Initial ellipse fit with minial parameter constraint**. As part of this first stage, the boundary of pupil is refined through the application of "cuts" of different angles or extent. The minimum cut that provides an acceptable ellipse fit is retained. This addresses obscuration of the pupil border by the eye lids or by non-uniform IR illumination of the pupil.
- **Estimation of secne geometry**. 




Notably, the non-linear search places boundaries upon the maximum allowed eccentricity, forcing even a partial arc of pupil boundary points to be fit with a relatively circular ellipse. Upper and lower boundary constraints are placed on the other parameters of the ellipse within the non-linear search, reflecting anatomic and device limits. In addition to identifying the parameters of the best fitting elipse, the boundary points are also subject to repeated division and then re-fit. The variation in the ellipse parameters that result from fitting these sub-divisions of the boundary points estimate the confidence with which the best fit parameters are known for this image.

Ellipse fitting to the initial and refined boundary is conducted within a constrained, non-linear search. 

An . 

After an initial ellipse fit to the pupil boundary in each frame of the video, the set of ellipses is used to estimate the geometry of the scene. Specifically, an estimate is made of the center of rotation of the eye, with the assumption that the eye is spherical and tht the image of the pupil is the orthogonal projection of the circular pupil to the image plane.



Following an initial ellipse fit, the parameters are subject to an empirical Bayes smoothing step. The initial ellipse fit (and the confidence in the parameters across sub-divisions) is treated as the likelihood. A temporal prior is constructed from a weighted, non-causal average of the values for each parameter. The relative influence of past and future timepoints is subject to weighting by a decaying exponential, with different integration times for the different ellipse parameters (e.g., we encourage the area of the pupil to change more gradually than its position).

There are many software parameters that control the behavior of the routines. While the default settings work well for some videos (including those that are part of the sandbox demo included in this repository), other parameter settings may be needed for videos with different qualities.

To install and configure transparentTrack, first install toolboxToolbox (tBtB), which provides for declarative dependency management for Matlab: https://github.com/ToolboxHub/ToolboxToolbox

Once tBtB is installed, transparentTrack (and all its dependencies) can be installed and readied for use with the command `tbUse('transparentTrack')`.
