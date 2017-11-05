# transparentTrack
Code to analyze pupil size and gaze location in IR video of the eye.

These MATLAB routines are designed to operate upon infra-red videos of the human eye and extract the elliptical boundary of the pupil and the location of the IR "glint" (first Purkinje image). Additional routines support calibration of absolute pupil size and gaze position, resulting in extracted time-series data that provide eye gaze in degrees of visual angle relative to a viewed screen and pupil size in mm^2.

Notably, this software is computationally intensive and slow, and thus designed to be run post-hoc upon videos collected during an eperimental session. A particular design goal is to provide good fitting of the pupil when it is partially obscured by the eyelid. This is particularly the case when the pupil is large, as is seen in data collected under low-light conditions or in people with retinal disease.

An initial intensity segmentation is used to extract the boundary of the pupil, assisted by a preliminary circle fit via the Hough transform. We observe that, when the pupil is large, this boundary is frequently impacted by the upper eyelid and by non-uniformity in the IR illumination of the pupil. To address this, the boundary is subject to a series of exploratory "cuts", in which boundary points are removed and the remaining points fit with an ellipse. The minimum cut that provides an acceptable ellipse fit is retained.

Ellipse fitting to the initial and refined boundary is conducted within a constrained, non-linear search. The parameters of the ellipse function are recast in "transparent" form, which define the ellipse by X center, Y center, area, eccentricity (e.g., aspect ratio), and theta (angle). Notably, the non-linear search places boundaries upon the maximum allowed eccentricity, forcing even a partial arc of pupil boundary points to be fit with a relatively circular ellipse. Upper and lower boundary constraints are placed on the other parameters of the ellipse within the non-linear search, reflecting anatomic and device limits. In addition to identifying the parameters of the best fitting elipse, the boundary points are also subject to repeated division and then re-fit. The variation in the ellipse parameters that result from fitting these sub-divisions of the boundary points estimate the confidence with which the best fit parameters are known for this image.

Following an initial ellipse fit, the parameters are subject to an empirical Bayes smoothing step. The initial ellipse fit (and the confidence in the parameters across sub-divisions) is treated as the likelihood. A temporal prior is constructed from a weighted, non-causal average of the values for each parameter. The relative influence of past and future timepoints is subject to weighting by a decaying exponential, with different integration times for the different ellipse parameters (e.g., we encourage the area of the pupil to change more gradually than its position).

There are many software parameters that control the behavior of the routines. While the default settings work well for some videos (including those that are part of the sandbox demo included in this repository), other parameter settings may be needed for videos with different qualities.

The transparentTrack toolbox has several dependencies, most notably the "quadfit" toolbox: https://www.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces

To properly install and configure these dependencies, install toolboxToolbox (tBtB), which provides for declarative dependency management for Matlab: https://github.com/ToolboxHub/ToolboxToolbox

Once tBtB is installed, transparentTrack can be readied for use with the command `tbUse('transparentTrack')`.
