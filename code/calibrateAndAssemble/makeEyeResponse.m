function makeEyeResponse(timebase,pupil,gaze,eyeResponseFileName,varargin)
% Assembles the eye response file in its standard format
% 
% Description:
%  This function assebles the calibrated data in a standardized format
%  for further analysis.
%  The eyeResponse mat file contains a structure with the following fields:
%         eyeResponse.timebase
%         eyeResponse.pupil.width
%         eyeResponse.pupil.height
%         eyeResponse.pupil.area
%         eyeResponse.gaze.X
%         eyeResponse.gaze.Y
%         eyeResponse.gaze.ecc
%         eyeResponse.gaze.pol
%         eyeResponse.meta
%
% Input (required)
%   timebase 			- structure with fields that contain the timebase information
%       				  in milliseconds, and a meta field.
%   pupil 				- structure with fields that contain the pupil size information.
%       				  As default, this will be calibrated values expressed in mm, but the
%       				  routine will also accept uncalibrated pupil structs.
%   gaze 				- structure with fields that contain the pupil size information.
%       				  As default, this will be calibrated values, but the routine will
%       				  also accept uncalibrated gaze vectors.
%   eyeResponseFileName - name of the file in which to save the
%       				  eyeResponse struct.
% 
% Optional key/value pairs (analysis)
%   'uncalibratedData'  - toggle to true to add a flag for uncalibrated data
%
% Optional key/value pairs (environment)
%  'tbSnapshot'         - This should contain the output of the
%                         tbDeploymentSnapshot performed upon the result
%                         of the tbUse command. This documents the state
%                         of the system at the time of analysis.
%  'timestamp'          - AUTOMATIC; The current time and date
%  'username'           - AUTOMATIC; The user
%  'hostname'           - AUTOMATIC; The host
%
% Output (saved to file)
%	eyeResponse 		- structure with fields that contain the eyeResponse parameters.
% 

%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('timebase',@isstruct);
p.addRequired('pupil',@(x)(isempty(x) | isstruct(x)));
p.addRequired('gaze',@(x)(isempty(x) | isstruct(x)));
p.addRequired('eyeResponseFileName', @ischar);

% Optional analysis parameters
p.addParameter('uncalibratedData',false, @islogical);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(timebase,pupil,gaze,eyeResponseFileName,varargin{:})


%% sanity check
% verify that all the structs are columns of the same length

if ~isempty(pupil) && ~isempty(gaze)
    try
        tmp = cat(2,timebase.timebase,pupil.width,gaze.X);
        clear tmp
    catch
        error('eye Response parameters are not of the same length')
    end
end
%% asseble eyeResponse struct

eyeResponse.timebase = timebase.timebase;

if ~isempty(pupil)
    eyeResponse.pupil.width = pupil.width;
    eyeResponse.pupil.height = pupil.height;
    eyeResponse.pupil.area = pupil.area;
end

if ~isempty(gaze)
    eyeResponse.gaze.X = gaze.X;
    eyeResponse.gaze.Y = gaze.Y;
    eyeResponse.gaze.ecc = gaze.ecc;
    eyeResponse.gaze.pol = gaze.pol;
    eyeResponse.gaze.viewingDist = gaze.viewingDist;
end


%% add metadata
eyeResponse.meta = p.Results;
eyeResponse.meta.timebase = timebase.meta;
if ~isempty(pupil)
    eyeResponse.meta.pupil = pupil.meta;
end

if ~isempty(gaze)
    eyeResponse.meta.gaze = gaze.meta;
end

%% save out eyeResponse file

 save(eyeResponseFileName, 'eyeResponse');