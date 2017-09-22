function makeEyeResponse(timebase,pupil,gaze,eyeResponseFileName,varargin)
% 
% makeEyeResponse(timebase,pupil,gaze,eyeResponseFileName)
% 
% header
% 
%  The eyeResponse struct includes the following fields:
%         eyeResponse.timebase
%         eyeResponse.pupil.width
%         eyeResponse.pupil.height
%         eyeResponse.pupil.area
%         eyeResponse.gaze.X
%         eyeResponse.gaze.Y
%         eyeResponse.gaze.ecc
%         eyeResponse.gaze.pol
%         eyeResponse.meta
% % Output (saved to file)
%	eyeResponse - structure with fields that contain the eyeResponse parameters.
%
% Input (required)
%   timebase - structure with fields that contain the timebase information
%       in milliseconds, and a meta field.
%   pupil - structure with fields that contain the pupil size information.
%       As default, this will be calibrated values expressed in mm, but the
%       routine will also accept uncalibrated pupil structs.
%   gaze - structure with fields that contain the pupil size information.
%       As default, this will be calibrated values, but the routine will
%       also accept uncalibrated gaze vectors.
%   eyeResponseFileName - name of the file in which to save the
%       eyeResponse struct.
% 
% Options (analysis)
%   uncalibratedData - 
%   
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
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