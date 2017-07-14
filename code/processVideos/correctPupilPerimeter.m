function correctPupilPerimeter(perimeterFileName, controlFileName, correctedPerimeterFileName, varargin)

% correctPupilPerimeterVideo applies the instructions from the control file
% on the pupil perimeter video. A new corrected perimeter video will be
% saved out in the specified file.
% 
% Each frame of the original perimeter video is loaded and elaborated
% according to the control file instructions (if any) for the given frame.
% 
% Currently available instructions include:
% - blink >> save out a black frame
% - bad >> save out a black frame
% - ellipse >> draw an ellipse with the given params into the frame
% - cut >> cut the perimeter using radius and theta given
% - error >> if any error occurred while compiling the automatic
%       instructions the frame won't be corrected, but an error flag will be
%       displayed for later inspection.
% 
% Note that each line of the control file is set of instructions for one
% specifical video frame, identified by the FrameNumber. If there is no
% instruction for a given frame, the frame will be saved as is. The control
% file may contain multiple instruction lines referred to the same frame
% (e.g. first do a vertical cut, then do an horizontal cut); in this case,
% the routine will process the instruction on the frame one after the other
% in the order they are presented.
% 
% Once the instructions have been applied the video will be saved out.
% 
% 
% Input (required)
%	perimeterFileName - path to the .mat file containing perimeter data.
%   controlFileName - path to the control file (with/without extestion)
%   correctedPerimeterFileName - path to the corrected perimeter data file
%
% Optional key/value pairs (verbosity and I/O)
%  'verbosity' - level of verbosity. [none, full]
%
% Optional key/value pairs (Environment parameters)
%  'tbSnapshot' - This should contain the output of the tbDeploymentSnapshot
%    performed upon the result of the tbUse command. This documents the
%    state of the system at the time of analysis.
%  'timestamp' - AUTOMATIC - The current time and date
%  'username' - AUTOMATIC - The user
%  'hostname' - AUTOMATIC - The host

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('perimeterFileName',@isstr);
p.addRequired('controlFileName',@isstr);
p.addRequired('correctedPerimeterFileName',@isstr);

% Optional display params
p.addParameter('verbosity','none',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(perimeterFileName, controlFileName, correctedPerimeterFileName, varargin{:})


%% Sanity check the inputs and load the files

% check controlFileName format
[~,~,ext] = fileparts(controlFileName);
if ~strcmp(ext,'.csv')
    error (' Only csv estension can be used')
end

% load control file
instructions = importControlFile(controlFileName);

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
originalPerimeter=dataLoad.perimeter;
clear dataLoad

% Set up some variables to guide the analysis and hold the result
nFrames=size(originalPerimeter.data,3);
blankFrame=originalPerimeter.data(:,:,1)*0;
perimeter=originalPerimeter;

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Correcting the perimeter file. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end


% loop through video frames
for ii = 1:nFrames

    % Update progress
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end

    % Proceed if there areinstructions for this frame
    instructionIdx = find ([instructions.frame] == ii);    
    if ~isempty(instructionIdx)
        
        % Obtain this frame
        img = squeeze(originalPerimeter.data(:,:,ii));

        for dd=1:length(instructionIdx)
            switch instructions(instructionIdx(dd)).type
                case 'blink'
                    img=blankFrame;
                case 'bad'
                    img=blankFrame;
                case 'error'
                    % do not modify img at all 
                case 'ellipse'
                    % get the instruction params
                    [cx, cy, a, b, phi] = parseInstructionParams(instructions(instructionIdx(dd)));
                    % start from back frame
                    img = blankFrame;
                    % find ellipse points
                    [Xe,Ye] = ellipse(N, cx, cy, a, b, phi);
                    Xe = round(Xe);
                    Ye = round(Ye);
                    % draw ellipse in frame
                    img(sub2ind(size(img),Ye(:),Xe(:))) = 1;
                case 'cut'
                    % get cut params
                    [radiusThresh,theta] = parseInstructionParams(instructions(instructionIdx(dd)));
                    % cut
                    [img] = cutPupil(img,radiusThresh,theta);
                otherwise
                    warning(['Instruction ' instructions(instructionIdx(dd)).type ' for frame ' num2str(ii) ' is unrecognized.']);
            end % switch instruction types
        end % loop over instructions

        % save modified frame
        perimeter.data(:,:,ii)=img;
    end % we have instructions for this frame
end % loop through frames

% save mat file with the video and analysis details
perimeter.meta = p.Results;
save(correctedPerimeterFileName,'perimeter');

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end


end % function