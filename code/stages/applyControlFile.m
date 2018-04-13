function perimeter = applyControlFile(perimeterFileName, controlFileName, correctedPerimeterFileName, varargin)
% Modify a perimeter file by applying instructions from a control File
%
% Syntax:
%  perimeter = applyControlFile(perimeterFileName, controlFileName, correctedPerimeterFileName)
%
% Description:
%   This routine applies the instructions from the control file to the
%   pupil perimeter video. A new corrected perimeter video will be saved
%   out in the specified file.
% 
%   Each frame of the original perimeter video is loaded and elaborated
%   according to the control file instructions (if any) for the given
%   frame.
% 
%   Currently available instructions include:
%     - blink: save out a black frame
%     - bad: save out a black frame
%     - glintPatch: cut the portion of the perimeter that intersects the glint
%     - cut: cut the perimeter using radius and theta given
%     - ellipse: draw an ellipse with the given params into the frame
%     - error: if any error occurred while compiling the automatic
%       instructions the frame won't be corrected, but an error flag will be
%       displayed for later inspection.
% 
%   Note that each line of the control file is set of instructions for one
%   specifical video frame, identified by the FrameNumber. If there is no
%   instruction for a given frame, the frame will be saved as is. The
%   control file may contain multiple instruction lines referred to the
%   same frame (e.g. first do a vertical cut, then do an horizontal cut);
%   in this case, the routine will process the instruction on the frame one
%   after the other in the order they are presented.
% 
%   Once the instructions have been applied the video will be saved out.
% 
% Inputs:
%	controlFileName       - Full path (including csv extension) to the
%                           control file 
%   perimeterFileName     - Full path to the perimeter file
%   correctedPerimeterFileName - Full path to the corrected perimeter file
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Boolean. Default false.
%
% Optional key/value pairs (environment)
%  'tbSnapshot'           - This should contain the output of the
%                           tbDeploymentSnapshot performed upon the result
%                           of the tbUse command. This documents the state
%                           of the system at the time of analysis.
%  'timestamp'            - AUTOMATIC; The current time and date
%  'username'             - AUTOMATIC; The user
%  'hostname'             - AUTOMATIC; The host
%
% Outputs:
%   None. Saves the correctedPerimeterFile
%

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@isstr);
p.addRequired('controlFileName',@isstr);
p.addRequired('correctedPerimeterFileName',@isstr);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional environment params
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
instructions = loadControlFile(controlFileName);

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
originalPerimeter=dataLoad.perimeter;
clear dataLoad

% Set up some variables to guide the analysis and hold the result
nFrames=size(originalPerimeter.data,1);
perimeter = struct();
perimeter.size = originalPerimeter.size;
perimeter.data = cell(nFrames,1);
blankFrame = uint8(zeros(perimeter.size));

% alert the user
if p.Results.verbose
    tic
    fprintf(['Correcting the perimeter file. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end


% loop through video frames
for ii = 1:nFrames

    % Update progress
    if p.Results.verbose && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % get the data frame
    thisFrame = uint8(zeros(originalPerimeter.size));
    thisFrame(sub2ind(originalPerimeter.size,originalPerimeter.data{ii}.Yp,originalPerimeter.data{ii}.Xp))=255;

    % Proceed if there are instructions for this frame
    instructionIdx = find ([instructions.frame] == ii);    
    if ~isempty(instructionIdx)
        
        for dd=1:length(instructionIdx)
            switch instructions(instructionIdx(dd)).type
                case 'blink'
                    thisFrame=blankFrame;
                case 'bad'
                    thisFrame=blankFrame;
                case 'error'
                    thisFrame=blankFrame;
                case 'ellipse'
                    % get the instruction params
                    [p1, p2, p3, p4, p5] = parseControlInstructions(instructions(instructionIdx(dd)));
                    pE = ellipse_transparent2ex(p1, p2, p3, p4, p5);
                    % start from blank frame
                    thisFrame = blankFrame;
                    % find ellipse points
                    [Xe,Ye] = ellipse(pE(1),pE(2),pE(3),pE(4),pE(5));
                    Xe = round(Xe);
                    Ye = round(Ye);
                    % draw ellipse in frame
                    thisFrame(sub2ind(size(thisFrame),Ye(:),Xe(:))) = 1;
                case 'cut'
                    % get cut params
                    [radiusThresh,theta] = parseControlInstructions(instructions(instructionIdx(dd)));
                    [thisFrame] = applyPupilCut(thisFrame,radiusThresh,theta);
                case 'glintPatch'
                    % get cut params
                    [glintX,glintY,glintPatchRadius] = parseControlInstructions(instructions(instructionIdx(dd)));
                    % apply patch
                    glintPatch = ones(size(thisFrame));
                    glintPatch = insertShape(glintPatch,'FilledCircle',[glintX glintY glintPatchRadius],'Color','black');
                    glintPatch = im2bw(glintPatch);
                    thisFrame = immultiply(thisFrame,glintPatch);
                otherwise
                    warning(['Instruction ' instructions(instructionIdx(dd)).type ' for frame ' num2str(ii) ' is unrecognized.']);
            end % switch instruction types
        end % loop over instructions
    end % we have instructions for this frame
        
    % save the frame, which may include modifications
    [perimeter.data{ii}.Yp, perimeter.data{ii}.Xp] = ind2sub(size(thisFrame),find(thisFrame));

end % loop through frames

% save mat file with the video and analysis details
perimeter.meta = p.Results;
save(correctedPerimeterFileName,'perimeter');

% report completion of analysis
if p.Results.verbose
    fprintf('\n');
    toc
    fprintf('\n');
end


end % function