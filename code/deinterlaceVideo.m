function deinterlaceVideo (params, dropboxDir, bobMode)

% This function allows to deinterlace NTSC DV 30Hz videos, saving out
% progressive 60 Hz videos, using a "bob deinterlacing" strategy.
%
% These deinterlace strategies are available (bobMode):
% 'Raw'    =  extract 2 fields for every frame. Save progressive video.
%             Final spatial resolution is half the original resolution.
% 'Zero'   =  extract 2 fields for every frame. Alternate every row with a
%             row of zeros to preserve aspect ratio.
% 'Double' =  extract 2 fields for every frame. Duplicate each raw to
%             preserve aspect ratio.
% 'Mean'   =  extract 2 fields for every frame. Add a row with the mean of
%             two consecutive rows to preserve aspect ratio.

% References on bob technique for deinterlacing and on deinterlacing in
% general:
% https://www.altera.com/content/dam/altera-www/global/en_US/pdfs/literature/wp/wp-01117-hd-video-deinterlacing.pdf
% http://www.100fps.com/
%
%   Usage:
%       deinterlaceVideo (params, dropboxDir, bobMode)
%
%   Required inputs:
%       dbDir
%       params.outputDir
%       params.projectFolder
%       params.projectSubfolder
%       params.eyeTrackingDir
%
%       params.subjectName
%       params.sessionDate
%       params.runName
%
% Note that the params field are the same as the metaData fields for a
% standard pupilResponse struct, so this function can also be used like
% this:
%       deinterlaceVideo (metadata, dropboxDir, bobMode)
%
%   Written by Giulia Frazzetta - Nov.2016


%% Set session and file names
if ~exist('bobMode', 'var')
    bobMode = 'Mean';
end

if isfield(params,'projectSubfolder')
    sessDir = fullfile(dropboxDir,params.projectFolder, params.projectSubfolder, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir);
    
    outDir = fullfile(dropboxDir,params.outputDir, params.projectSubfolder, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir);
else
    sessDir = fullfile(dropboxDir,params.projectFolder, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir);
    
    outDir = fullfile(dropboxDir,params.outputDir, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir);
end

if ~exist ('outDir', 'dir')
    mkdir (outDir)
end
%% Load video to deinterlace and set parameters for output video file.

inFile = fullfile(sessDir,[params.runName '_raw.mov']);
if ~exist(inFile,'file')
    inFile = fullfile(sessDir,[params.runName '.mov']);
end
inObj = VideoReader(inFile);
nFrames = floor(inObj.Duration*inObj.FrameRate);

Bob = VideoWriter(fullfile(outDir,[params.runName '_60hz.avi']));
Bob.FrameRate = inObj.FrameRate * 2;
Bob.Quality = 100;


%%
progBar = ProgressBar(nFrames,'Deinterlacing video...');
open(Bob)

switch bobMode
    case 'Raw'
        for i = 1:nFrames
            tmp = readFrame(inObj);
            thisFrame = rgb2gray(tmp);
            oddFields = thisFrame(1:2:end,:);
            evenFields = thisFrame(2:2:end,:);
            % shift the even lines to avoid "jumping" from frame to frame. (i.e
            % align the two fields)
            evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
            writeVideo(Bob,oddFields);
            writeVideo(Bob,evenFields);
            if ~mod(i,10);progBar(i);end
        end
        
    case 'Zero'
        for i = 1:nFrames
            tmp = readFrame(inObj);
            thisFrame = rgb2gray(tmp);
            oddFields = thisFrame(1:2:end,:);
            evenFields = thisFrame(2:2:end,:);
            % put zero rows in
            m = 1;
            k = 1;
            n = size(oddFields);
            oddFields = reshape([reshape(oddFields,m,[]);zeros(k,n(1)/m*n(2))],[],n(2));
            evenFields = reshape([reshape(evenFields,m,[]);zeros(k,n(1)/m*n(2))],[],n(2));
            evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
            writeVideo(Bob,oddFields)
            writeVideo(Bob,evenFields)
            if ~mod(i,10);progBar(i);end
        end
        
    case 'Double'
        for i = 1:nFrames
            tmp = readFrame(inObj);
            thisFrame = rgb2gray(tmp);
            oddFields = thisFrame(1:2:end,:);
            evenFields = thisFrame(2:2:end,:);
            % duplicate each row
            oddFields = repelem(oddFields, 2, 1);
            evenFields = repelem(evenFields, 2, 1);
            evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
            writeVideo(Bob,oddFields)
            writeVideo(Bob,evenFields)
            if ~mod(i,10);progBar(i);end
        end
        
    case 'Mean'
        for i = 1:nFrames
            tmp             = readFrame(inObj);
            thisFrame       = rgb2gray(tmp);
            oddFields = thisFrame(1:2:end,:);
            evenFields = thisFrame(2:2:end,:);
            % put means in between rows (odd fields)
            tmp = [oddFields(1,:); ((oddFields(1,:)+oddFields(2,:))/2);oddFields(2,:)];
            for jj = 2 : size(oddFields,1)-1
                newLines = [mean([oddFields(jj,:);oddFields(jj+1,:)],1);oddFields(jj+1,:)];
                tmp = cat(1,tmp,newLines);
            end
            oddFields = cat(1,tmp,oddFields(end,:));
            clear tmp
            clear newLines
            % put means in between rows (even fields)
            tmp = [evenFields(1,:); ((evenFields(1,:)+evenFields(2,:))./2);evenFields(2,:)];
            for jj = 2 : size(evenFields,1)-1
                newLines = [mean([evenFields(jj,:);evenFields(jj+1,:)],1);evenFields(jj+1,:)];
                tmp = cat(1,tmp,newLines);
            end
            evenFields = cat(1,evenFields(1,:),tmp);
            clear tmp
            clear newLines
            writeVideo(Bob,oddFields)
            writeVideo(Bob,evenFields)
            if ~mod(i,10);progBar(i);end
        end
end
close (Bob)