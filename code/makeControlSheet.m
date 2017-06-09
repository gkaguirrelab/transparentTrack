function makeControlSheet( params )
%makeControlSheet: produces .xls and csv files from the cutting params
%   params are: savepath, trackfile

params.trackfile = '~/Desktop/GoodTrack/tfMRI_RETINO_PA_run01_testHC.mat';
params.runName = 'tfMRI_RETINO_PA_run01';
params.savePath = '~/Desktop/GoodTrack/';
%% load trackfile
load(params.trackfile);

%% make excel file
xlsFile = fullfile (params.savePath, [params.runName '.xlsx']);

% SHEET ONE - DATA

% make first line as the header:
% FRAME | CUT | DISTANCE METRIC | CIRCULARITY | AREA |
pupilHeader = {'FRAME' 'CUT'  'DISTANCE METRIC' 'AXES RATIO' 'AREA'};

% initialize data cell
pupilCell = NaN(length(pupil.X),5);

% fill in data cell
pupilCell(:,1) = 0:1:length(pupil.X)-1; % note that first frame is counted as zero on quicktime 7

for ff =1:length(pupil.X)
    if ~isempty (pupil.bestCut{ff})
        % get the error parameters
        pupilCell(ff,3) = pupil.(pupil.bestCut{ff}).distanceError(ff);
        pupilCell(ff,4) = pupil.(pupil.bestCut{ff}).axRatio(ff);
        pupilCell(ff,5) = pupil.(pupil.bestCut{ff}).area(ff);
        % assign a code to the cut
        if strcmp (pupil.bestCut(ff), 'full')
            pupilCell(ff,2) = 1;
        elseif strcmp (pupil.bestCut(ff), 'cut75up')
            pupilCell(ff,2) = 2;
        elseif strcmp (pupil.bestCut(ff), 'cut75right')
            pupilCell(ff,2) = 3;
        elseif strcmp (pupil.bestCut(ff), 'cutDiagonal4075')
            pupilCell(ff,2) = 4;
        else
            pupilCell(ff,2) = nan;
        end
    else
        continue
    end
end




end