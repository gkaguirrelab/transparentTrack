function [isGood] = checkGazeCal(LTdatFile, LTcalFile)
% checks that the calibration matrix is the result of a good Fixation test.
% if there are more than 7 measured points and the average error is less than
% 10 mm, we consider this a good calibration.

load (LTdatFile)
load (LTcalFile)

data = crsLiveTrackCalibrateRawData(CalMat, Rpc, pupil, glint);

errors = sqrt((targets(:,1)-data(:,1)).^2+(targets(:,2)-data(:,2)).^2); 
meanError = mean(errors(~isnan(errors))); %[mm]

if length(errors(~isnan(errors))) > 6 && meanError < 11
    isGood = 1;
else
    isGood = 0;
end