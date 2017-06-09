function [trackedSamples,trackedPct] = CountTrackedSamples (Report)
% Calculates how many of the acquired samples have been tracked by the
% LiveTrack algorithm. A sample is considered tracked if both the glint and
% the pupil position are acquired.

Ch01Samples = [Report.PupilTracked_Ch01];
Ch02Samples = [Report.PupilTracked_Ch02];
allSamples = Ch01Samples + Ch02Samples;
trackedSamples = find (allSamples);
trackedPct = (length(trackedSamples)/length(allSamples))*100;
