function [trackedSamples,trackedPct] = countTrackedSamples (Report)
% [trackedSamples,trackedPct] = countTrackedSamples (Report)
% 
% Calculates how many of the acquired samples have been tracked by the
% LiveTrack algorithm. A sample is considered tracked if both the glint and
% the pupil position are acquired. 
% This small tool is useful to verify the quality of the LiveTrack tracking
% and if it can be used to align the transparentTrack tracked data if the
% eye video was acquired in the LiveTrack+VTOP setup.
% 
% Input
% Report - struct with data as it was tracked by the liveTrack, in the
% standard format returned by the liveTrack itself.
% 
% Outputs
% trackedSamples - number of the tracked samples
% trackedPct - percentage of the tracked samples in the video

Ch01Samples = [Report.PupilTracked_Ch01];
Ch02Samples = [Report.PupilTracked_Ch02];
allSamples = Ch01Samples + Ch02Samples;
trackedSamples = find (allSamples);
trackedPct = (length(trackedSamples)/length(allSamples))*100;
