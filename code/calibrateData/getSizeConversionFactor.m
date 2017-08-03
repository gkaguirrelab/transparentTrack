function sizeConversionFactor = getSizeConversionFactor(dotDataFilesNames,sizeGroundTruths, sizeConversionFileName, varargin)
% sizeConversionFactor = getSizeConversionFactor(sizeVideosNames,sizeGroundTruths)
%
% This routine will compute the size conversion factor from px to mm.
% If more than a single size video is used, the size conversion factor will
% be the median of the factors obtained by each of the size videos.
%
% OUTPUTS:
%   sizeConversionFactor: conversion factor for pupil size expressed in
%   pixels per mm [px/mm]
% INPUTS:
%   dotDataFilesNames: cell array containing the names of the dot data
%       files to be used.
%   sizeGroundTruths: array containing the ground truth for the dot size in
%       mm.
%   sizeConversionFileName: name of the mat file to save the size
%       conversion factor.
%
% Optional key/value pairs (analysis)
% 'calibrationType' - how to calibrate data (default: `useMajorAxis`)
% 
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%
%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('dotDataFilesNames',@(x) (iscell(x) | ischar(x)));
p.addRequired('sizeGroundTruths', @isnumeric);
p.addRequired('sizeConversionFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none', @ischar);


% parse
p.parse(dotDataFilesNames, sizeGroundTruths, varargin{:})

%% Check input

if iscell(dotDataFilesNames) && length(dotDataFilesNames)~=length(sizeGroundTruths)
    error('The number of data files and ground truth values does not match')
elseif ~iscell(dotDataFilesNames) && length(sizeGroundTruths) ~=1
    error('The number of data files and ground truth values does not match')
end

%% load data from each file



%% Derive conversion factor
for rr = 1: length(sizeGroundTruths)
    switch p.Results.calibrationType
        case 'horizontalOnly'
%             PXperMM.horizontal(rr) = nanmedian(dotsPX(rr).size) / sizeGroundTruths(rr);
            case 'verticalOnly'
%             PXperMM.vertical(rr) = nanmedian(dotsPX(rr).size) / sizeGroundTruths(rr);
            case 'horizontalAndVertical'
%             PXperMM.horizontal(rr) = nanmedian(dotsPX(rr).size) / sizeGroundTruths(rr);
%             PXperMM.vertical(rr) = nanmedian(dotsPX(rr).size) / sizeGroundTruths(rr);
    end
end
    sizeConversionFactor = [nanmedian(PXperMM.horizontal) nanmedian(PXperMM.vertical)];