function params = thresholdPupilVideo(grayI, params)

% this function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only

%% set default params

% params for circleFit (always needed)
if ~isfield(params,'rangeAdjust')
    params.rangeAdjust = 0.05;
end
if ~isfield(params,'circleThresh')
    params.circleThresh = [0.06 0.999];
end
if ~isfield(params,'pupilRange')
    params.pupilRange   = [20 90];
end
if ~isfield(params,'glintRange')
    params.glintRange   = [10 30];
end
if ~isfield(params,'glintOut')
    params.glintOut     = 0.1;
end
if ~isfield(params,'sensitivity')
    params.sensitivity  = 0.99;
end

% params for ellipse fit
if ~isfield(params,'ellipseThresh')
    params.ellipseThreshPupil = 0.97;
end
if ~isfield(params,'maskBox')
    params.maskBox   = [4 30];
end
if ~isfield(params,'gammaCorrection')
    params.gammaCorrection   = 1;
end



%% extract pupil perimeter
ih = figure;

for i = 1:numFrames
    
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
    % Track with circles
    % we assign a starting glint range for the circle fit routine. Glint
    % data won't be stored at this point.
    [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange);
    
    
    if isempty(pCenters) %no pupil circle
        % just save frame
        if isfield(params,'outVideo')
            frame   = getframe(ih);
            writeVideo(outObj,frame);
        end
        if ~mod(i,10);progBar(i);end;
        continue
    else
        % get pupil perimeter
        [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);
        
        imshow(binP)
        
        frame   = getframe(ih);
        writeVideo(outObj,frame);
    end
end