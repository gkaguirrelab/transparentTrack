function saveRelCameraPos(obj,fileNameSuffix)
% Saves the current state of the relative camera position in the object
%
% Syntax:
%  obj.saveRelCameraPos()
%
% Description:
%   Save the updated relativeCameraPosition with the current state of the
%   parameters to disk, using the filepath provided by the videoStemName
%   input to the sceneObj class.
%
% Inputs:
%   none
%
% Outputs:
%   none
%


% Load or create a new structure
if exist([obj.videoStemName '_relativeCameraPosition.mat'], 'file') == 2
    load([obj.videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');
else
    return
end

% Update the values
relativeCameraPosition.values = obj.relCamPos;

% Update the meta data
relativeCameraPosition.meta.estimateSceneParams.p = obj.meta;
relativeCameraPosition.meta.estimateSceneParams.xHead = obj.x(obj.model.func.fieldSetIdx('head','all'));
relativeCameraPosition.meta.estimateSceneParams.model = obj.model;

% Save the  file
fileName = [obj.videoStemName '_relativeCameraPosition' fileNameSuffix '.mat'];
save(fileName,'relativeCameraPosition');


end




