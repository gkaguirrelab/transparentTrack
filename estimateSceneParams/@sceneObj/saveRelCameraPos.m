function saveRelCameraPos(obj)
% Saves the current state of the relative camera position in the object
%
% Syntax:
%  obj.saveRelCameraPos()
%
% Description:
%   Save the relativeCameraPosition with the current state of the
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
if exist([videoStemName '_relativeCameraPosition.mat'], 'file') == 2
    load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');
else
    relativeCameraPosition.values = zeros(3,max(frameSet));
end

% Update the values
relativeCameraPosition.values = obj.relCamPos;

% Update the meta data
relativeCameraPosition.meta.estimateSceneParams.p = obj.meta;
relativeCameraPosition.meta.estimateSceneParams.xHead = obj.x(model.func.fieldSetIdx('head','all'));
relativeCameraPosition.meta.estimateSceneParams.model = obj.model;

% Save the  file
fileName = [obj.videoStemName '_relativeCameraPosition.mat'];
save(fileName,'relativeCameraPosition');


end




