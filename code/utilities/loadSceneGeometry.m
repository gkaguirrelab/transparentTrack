function sceneGeometry = loadSceneGeometry(sceneGeometryFileName, verbosity)
% Load sceneGeometry and instantiate the ray tracing function
%
% Syntax:
%  sceneGeometry = loadSceneGeometry(sceneGeometryFileName, verbosity)
%
% Description:
%   Loads the sceneGeometry file from file and ensures that a valid
%   virtual image ray-tracing function is available.
%
% Inputs:
%   sceneGeometryFileName - Full path to the sceneGeometry file. If left
%                           empty, then an empty variable will be returned.
%   verbosity             - Character vector. A valid value is "full".
%
% Outputs:
%   sceneGeometry         - The sceneGeometry structure
%

if isempty(sceneGeometryFileName)
    sceneGeometry=[];
else
    % load the sceneGeometry structure
    dataLoad=load(sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    % instantiate the ray-tracing function
    if ~isempty(sceneGeometry.virtualImageFunc)
        % The field is not empty, so we should have a ray tracing function.
        try
            funResult = functions(sceneGeometry.virtualImageFunc.handle);
        catch
            error('Invalid function definition in sceneGeometry.virtualImageFunc');
        end
        switch funResult.type
            case 'anonymous'
                % No further action needed
                if strcmp(verbosity,'full')
                    fprintf('Anonymous ray trace function available\n');
                end
            case 'simple'
                % Determine if the function exists
                if exist(func2str(sceneGeometry.virtualImageFunc.handle))==0
                    % We need to add the function back to the path
                    addpath(fileparts(sceneGeometry.virtualImageFunc.path),'-BEGIN')
                    % Check to make sure that it is now available
                    if exist(func2str(sceneGeometry.virtualImageFunc.handle))==0
                        error('Unable to re-instantiate the ray tracing function')
                    end
                else
                    % If we have a compiled function, make sure that it is
                    % the right compiled function.
                    if exist(func2str(sceneGeometry.virtualImageFunc.handle))==3
                        if ~strcmp(which(func2str(sceneGeometry.virtualImageFunc.handle)), sceneGeometry.virtualImageFunc.path)
                                error('The available virtualImageFunc does not match that in sceneGeometry')
                        end
                    end
                end
                if strcmp(verbosity,'full')
                    fprintf('Compiled ray trace function available\n');
                end
            otherwise
                error('Unrecognized function type');
        end
    end
end

end