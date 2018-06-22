function sceneGeometry = loadSceneGeometry(sceneGeometryFileName, verbose)
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
%   verbose               - Logical. Defaults false if not passed.
%
% Outputs:
%   sceneGeometry         - The sceneGeometry structure
%

% Set the verbose flag to false if not passed
if nargin==1
    verbose = false;
end


if isempty(sceneGeometryFileName)
    sceneGeometry=[];
else
    % load the sceneGeometry structure
    dataLoad=load(sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    % instantiate the ray-tracing function
    if ~isempty(sceneGeometry.refraction)
        % The field is not empty, so we should have a ray tracing function.
        try
            funResult = functions(sceneGeometry.refraction.handle);
        catch
            error('Invalid function definition in sceneGeometry.refraction');
        end
        switch funResult.type
            case 'anonymous'
                % No further action needed
                if verbose
                    fprintf('Anonymous ray trace function available\n');
                end
            case 'simple'
                % Determine if the function exists
                if exist(func2str(sceneGeometry.refraction.handle))==0
                    % This might be a sceneGeometry structure created on a
                    % different computer, with a different path to the
                    % compiled virtual image function.
                    functionName = func2str(sceneGeometry.refraction.handle);
                    sceneGeometry.refraction.handle = eval(['@' functionName]);
                    sceneGeometry.refraction.path = eval(['which ' functionName]);
                    % Check to make sure that it is now available
                    if exist(func2str(sceneGeometry.refraction.handle))==0
                        error('Unable to re-instantiate the ray tracing function')
                    end
                end
                if verbose
                    fprintf('Compiled ray trace function available\n');
                end
            otherwise
                error('Unrecognized function type');
        end
    end
end

end