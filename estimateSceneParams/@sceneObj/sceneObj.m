classdef sceneObj < handle
    
    properties (Constant)
                
        % foo
        foo = 3;

    end
    
    % Private properties
    properties (GetAccess=private)        

        % bar 
        bar = 4;

    end
    
    % Fixed after object creation
    properties (SetAccess=private)

        % Stored inputs to the object
        videoStemName
        frameSet
        
        % The arguments for the objective function, consisting of
        % 	[perimeter, glintData, ellipseRMSE, gazeTargets]
        args
        
        % These key values are passed to calcGlintGazeError
        keyVals
                        
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
                
        % The sceneGeometry that is being modeled
        sceneGeometry

        % The model parameters
        x
        
        % The fVal for the model
        fVal
        
        % The set of model components returned by calcGlintGazeError
        modelEyePose
        modelPupilEllipse
        modelGlintCoord
        modelPoseGaze
        modelVecGaze
        poseRegParams
        vectorRegParams
        rawErrors
        
        fixationEyePose
        screenTorsion
        screenRotMat
        
        % Verbosity
        verbose

    end
    
    methods

        % Constructor
        function obj = sceneObj(videoStemName, frameSet, gazeTargets, setupArgs, keyVals, varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('videoStemName',@ischar);
            p.addRequired('frameSet',@isnumeric);
            p.addRequired('gazeTargets',@isnumeric);
            p.addRequired('setupArgs',@iscell);
            p.addRequired('keyVals',@iscell);
            
            p.addParameter('verbose',true,@islogical);
        
            % parse
            p.parse(videoStemName, frameSet, gazeTargets, setupArgs, keyVals, varargin{:})
                        
            
            %% Store inputs in the object
            obj.videoStemName = videoStemName;
            obj.frameSet = frameSet;
            
            
            %% Create initial sceneGeometry structure
            obj.sceneGeometry = createSceneGeometry(setupArgs{:});
            
            
            %% Load the materials
            load([videoStemName '_correctedPerimeter.mat'],'perimeter');
            load([videoStemName '_glint.mat'],'glintData');
            load([videoStemName '_pupil.mat'],'pupilData');
                        
            % Extract the frames we want
            perimeter.data = perimeter.data(frameSet);
            glintData.X = glintData.X(frameSet); glintData.Y = glintData.Y(frameSet);
            ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSet);
            
            % Assemble these components into the args variable
            obj.args = {perimeter, glintData, ellipseRMSE, gazeTargets};
            
            % Store the keyVals
            obj.keyVals = keyVals;            
            
            % Done with these big variables
            clear perimeter glintData pupilData
                                    
        end
        
        % Required methds
        fVal = calcError(obj, x)
        saveEyeModelMontage(obj)
        saveModelFitPlot(obj)
        saveSceneGeometry(obj)

    end
end