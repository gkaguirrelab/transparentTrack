classdef sceneObj < handle
    
    properties (Constant)
                
        % None

    end
    
    % Private properties
    properties (GetAccess=private)        

        % None 

    end
    
    % Fixed after object creation
    properties (SetAccess=private)

        % Stored inputs to the object
        videoStemName
        frameSet
        meta
        
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

        % The current model parameters, and the best parameters seen
        x
        xBest
        
        % The fVal for the current model, and the model performance
        fVal
        fValBest
        
        % The set of model components returned by calcGlintGazeError
        modelEyePose
        modelPupilEllipse
        modelGlintCoord
        modelPoseGaze
        modelVecGaze
        poseRegParams
        vectorRegParams
        rawErrors
        
        % These three are derived from the poseRegParams
        fixationEyePose
        screenTorsion
        screenRotMat
        
        % Verbosity
        verbose

        % The multi-scene objective can stash values here related to the
        % search across all scene objects
        multiSceneMeta
        multiSceneIdx
        
    end
    
    methods

        % Constructor
        function obj = sceneObj(videoStemName, frameSet, gazeTargets, setupArgs, keyVals, meta, varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('videoStemName',@ischar);
            p.addRequired('frameSet',@isnumeric);
            p.addRequired('gazeTargets',@isnumeric);
            p.addRequired('setupArgs',@iscell);
            p.addRequired('keyVals',@iscell);
            p.addRequired('meta',@isstruct);
            
            p.addParameter('verbose',false,@islogical);
        
            % parse
            p.parse(videoStemName, frameSet, gazeTargets, setupArgs, keyVals, meta, varargin{:})
                        
            
            %% Store inputs in the object
            obj.videoStemName = videoStemName;
            obj.frameSet = frameSet;
            obj.meta = meta;
            obj.verbose = p.Results.verbose;
            
            
            %% Initialize some properties
            obj.fValBest = Inf;
            obj.multiSceneMeta = [];
            
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
        updateScene(obj, x)
        fVal = calcError(obj, x)
        saveEyeModelMontage(obj,fileNameSuffix)
        saveModelFitPlot(obj,fileNameSuffix)
        saveSceneGeometry(obj,fileNameSuffix)
    
    end
end