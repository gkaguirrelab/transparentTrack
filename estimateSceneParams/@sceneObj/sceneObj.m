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
        model
        videoStemName
        frameSet
        gazeTargets
        setupArgs
        meta
        verbose
        keyVals
            
        % Fixed data used to guide the search
        perimeter
        glintDataX
        glintDataY
        ellipseRMSE        

        % The origRelCamPos vector (derived from the relativeCameraPosition
        % file
        origRelCamPos
                                
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
                
        % The sceneGeometry that is being modeled
        sceneGeometry

        % The relCameraPos, which is updated based upon search params
        relCamPos        
        
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
        
        % The multi-scene objective can stash values here related to the
        % search across all scene objects
        multiSceneMeta
        multiSceneIdx
        
    end
    
    methods

        % Constructor
        function obj = sceneObj(model, videoStemName, frameSet, gazeTargets, setupArgs, keyVals, meta, varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('model',@isstruct);
            p.addRequired('videoStemName',@ischar);
            p.addRequired('frameSet',@isnumeric);
            p.addRequired('gazeTargets',@isnumeric);
            p.addRequired('setupArgs',@iscell);
            p.addRequired('keyVals',@iscell);
            p.addRequired('meta',@isstruct);
            
            p.addParameter('verbose',false,@islogical);
        
            % parse
            p.parse(model, videoStemName, frameSet, gazeTargets, setupArgs, keyVals, meta, varargin{:})
                        
            
            %% Store inputs in the object
            obj.model = model;
            obj.videoStemName = videoStemName;
            obj.frameSet = frameSet;
            obj.gazeTargets = gazeTargets;
            obj.setupArgs = setupArgs;
            obj.meta = meta;
            obj.verbose = p.Results.verbose;
            obj.keyVals = keyVals;            

            
            %% Initialize some properties
            obj.fValBest = Inf;
            obj.multiSceneMeta = [];
            

            %% Create initial sceneGeometry structure
            obj.sceneGeometry = createSceneGeometry(setupArgs{:});
            
            
            %% Load the materials
            load([videoStemName '_correctedPerimeter.mat'],'perimeter');
            load([videoStemName '_glint.mat'],'glintData');
            load([videoStemName '_pupil.mat'],'pupilData');
            if exist([videoStemName '_relativeCameraPosition.mat'], 'file') == 2
                load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');                
            else
                relativeCameraPosition.values = zeros(3,max(frameSet));
            end
            
            % Extract data for the frames we want and store in the object
            obj.perimeter = perimeter.data(frameSet);
            obj.glintDataX = glintData.X(frameSet);
            obj.glintDataY = glintData.Y(frameSet);
            obj.ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSet);
            
            % We store the entire relative camera position vector, as we
            % will be shifting and interpolating this.
            obj.origRelCamPos = relativeCameraPosition.values;
            obj.relCamPos = obj.origRelCamPos;
            
            % Done with these big variables
            clear perimeter glintData pupilData relativeCameraPosition
                                    
        end
        
        % Required methds
        updateScene(obj)
        updateHead(obj)
        updateError(obj, varargin)
        fVal = returnError(obj, x)
        saveEyeModelMontage(obj,fileNameSuffix)
        saveModelFitPlot(obj,fileNameSuffix)
        saveSceneGeometry(obj,fileNameSuffix)
        
        % Private methods
        
    
    end
end