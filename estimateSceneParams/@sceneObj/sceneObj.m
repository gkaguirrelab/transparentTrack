classdef sceneObj < handle
% Class definition for sceneObj
%
% Syntax:
%  obj = sceneObj(model, videoStemName, frameSet, gazeTargets, setupArgs, meta)
%
% Description:
%   An object that holds information regarding a scene. A scene is a
%   collection of measurements made from the eye during one video
%   acquisition. The inputs specify a video, a set of frames from the video
%   during which measurements were made, and a set of positions on a screen
%   that the eye was fixated upon for each of the frames. A set of
%   measurements for each frame, derived from the video, are loaded from
%   disk (including the perimeter of the pupi, the ellipse fit to the pupil
%   boundary, and the position of one or more glints from an active light
%   source.
%
% Inputs:
%   model                 - Structure. The model search parameters and
%                           properties created by defineModelParams.m
%	videoStemName         - Char vector, or cell array of n char vectors.
%                           Full path to the n video files from which the
%                           scene observations have been derived. The stem
%                           name should omit the "_gray.avi" suffix that is
%                           usually present in the names of these video
%                           files. The routine assumes that files to be
%                           used in the analysis have the suffixes:
%                               {'_correctedPerimeter.mat','_glint.mat',...
%                                '_pupil.mat'}
%                           and optionally {'_headMotion.mat'}. The
%                           sceneGeometry output files will be saved at the
%                           same path location with the suffix:
%                               '_sceneGeometry.mat'
%   frameSet              - A 1xm vector that specifies the m frame indices
%                           (indexed from 1) which identify the set of
%                           frames from the scene to guide the search.
%   gazeTargets           - A 2xm matrix that provides the positions, in 
%                           degrees of visual angle, of fixation targets
%                           that correspond to each of the frames. The
%                           visual angle of the stimuli should be adjusted
%                           for min/magnification produced by spectacle
%                           lenses worn by the subject prior to being
%                           passed to this routine.
%   setupArgs             - Cell array. Key-value pairs used in the call
%                           to createSceneGeometry.
%   errorArgs             - Cell array. These are key-value pairs that are
%                           passed to method 'updateError'.
%   meta                  - Structure. Usually the p.Results output of the
%                           parameter parser from the calling function.
%   
% Outputs:
%   obj                   - Handle to the object.
%
    
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
            
        % Fixed data used to guide the search
        perimeter
        glintDataX
        glintDataY

        % The origRelCamPos vector (derived from the relativeCameraPosition
        % file
        origRelCamPos
                                
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
                
        % The sceneGeometry that is being modeled
        sceneGeometry
        
        % The aguments passed to updateError
        errorArgs

        % The relCameraPos, which is updated based upon search params
        relCamPos        
        
        % The current model parameters, the last, and best parameters seen
        x
        xLast
        xBest
        xHist
        
        % The fVal for the current model, and the model performance
        fVal
        fValBest
        fValHist
        
        % The set of model components returned by calcGlintGazeError
        modelEyePose
        modelCameraTrans
        modelPupilEllipse
        modelGlintCoord
        modelPoseGaze
        modelVecGaze
        poseRegParams
        vecRegParams
        rawErrors
                
        % The multi-scene objective can stash values here related to the
        % search across all scene objects
        multiSceneMeta
        multiSceneIdx
        
    end
    
    methods

        % Constructor
        function obj = sceneObj(model, videoStemName, frameSet, gazeTargets, setupArgs, meta, varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('model',@isstruct);
            p.addRequired('videoStemName',@ischar);
            p.addRequired('frameSet',@isnumeric);
            p.addRequired('gazeTargets',@isnumeric);
            p.addRequired('setupArgs',@iscell);
            p.addRequired('meta',@isstruct);
            
            % Optional
            p.addParameter('glintData',[],@(x)(isempty(x) || isstruct(x)));
            p.addParameter('perimeter',[],@(x)(isempty(x) || isstruct(x)));
            p.addParameter('relativeCameraPosition',[],@(x)(isempty(x) || isstruct(x)));
            p.addParameter('verbose',false,@islogical);
        
            % parse
            p.parse(model, videoStemName, frameSet, gazeTargets, setupArgs, meta, varargin{:})
                        
            
            %% Store inputs in the object
            obj.model = model;
            obj.videoStemName = videoStemName;
            obj.frameSet = frameSet;
            obj.gazeTargets = gazeTargets;
            obj.setupArgs = setupArgs;
            obj.meta = meta;
            obj.verbose = p.Results.verbose;

            
            %% Initialize some properties
            obj.fValBest = Inf;
            obj.multiSceneMeta = [];
            obj.errorArgs = {};
            

            %% Create initial sceneGeometry structure
            obj.sceneGeometry = createSceneGeometry(setupArgs{:});
            
            
            %% Load the materials
            if isempty(p.Results.perimeter)
                load([videoStemName '_correctedPerimeter.mat'],'perimeter');
            else
                perimeter = p.Results.perimeter;
            end
            if isempty(p.Results.glintData)
                load([videoStemName '_glint.mat'],'glintData');
            else
                glintData = p.Results.glintData;
            end

            % Sanity check the intputs, derive the total number of frames
            if length(perimeter.data) ~= length(glintData.X)
                error('The perimeter and gint input data are unequal in length');
            else
                nFrames = length(perimeter.data);
            end
            if max(frameSet)>nFrames
                error('The FrameSet contains one or more indices that exceed the number of available frames');
            end
            if isempty(p.Results.relativeCameraPosition)
                if exist([videoStemName '_relativeCameraPosition.mat'], 'file') == 2
                    load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');
                    if ~isfield(relativeCameraPosition,'initial')
                        relativeCameraPosition.initial.values = zeros(3,nFrames);
                        relativeCameraPosition.currentField = 'initial';
                    end
                else
                    relativeCameraPosition.initial.values = zeros(3,nFrames);
                    relativeCameraPosition.currentField = 'initial';
                end
            else
                relativeCameraPosition = p.Results.relativeCameraPosition;
            end
            
            % Extract data for the frames we want and store in the object
            obj.perimeter = perimeter.data(frameSet);
            obj.glintDataX = glintData.X(frameSet);
            obj.glintDataY = glintData.Y(frameSet);
            
            % We store the entire relative camera position vector, as we
            % will be shifting and interpolating this.
            obj.origRelCamPos = relativeCameraPosition.initial.values;
            obj.relCamPos = obj.origRelCamPos;
            
            % Done with these big variables
            clear perimeter glintData pupilData relativeCameraPosition
                                    
        end
        
        % Required methds
        updateScene(obj)
        updateHead(obj)
        updateError(obj, varargin)
        fVal = returnError(obj, x, errorArgs)
        saveEyeModelMontage(obj,fileNameSuffix,showImageFeatures,showModelEye,showAllElements)
        saveModelFitPlot(obj,fileNameSuffix)
        saveSceneGeometry(obj,fileNameSuffix)
        saveRelCameraPos(obj,fileNameSuffix)
        saveRelCameraPosPlot(obj,fileNameSuffix)
        
        % Private methods
        
    
    end
end