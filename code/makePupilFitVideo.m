function [videoFileNameOut] = makePupilFitVideo(videoInFileName,videoOutFileName,ellipseFitDataFileName, varargin)
% [videoFileNameOut] = makePupilFitVideo(perimeterVideoFileName, varargin)
%
% This routine will superimpose ellipses on frames of a video file.
%
% INPUTS:
%   videoInFileName: full path to an .avi file. Points on the
%     boundary of the pupil should have a value of unity, and the frame
%     should be otherwise zero-filled. A frame that has no information
%     regarding the pupil (e.g., during a blink) should be zero-filled.
%   videoOutFileName: full path to an .avi file. Points on the
%     boundary of the pupil should have a value of unity, and the frame
%     should be otherwise zero-filled. A frame that has no information
%     regarding the pupil (e.g., during a blink) should be zero-filled.
%  ellipseFitDataFileName: full path to the .mat file that contains the
%     pupil tracking information.
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'display' - controls a display of the fitting outcome. [none, full]
%  'videoOutFrameRate' - frame rate (in Hz) of saved video. Default 30.
%  'whichFieldToPlot' - The name of the field within the ellipseFitData
%     struct that is to be plotted
%

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('perimeterVideoFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);
p.addParameter('display','none',@ischar);
p.addParameter('videoOutFrameRate',30,@isnumeric);
p.addParameter('whichFieldToPlot','whichFieldToPlot',@ischar);


%% Create a fit video if requested


    % Alert the user
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Creating and saving fit video. Started ' char(datetime('now')) '\n']);
        fprintf('| 0                      50                   100%% |\n');
        fprintf('.');
    end

    % Create a figure
    if strcmp(p.Results.display,'full')
        frameFig = figure( 'Visible', 'on');
    else
        frameFig = figure( 'Visible', 'off');
    end
    
    % Open the video out object
    outVideoObj = VideoWriter(p.Results.finalFitVideoOutFileName);
    outVideoObj.FrameRate = p.Results.videoOutFrameRate;
    open(outVideoObj);
    
    % Loop through the frames
    for ii=1:nFrames

        % Update the progress display
        if strcmp(p.Results.verbosity,'full')
            if mod(ii,round(nFrames/50))==0
                fprintf('.');
            end
        end
            
        % Plot the pupil boundary data points
        imshow(squeeze(pupilBoundaryData(:,:,ii)))
        
        if ~isnan(ellipseFitData.pPosteriorMeanTransparent(ii,1))
            % build ellipse impicit equation
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(ellipseFitData.pPosteriorMeanTransparent(ii,:)));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            
            % superimpose the ellipse using fimplicit
            hold on
            fimplicit(fh,[1, videoSizeY, 1, videoSizeX],'Color', 'green');
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis off;
        end
        
        % Write the frame to the file
        writeVideo(outVideoObj,getframe(frameFig));
        
    end
    
    % close the video object
    close(outVideoObj);
    
    % close the figure
    close(frameFig)
    
    % report completion of fit video generation
    if strcmp(p.Results.verbosity,'full')
        fprintf('\n');
        toc
    end


end

