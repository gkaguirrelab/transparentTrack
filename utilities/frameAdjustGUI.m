startPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session1_restAndStructure/TOME_3002/082616/EyeTracking/dMRI_dir98_AP_sceneGeometry.mat';

if exist(startPath)
    [path,file,suffix]=fileparts(startPath);
    file=[file suffix];
else
    [file,path] = uigetfile(fullfile(startPath,'*_sceneGeometry.mat'),'Choose a sceneGeometry file');
end

% Load the selected sceneGeometry file
sceneGeometryIn = fullfile(path,file);
dataLoad=load(sceneGeometryIn);
sceneGeometrySource=dataLoad.sceneGeometry;
clear dataLoad


fileStem = strsplit(file,'_sceneGeometry.mat');
fileStem = fileStem{1};
videoInFileName = fullfile(path,[fileStem '_gray.avi']);
fixedFrame = makeMedianVideoImage(videoInFileName);
blankFrame = ones(size(fixedFrame))*128;

fileList = dir(fullfile(path,'*_gray.avi'));
keep=cellfun(@(x) ~strcmp(x,[fileStem '_gray.avi']),extractfield(fileList,'name'));
fileList = fileList(keep);

fprintf('\n\nSelect the acquisition to adjust:\n')
for pp=1:length(fileList)
    optionName=['\t' num2str(pp) '. ' fileList(pp).name '\n'];
    fprintf(optionName);
end
fprintf('\nYou can enter a single acquisition number (e.g. 4),\n  a range defined with a colon (e.g. 4:7),\n  or a list within square brackets (e.g., [4 5 7]):\n')
choice = input('\nYour choice: ','s');
fileList = fileList(eval(choice));

figHandle = figure();
imshow(fixedFrame,[]);
hold on
title('\color{green}\fontsize{16}FIXED -- define canthus');
fprintf('Define the medial canthus triangle for the fixed image (lower, nasal, upper)\n');
[xF,yF] = ginput(3);

% Provide some instructions for the operator
fprintf('Adjust horizontal /vertical camera translation with the arrow keys.\n');
fprintf('Switch between moving and fixed image by pressing a.\n');
fprintf('Press esc to exit.\n\n');
fprintf([path '\n']);

% Loop over the selected acquisitions
for ff=1:length(fileList)
    
    videoInFileName = fullfile(path,fileList(ff).name);
    movingFrame = makeMedianVideoImage(videoInFileName);
    
    fprintf(fileList(ff).name);
    
    % Define the medial canthus for the moving image
    hold off
    imshow(movingFrame,[]);
    hold on
    title('\color{red}\fontsize{16}MOVING -- define canthus');
    [xM,yM] = ginput(3);
    
    % Enter a while loop
    showMoving = true;
    x = [0 0];
    notDoneFlag = true;
    while notDoneFlag
        hold off
        if showMoving
            movingImHandle = imshow(imtranslate(movingFrame,x,'method','cubic'),[]);
            hold on
            title('\color{red}\fontsize{16}MOVING');
        else
            fixedImHandle = imshow(fixedFrame,[]);
            hold on
            title('\color{green}\fontsize{16}FIXED');
        end
        
        % Plot the canthi
        triFixedHandle = plot(xF,yF,'-g');
        triMoveHandle = plot(xM+x(1),yM+x(2),'-r');
        
        keyAction = waitforbuttonpress;
        if keyAction
            keyChoiceValue = double(get(gcf,'CurrentCharacter'));
            switch keyChoiceValue
                case 28
                    text_str = 'translate left';
                    x(1)=x(1)-1;
                case 29
                    text_str = 'translate right';
                    x(1)=x(1)+1;
                case 30
                    text_str = 'translate up';
                    x(2)=x(2)-1;
                case 31
                    text_str = 'translate down';
                    x(2)=x(2)+1;
                case 97
                    text_str = 'swap image';
                    showMoving = ~showMoving;
                case 27
                    notDoneFlag = false;
                otherwise
                    text_str = 'unrecognized command';
            end
        end
    end
    
    % We are done. Update the figure window
    imshow(blankFrame)
    title('\color{black}\fontsize{16}Calculating camera translation');
    drawnow
    
    % Find the pupil center for the eye model in the fixed image
    eyePose = [0 0 0 3];
    pupilEllipse = pupilProjection_fwd(eyePose,sceneGeometrySource);
    targetPupilCenter = pupilEllipse(1:2)-x;
    
    % Now find the change in the extrinsic camera translation needed to
    % shift the eye model the observed number of pixels
    p0 = sceneGeometrySource.cameraPosition.translation;
    ub = sceneGeometrySource.cameraPosition.translation + [10; 10; 0];
    lb = sceneGeometrySource.cameraPosition.translation - [10; 10; 0];
    place = {'cameraPosition' 'translation'};
    mySG = @(p) setfield(sceneGeometrySource,place{:},p);
    pupilCenter = @(k) k(1:2);
    myError = @(p) norm(targetPupilCenter-pupilCenter(pupilProjection_fwd(eyePose,mySG(p))));
    options = optimoptions(@fmincon,'Diagnostics','off','Display','off');
    p = fmincon(myError,p0,[],[],[],[],lb,ub,[],options);
    
    % Report the value    
    fprintf(': camera translation [x,y,z] = [%2.2f; %2.2f; %2.2f] \n',p(1),p(2),p(3));
    
    % Save a sceneGeometry file with the adjusted extrinsic camera
    % translation
    fileStem = strsplit(fileList(ff).name,'_gray.avi');
    fileStem = fileStem{1};
    outFileName =  fullfile(path,[fileStem '_sceneGeometry.mat']);
    sceneGeometry = mySG(p);
    sceneGeometry.meta.frameAdjust.sourceSceneGeometry = sceneGeometryIn;
    sceneGeometry.meta.frameAdjust.xAdjust = x;
    save(outFileName,'sceneGeometry');
    
end
close(figHandle);
