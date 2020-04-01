function updateHead(obj)
% Update the description of the effect of head motion on camera position
%
% Syntax:
%  obj.updateHead()
%
% Description:
%   The sceneObj will attempt to load a representation of the position of
%   the head over the time of the acquisition of the frames of the scene.
%   The vectors are in the X, Y, Z world coordinate frame, and express head
%   motion as the position of the camera relative to the fixed head, and
%   relative to the start of the acquisition. These vectors of head motion
%   (and thus relative camera position) may be measured, or synthesized
%   (e.g., a set of linear drift vectors). The estimateSceneParams function
%   supports searching over parameters that transform the head motion
%   vectors to better account for the eye appearance in the frames.
%   Specifically, the coordinate frame in which head motion was measured
%   is not necessarily the same as the coordinate frame in which camera
%   position is specified relative to the head. This is certainly true for
%   the case in which head motion is derived from fMRI measures of brain
%   position over time--in which the coordinate frame is relative to the
%   axis of the scanner bore.
%
%   The parameters for this transformation are three Euler angles of
%   rotation (termed here alpha, beta, gamma) in units of degrees, and a
%   temporal phase shift in units of frames. The phase shift is required as
%   the measurement of headmotion may not be temporally locked with the
%   video of the eye.
%
% Inputs:
%   none
%
% Outputs:
%   none
%

% Get the x and model
x = obj.x;
model = obj.model;

% Obtain the head position parameters by reference to the model structure.
xHead = x(model.func.fieldSetIdx('head','all'));

% The matrix A contains the X, Y, Z vectors of the original relative camera
% position over time that was loaded at the time of instantiation of the
% object. We set any nan values to zero.
A = obj.origRelCamPos;
A(isnan(A))=0;

% Apply the temporal shift and rotation
B = (censorShift(A,xHead(1))'*rotMat(xHead(2),xHead(3),xHead(4))')';

% Store the updated version of the relative camera position.
obj.relCamPos = B;

end


%% LOCAL FUNCTIONS

function R = rotMat(alpha,beta,gamma)
% Return a rotation matrix for the passed Euler angles

Rz = [cosd(alpha) -sind(alpha) 0; sind(alpha) cosd(alpha) 0; 0 0 1];
Ry = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
Rx = [1 0 0; 0 cosd(gamma) -sind(gamma); 0 sind(gamma) cosd(gamma)];

R = Rz*Ry*Rx;

end

function vecOut = censorShift(vecIn,frameShift)
% For those points where we have shifted the data out of frame, set the
% head position correction to zero.
for jj = 1:size(vecIn,1)
    vecOut(jj,:) = fshift(vecIn(jj,:),frameShift);
    if frameShift<0
        vecOut(jj,end+round(frameShift):end)=0;
    end
    if frameShift>0
        vecOut(jj,1:round(frameShift))=0;
    end
end
end