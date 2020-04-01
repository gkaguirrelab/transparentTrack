function updateHead(obj)

% Get the x and model
x = obj.x;
model = obj.model;

% Obtain the head position parameters
xHead = x(model.func.fieldSetIdx('head','all'));

% Create a matrix of x'y'z' locations of the camera 
A = obj.origRelCamPos;
A(isnan(A))=0;

% Shift and rotate the vector
B = (censorShift(A,xHead(1))'*rotMat(xHead(2),xHead(3),xHead(4))')';

obj.relCamPos = B;

end


function R = rotMat(alpha,beta,gamma)

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