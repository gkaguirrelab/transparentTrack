function eyeCenterOfRotation = calcEyeCenterOfRotation(transparentEllipse, eyeRadiusInPixels, varargin)
% candidatesEyeball = findCandidatesEyeball(transparentEllipse,eyeballRadius)
%
% this function returns the coordinates of the candidate eyeball centers
% for a given ellipse in transparent form. Note that if the ellipse has an
% eccentricity other than zero, 2 candidates Center of Projection will be
% returned. Further steps and/or consideration will help to determine which
% of the candidate is the "real" center of projection on the scene.
%
% Outputs:
%   eyeCenterOfRotation - matrix with the [X Y Z radius] coordinates of a
%       possible eyeball center on each row. Note that if only 1 couple of
%       coordinates is returned, the ellipse eccentricity is zero. In case
%       of orthogonal projection, the Z coordinate will be returned as
%       zero. The units are the same as the linear units used for the
%       transparent ellipse and the eyeball radius.
%
% Inputs:
%   transparentEllipse - ellipse in transparent form
%   eyeRadiusInPixels - radius of the eyeball, in the same linear coordinates
%       as in the transparent ellipse.
%
%% input parser
p = inputParser;

% required input
p.addRequired('transparentEllipse',@isnumeric);
p.addRequired('eyeRadiusInPixels',@isnumeric);

% Analysis parameters
p.addParameter('projectionModel','orthogonal',@ischar);


% parse
p.parse(transparentEllipse, eyeRadiusInPixels, varargin{:})


%% main


% find candidates center of projection for the ellipse
switch p.Results.projectionModel
    case 'orthogonal'
        if transparentEllipse(4) == 0 && (transparentEllipse(5) == 0 || transparentEllipse(5) == pi) % clear the easiest case
            eyeCenterOfRotation = [ transparentEllipse(1) transparentEllipse(2) 0 eyeRadiusInPixels];
            
        else
            % Find candidate azimuth and elevation values (unless the eccentricity
            % is zero, this will return 2 values for each angle).
            centerOfProjection = nan;
            [reconstructedPupilAzi, reconstructedPupilEle, ~] = pupilProjection_inv(transparentEllipse, centerOfProjection, p.Results.projectionModel);
            
            % if the pupilAzi or the pupilEle is zero
            if ~any(reconstructedPupilAzi) || ~any(reconstructedPupilEle)
                anglePairs = [reconstructedPupilAzi' reconstructedPupilEle'];
            else
                % look at the ellipse theta (assumed to vary from 0 to pi/2)  and
                % derive the plausible couples of angles for the center of projection
                if (transparentEllipse(5) <= pi/2 && transparentEllipse(5) >= 0)
                    for jj = 1: length(reconstructedPupilAzi)
                        anglePairs(jj,1) = reconstructedPupilAzi(jj);
                        eleIDX = find(sign(reconstructedPupilEle)~=sign(reconstructedPupilAzi(jj)));
                        anglePairs(jj,2) =  reconstructedPupilEle(eleIDX);
                    end
                else
                    for jj = 1: length(reconstructedPupilAzi)
                        anglePairs(jj,1) = reconstructedPupilAzi(jj);
                        eleIDX = find(sign(reconstructedPupilEle)==sign(reconstructedPupilAzi(jj)));
                        anglePairs(jj,2) =  reconstructedPupilEle(eleIDX);
                    end
                end
            end
            % find the eyeball center coordinates for this ellipse
            for ii = 1: length(anglePairs)
                delta(ii,:) = eyeRadiusInPixels * [ sind(anglePairs(ii,1)) sind(anglePairs(ii,2))];
                eyeCenterOfRotation (ii,:) = [transparentEllipse(1)-delta(ii,1) transparentEllipse(2)-delta(ii,2) 0];
            end
        end
    case 'perspective'
        error('not implemented yet');
end % switch

end % function


