function candidatesCOP = findCandidatesCenterOfProjection(transparentEllipse,eyeballRadius, varargin)
% candidatesCOP = findCandidatesCenterOfProjection(transparentEllipse,eyeballRadius)
%
% this function returns the candidate center of projections for a given
% ellipse in transparent form. Note that if the ellipse has an eccentricity
% other than zero, 2 candidates Center of Projection will be returned.
% Further steps and/or consideration will help to determine which of the
% candidate is the "real" center of projection on the scene.
%
% Outputs:
%   candidatesCOP - matrix with the [X Y] coordinates of a center of
%       projection candidates on each row. Note that if only 1 coupe of
%       coordinates is return, the ellipse eccentricity is zero. The units
%       are the same as the linear units used for the transparent ellipse
%       and the eyeball radius.
%
% Inputs:
%   transparentEllipse - ellipse in transparent form
%   eyeballRadius - radius of the eyeball, in the same linear coordinates
%       as in the transparent ellipse.
%
%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('transparentEllipse',@isnumeric);
p.addRequired('eyeballRadius',@isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(transparentEllipse, eyeballRadius, varargin{:})

%% find candidates center of projection for the ellipse

if transparentEllipse(4) == 0 && (transparentEllipse(5) == 0 || transparentEllipse(5) == pi) % clear the easiest case
    candidatesCOP = [ transparentEllipse(1) transparentEllipse(2)];
    
else
    % Find candidate azimut and elevation values (unless the eccentricity
    % is zero, this will return 2 values for each angle).
    centerOfProjection = nan;
    [reconstructedPupilAzi, reconstructedPupilEle, ~] = pupilProjection_inv(transparentEllipse,centerOfProjection);
    
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
    % find the candidates center of projection for this ellipse
    for ii = 1: length(anglePairs)
        delta(ii,:) = eyeballRadius * [ sind(anglePairs(ii,1)) sind(anglePairs(ii,2))];
        candidatesCOP (ii,:) = [transparentEllipse(1)-delta(ii,1) transparentEllipse(2)-delta(ii,2)];
    end
end


