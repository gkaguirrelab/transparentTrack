function candidatesEyeball = findCandidatesEyeball(transparentEllipse,eyeballRadius, varargin)
% candidatesEyeball = findCandidatesEyeball(transparentEllipse,eyeballRadius)
%
% this function returns the coordinates of the candidate eyeball centers
% for a given ellipse in transparent form. Note that if the ellipse has an
% eccentricity other than zero, 2 candidates Center of Projection will be
% returned. Further steps and/or consideration will help to determine which
% of the candidate is the "real" center of projection on the scene.
%
% Outputs:
%   candidatesEyeball - matrix with the [X Y Z radius] coordinates of a
%       possible eyeball center on each row. Note that if only 1 couple of
%       coordinates is returned, the ellipse eccentricity is zero. In case
%       of orthogonal projection, the Z coordinate will be returned as
%       zero. The units are the same as the linear units used for the
%       transparent ellipse and the eyeball radius.
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
p.addParameter('orthogonalProjection',true,@islogical)
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(transparentEllipse, eyeballRadius, varargin{:})

%% find candidates center of projection for the ellipse
if p.Results.orthogonalProjection
    if transparentEllipse(4) == 0 && (transparentEllipse(5) == 0 || transparentEllipse(5) == pi) % clear the easiest case
        candidatesEyeball = [ transparentEllipse(1) transparentEllipse(2) 0 eyeballRadius];
        
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
        % find the eyeball center coordinates for this ellipse
        for ii = 1: length(anglePairs)
            delta(ii,:) = eyeballRadius * [ sind(anglePairs(ii,1)) sind(anglePairs(ii,2))];
            candidatesEyeball (ii,:) = [transparentEllipse(1)-delta(ii,1) transparentEllipse(2)-delta(ii,2) 0 eyeballRadius];
        end
    end  
else
    % DEV placeholder - implement perspective correction case
end

end % function


