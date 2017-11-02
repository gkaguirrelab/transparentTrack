function [distances, ellipsesCoPs] = distanceToCandidateEyeBallCenter(ellipses,candidateEyeball,varargin)
%  [distances,ellipsesCoPs] = distanceToCandidateEyeBallCenter(ellipses,candidateEyeball,varargin)
%
% this function finds the distance of candidate center of projection for
% each ellipses and the center of projection on the scene. Ellipses and
% candidate Eyeball must be in the same units.
% (e.g. Pixels)
%
% input 
%  ellipses - set of ellipses in transparent form.
%  candidateEyeball -[X Y Z R] coordinates of the center and radius of the eyeball. If
%       orthogonal projection, the Z coordinate (distance from the scene plane,
%       will be set to zero).

%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipses',@isnumeric);
p.addRequired('candidateEyeball',@isnumeric);

% Analysis parameters
p.addParameter('orthogonalProjection',true,@islogical);
% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(ellipses,candidateEyeball,varargin{:})


if p.Results.orthogonalProjection
%% initialize variables
distances = nan(size(ellipses,1),2);
ellipsesCoPs = nan(size(ellipses,1),4);
%% loop through ellipses, find all distances and pick center of projection
for ii = 1:size(ellipses,1)
    if ~any(isnan(ellipses(ii,:)))
        % find candidates Center of rotation
        ellipseCoRCandidates = findCandidatesEyeball(ellipses(ii,:),candidateEyeball(4));
        
        % find eucledian distance for each ellipse COP
        for jj = 1:size(ellipseCoRCandidates,1)
            distances(ii,jj) = sqrt((candidateEyeball(1) - ellipseCoRCandidates(jj,1))^2 + (candidateEyeball(2) - ellipseCoRCandidates(jj,2))^2 + (candidateEyeball(3) - ellipseCoRCandidates(jj,3))^2);
        end
        % select ellipseCoR with the min distance from the scene CoR and store
        % it
        [~,minDistIDX] =min(distances(ii,:));
        ellipsesCoPs(ii,:) = ellipseCoRCandidates(minDistIDX,:);
    else
        continue
    end
end

distances = min(distances');

else
    % Dev PLACEHOLDER - develop perspective correction case
end

end % function