function [distances, ellipsesCOPs] = distanceToCandidateEyeBallCenter(ellipses,candidateEyeBallCenterInImagePlanePixels,varargin)
%  [distances,ellipsesCOPs] = findCOPDistance(ellipses,sceneCOP)
%
% this function finds the distance of candidate center of projection for
% each ellipses and the center of projection on the scene.
%
%

%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipses',@isnumeric);
p.addRequired('candidateEyeBallRotationCenterInImagePlanePixels',@isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(ellipses,candidateEyeBallCenterInImagePlanePixels,varargin{:})

%% initialize variables
distances = nan(size(ellipses,1),2);
ellipsesCOPs = nan(size(ellipses,1),2);
%% loop through ellipses, find all distances and pick center of projection
for ii = 1:size(ellipses,1)
    if ~any(isnan(ellipses(ii,:)))
        % find candidates COP for this ellipse
        ellipseCOPCandidates = findCandidatesCenterOfProjection(ellipses(ii,:),candidateEyeBallCenterInImagePlanePixels(3));
        
        % find eucledian distance for each ellipse COP
        for jj = 1:length(ellipseCOPCandidates)
            distances(ii,jj) = sqrt((candidateEyeBallCenterInImagePlanePixels(1) - ellipseCOPCandidates(jj,1))^2 + (candidateEyeBallCenterInImagePlanePixels(2) - ellipseCOPCandidates(jj,2))^2);
        end
        % select ellipseCOP with the min distance from the scene COP and store
        % it
        [~,minDistIDX] =min(distances(ii,:));
        ellipsesCOPs(ii,:) = ellipseCOPCandidates(minDistIDX,:);
    else
        continue
    end
end

distances = min(distances');

end % function