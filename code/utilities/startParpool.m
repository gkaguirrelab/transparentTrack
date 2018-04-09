function [ nWorkers ] = startParpool( nWorkers, verbosity )
% Open and configure the parpool
%
% Syntax:
%  [ nWorkers ] = startParpool( nWorkers, verbosity )
%
% Description:
%   Several stages of transparentTrack make use of the parpool. This
%   routine opens the parpool (if it does not currently exist) and returns
%   the number of available workers.
%
% Inputs:
%   nWorkers              - Scalar. The number of workers requested.
%   verbosity             - String. How verbose to be.
%
% Outputs:
%   nWorkers              - Scalar. The number of workers available.
%

% Silence the timezone warning
warningState = warning;
warning('off','MATLAB:datetime:NonstandardSystemTimeZoneFixed');
warning('off','MATLAB:datetime:NonstandardSystemTimeZone');

% If a parallel pool does not exist, attempt to create one
poolObj = gcp('nocreate');
if isempty(poolObj)
    if strcmp(verbosity,'full')
        tic
        fprintf(['Opening parallel pool. Started ' char(datetime('now')) '\n']);
    end
    if isempty(nWorkers)
        parpool;
    else
        parpool(nWorkers);
    end
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
    end
    if strcmp(verbosity,'full')
        toc
        fprintf('\n');
    end
else
    nWorkers = poolObj.NumWorkers;
end

% Restore the warning state
warning(warningState);

end % function -- startParpool

