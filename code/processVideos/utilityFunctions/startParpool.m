function [ nWorkers ] = startParpool( nWorkers, tbtbRepoName, verbosity )
% Open and configure the parpool
%
% Description
%   Several stages of transparentTrack make use of the parpool. This
%   routine opens the par pool if it does not exit, configures it with
%   tbBtb, and returns the number
%
% Inputs:
%   nWorkers              - Scalar. The number of workers requested.
%   tbtbRepoName          - String, The repository to call for tBtB.
%   verbosity             - String. How verbose to be.
%
% Outputs:
%   nWorkers              - Scalar. The number of workers available.
%

% Silence the timezone warning
warningState = warning;
warning('off','MATLAB:datetime:NonstandardSystemTimeZoneFixed');

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
        % Use TbTb to configure the workers.
        if ~isempty(tbtbRepoName)
            spmd
                tbUse(tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(verbosity,'full')
                fprintf('CAUTION: Any TbTb messages from the workers will not be shown.\n');
            end
        end
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

