function cleanupMatlabPrefs    
% this small function closes any open parallel pool and removes the
% matlabprefs.mat file in case it is corrupt.
% This file can get corrupt while using parallel pools because all of the
% pool worker try to read/write it at the same time. MathWorks has been
% contacted about this issue.

% close parallel pool if any is open
poolobj = gcp('nocreate');
if ~isempty (poolobj)
    delete(poolobj);
end

% clean up the corrupt matlabprefs file
prefFile = fullfile(prefdir,'matlabprefs.mat');
system(['rm -rf "' prefFile '"'])
