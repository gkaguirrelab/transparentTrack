function gkaModelEyeLocalHook
% 
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'transparentTrack'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
toolboxName = 'transparentTrack';
 
%% Clear out old preferences
if (ispref(toolboxName))
    rmpref(toolboxName);
end

%% Check for required Matlab toolboxes
% These are the toolboxes required for transparentTrack that are not
% otherwise required as part of gkaModelEye. Below is a list of toolbox
% license names. In the following assignment, the license name is given in
% the comment string after the matching version name for each toolbox.
requiredAddOns = {...
    'Computer Vision Toolbox',...                  % video_and_image_blockset
    'Image Processing Toolbox',...                 % image_toolbox
    'Parallel Computing Toolbox',...               % distrib_computing_toolbox
    };
% Given this hard-coded list of add-on toolboxes, we then check for the
% presence of each and issue a warning if absent.
V = ver;
VName = {V.Name};
warnState = warning();
warning off backtrace
for ii=1:length(requiredAddOns)
    if ~any(strcmp(VName, requiredAddOns{ii}))
        warnString = ['The Matlab ' requiredAddOns{ii} ' is missing. ' toolboxName ' may not function properly.'];
        warning('localHook:requiredMatlabToolboxCheck',warnString);
    end
end
warning(warnState);

end
