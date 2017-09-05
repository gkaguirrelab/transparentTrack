function deriveTimebaseFromLTData(glintFileName,ltReportFilename,timebaseFilename,varargin)

% deriveTimebaseFromLTData(glintFileName,ltReportFilename,timebaseFilename

% header

%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('glintFileName',@ischar);
p.addRequired('ltReportFilename',@ischar);
p.addRequired('timebaseFilename',@ischar);

% Optional analysis parameters


% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(glintFileName,ltReportFilename,timebaseFilename, varargin{:})
