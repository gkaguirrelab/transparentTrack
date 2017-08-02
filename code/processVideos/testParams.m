function testParams( pathParams, varargin )

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

p.addParameter('stageToTest', 'glint', @ischar);

% parse
p.parse(pathParams, varargin{:})
pathParams=p.Results.pathParams;

% get the list of available gray files
pathParams.dataOutputDirFull = fullfile(pathParams.dataOutputDirRoot, pathParams.projectSubfolder, ...
    pathParams.subjectID, pathParams.sessionDate, pathParams.eyeTrackingDir);

FileList=dir(fullfile(pathParams.dataOutputDirFull, '*gray.avi'));

% Ask the user which file to work on 
fprintf('\nSelect a gray file:\n')
for pp=1:length(FileList)
    optionName=['\t' char(pp+96) '. ' FileList(pp).name '\n'];
    fprintf(optionName);
end

choice = input('\nYour choice (return for none): ','s');
if ~isempty(choice)
    choice = int32(choice);
    if choice >= 97 && choice <= 122
        grayVideoName=fullfile(pathParams.dataOutputDirFull,FileList(choice-96).name);
    end
end

switch p.Results.stageToTest
    case glint
        trackGlint(grayVideoName, '', 'displayMode', true, varargin{:});
    case pupil
        findPupilPerimeter(grayVideoName, '', 'displayMode', true, varargin{:});
end

end % function

