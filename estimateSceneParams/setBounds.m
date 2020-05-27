function [x,lb,ub,lbp,ubp] = setBounds(x,model,stage,strategy)
% Provides the bounds and adjusts x as necessary for estimateSceneParams
%
% Syntax:
%  [x,lb,ub,lbp,ubp] = setBounds(x,model,stage,strategy)
%
% Description
%   Uses the values in the model structure for the specified strategy and
%   stage to create upper and lower bounds for the BADS search.
%
% Inputs:
%   x                     - 1xp vector of parameters.
%   model                 - Struct. Describes parameters of the search.
%   stage                 - Scalar. Which stage of the search we are
%                           currently conducting.
%   strategy              - Char vector. Directs the function to the field
%                              model.strategy.(strategy)
%
% Outputs:
%   x                     - 1xp vector of parameters, potentially adjusted
%                           to be within the upper and lower bounds.
%   lb, ub                - 1xp vectors that specify the hard upper and
%                           lower bounds.
%   lbp, ubp              - 1xp vectors that specify the probable bounds
%                           within which BADS will concentrate its search.
%


%% Touch the x vector
% We add 1e-10 to the corneal torsion x value so that the x vector has
% changed, which will force the objective function to recompute.
corneaTorsionIdx = find(model.eye.idxMultiScene(model.func.fieldParamIdx('eye','torsion')));
x(corneaTorsionIdx) = x(corneaTorsionIdx)+1e-10;


%% Construct the param search set for this strategy and stage

% A blank search set
searchSet = zeros(1,model.nParams);

% The sets of fields and parameters to be processed, identified by indirect
% reference to the model structure using the strategy variable
sets = model.strategy.(strategy).stages{stage};

% Loop through the number of parameter sets that define this search
for ii = 1:length(sets)
    
    % Each set entry specifies the field (e.g., "head" or "eye") and the
    % set name (e.g., "cameraPosition") in the form "field.set"
    fields = split(sets{ii},'.');
    
    % The function fieldSetIdx maps from the a parameter index for a field
    % to the entire x parameter vector for a given scene
    idx = model.func.fieldSetIdx(fields{1},fields{2});
    
    % The idxMultiScene function maps the x vector for one scene to a full,
    % concatenated vector of paramaters across all scenes in the search.
    idx = model.(fields{1}).idxMultiScene(idx);
    
    % Indicate that we will search the identified parameters
    searchSet(idx) = 1;
end

% Apply the bounds. The plausible bounds are simply half the full bounds.
lb = x - model.bounds;
lbp = x - model.bounds./2;
ubp = x + model.bounds./2;
ub = x + model.bounds;

% Lock params
lb(~searchSet) = x(~searchSet);
ub(~searchSet) = x(~searchSet);
lbp(~searchSet) = x(~searchSet);
ubp(~searchSet) = x(~searchSet);

% A special case issue is that the axis of corneal torsion must be locked
% between ±90 degrees.
lb(corneaTorsionIdx) = max([lb(corneaTorsionIdx) -90]);
ub(corneaTorsionIdx) = min([ub(corneaTorsionIdx) 90]);
lbp(corneaTorsionIdx) = max([lbp(corneaTorsionIdx) -90]);
ubp(corneaTorsionIdx) = min([ubp(corneaTorsionIdx) 90]);

% Ensure that x is within bounds
x=min([ub; x]);
x=max([lb; x]);

end