function [x,lb,ub,lbp,ubp] = setBounds(x,model,stage,searchStrategy)


% Construct the param search set for this strategy and stage
searchSet = zeros(1,model.nParams);
sets = model.stages.(searchStrategy){stage};
for ii = 1:length(sets)
   fields = split(sets{ii},'.');
   idx = model.func.fieldSetIdx(fields{1},fields{2});
   idx = model.(fields{1}).idxMultiScene(idx);
   searchSet(idx) = 1;
end

% Apply the bounds
lb = x - model.bounds;
lbp = x - model.bounds./2;
ubp = x + model.bounds./2;
ub = x + model.bounds;

% Lock params
lb(~searchSet) = x(~searchSet);
ub(~searchSet) = x(~searchSet);
lbp(~searchSet) = x(~searchSet);
ubp(~searchSet) = x(~searchSet);

% Keep the corneal axis (torsion) between +-90
lb(3) = max([lb(3) -90]);
ub(3) = min([ub(3) 90]);
lbp(3) = max([lbp(3) -90]);
ubp(3) = min([ubp(3) 90]);

% Ensure that x is within bounds
x=min([ub; x]);
x=max([lb; x]);

end