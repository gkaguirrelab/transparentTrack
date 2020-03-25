function [x,lb,ub,lbp,ubp] = setBounds(x,bb,searchSet)

% Apply the bounds
lb = x - bb;
lbp = x - bb./2;
ubp = x + bb./2;
ub = x + bb;

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