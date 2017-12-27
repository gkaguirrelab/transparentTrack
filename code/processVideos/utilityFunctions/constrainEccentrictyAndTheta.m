function [c, ceq] = constrainEccentrictyAndTheta( candidateEllipse, targetEllipse )


% First constraint - C reflects the mismatch in theta
c = abs(candidateEllipse(5) - targetEllipse(5));

% Second constraint - ceq reflects the mismatch in eccentricity
ceq = candidateEllipse(4) - targetEllipse(4);


end

