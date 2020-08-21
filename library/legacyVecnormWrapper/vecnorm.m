function n = vecnorm(A, p, dim)


if(nargin < 2)
    p = 2; dim = find(size(A)>1,1);
elseif (nargin < 3)
    dim = find(size(A)>1,1);
end
A2p = A.^p;
s   = sum(A2p, dim);
n   = s.^(1/p);

% If we are running 2017b or later, remove this function and use the built in.
if ~verLessThan('matlab','9.3')
    % The directory that contains this function
    enclosingDir = fileparts(mfilename('fullpath'));
    rmpath(enclosingDir)
end


end