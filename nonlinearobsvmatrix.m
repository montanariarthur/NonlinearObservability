function [Oc, Lieh] = nonlinearobsvmatrix(X,f,h,s)
% This code computes the nonlinear observability matrix of a nonlinear 
% system f(x) with measurement function h(x). Calculations are based on 
% symbolic computation using Lie derivatives.

% Arthur Montanari, 01/2022
% If you use this code, please reference:
% A. N. Montanari and L. A. Aguirre, “Observability of Network Systems: 
% A Critical Review of Recent Results,” Journal of Control, Automation and
% Electrical Systems, 31(6):1348–1374 (2020).
% DOI: 10.1007/s40313-020-00633-5.
% ResearchGate: https://bit.ly/39n08H8 (link to ResearchGate)

% Inputs:
%   X     -    variables (symbolic)
%   f     -    nonlinear system (symbolic)
%   h     -    measurement function (symbolic)
%   s     -    maximum Lie derivative (default is n-1)

% Outputs:
%   Oc    -    observability matrix (symbolic)
%   Lieh  -    Lie derivatives (symbolic)


% s = n - 1 is the default input
if nargin < 4
    n = size(f,1);
    s = n - 1;
end

% Lie derivatives
Lie_h{1} = h;
dLie_h{1} = simplify( jacobian(Lie_h{1},X) );
for i = 1:s
    Lie_h{i+1} = simplify( dLie_h{i}*f );
    dLie_h{i+1} = simplify( jacobian(Lie_h{i+1},X) );
end
    
% Observability matrix
Lieh = vertcat(Lie_h{1});
Oc = vertcat(dLie_h{1});
for i = 1:s
    Oc = vertcat(Oc,dLie_h{i+1});
    Lieh = vertcat(Lieh,Lie_h{i+1});
end
Oc = simplify(Oc);

end