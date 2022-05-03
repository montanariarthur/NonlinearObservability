%% Example of computation of the identifiability of nonlinear DAE systems,
%  using the pendulum equation as a test case.

% Code by Arthur Montanari, 01/2022
% If you use this code, please reference:
% A. N. Montanari and L. A. Aguirre, “Observability of Network Systems: 
% A Critical Review of Recent Results,” Journal of Control, Automation and
% Electrical Systems, 31(6):1348–1374 (2020).
% DOI: 10.1007/s40313-020-00633-5.
% ResearchGate: https://bit.ly/39n08H8 (link to ResearchGate)

clear all; close all; clc;

% Symbolic variables
n = 3;                               % system dimension
X = sym('x',[1 n],'real')';          % state variables (symbolic)
g = sym('g','real');                 % parameters
m = sym('m','real');
L = sym('L','real'); 

% Parameters to be identified
theta = symvar([m L g])'; 
np = length(theta);                  % number of parameters
N = n + np;                          % extended system size

% Nonlinear system function
f = [X(3);
     -3*g/(L^2)*(-X(1)*X(3)/(sqrt(L^2-X(1)^2)));
     -X(2)*X(1)/m;                   % pendulum equations (ODE)
     zeros(np,1)];                   % parameters treated as state vars

% Observation function
h = atan(-X(1)/X(2));                % output is the pendulum angle

% Observability matrix
Oc = obsvmatrix([X; theta],f,h)
if rank(Oc) == N
    disp('System is locally identifiable.')
else
    disp('System is locally unidentifiable.')
end

% PS: Note that the pendulum is identifiable for theta = [m, L], but not
% for theta = [m, L, g].
