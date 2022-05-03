%% Example of computation of the observability of nonlinear systems,
%  using the Rossler system as a test case.

% Code by Arthur Montanari, 01/2022
% If you use this code, please reference:
% A. N. Montanari and L. A. Aguirre, “Observability of Network Systems: 
% A Critical Review of Recent Results,” Journal of Control, Automation and
% Electrical Systems, 31(6):1348–1374 (2020).
% DOI: 10.1007/s40313-020-00633-5.
% ResearchGate: https://bit.ly/39n08H8 (link to ResearchGate)

clear all; close all; clc;

n = 3;                               % system dimension
X = sym('x',[1 n],'real')';          % state variables (symbolic)
a = sym('a','Real');                 % parameters
b = sym('b','Real');
c = sym('c','Real');

% Nonlinear system function
f = [- X(2) - X(3);
     X(1) + a*X(2);
     b + X(3)*(X(1) - c)]

% Observation function
h = X(2)

% Observability matrix
Oc = obsvmatrix(X,f,h)
if rank(Oc) == n
    disp('System is locally observable.')
else
    disp('System is locally unobservable.')
end

% PS: Note that the Rossler system is always observable for h(x) = X2 since
% det(Oc) is not equal to zero; however it is not locally observable for 
% h(x) = X3 at state X3 = 0 since det(Oc) = 0 for X3 = 0. For more details,
% check the reference above.

% Note that the implementation above is based on symbolic computation of
% Lie derivatives and it was not designed to be efficient for
% high-dimensional systems.

% For a more computationally efficient implementation, I recommend
% checking the following references:
% [1] J. D. Stigter, L. G. van Willigenburg, and J. Molenaar, “An Efficient 
%     Method to Assess Local Controllability and Observability for 
%     Non-Linear Systems,” IFAC-PapersOnLine, 51(2):535–540 (2018).
% [2]﻿J. D. Stigter and D. Joubert, “Computing measures of identifiability,
%     observability, and controllability for a dynamic system model with 
%     the StrucID app,” IFAC-PapersOnLine, 54(7):138–143 (2021).

