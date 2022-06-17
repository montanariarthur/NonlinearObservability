%% Complete and functional observability of nonlinear systems.
% This code evaluates the complete and functional observability of several
% dynamical nonlinear systems f(x), considering a measurement function h(x)
% and a functional g(x) sought to be reocnstructed.
% Calculations are based on symbolic computation using Lie derivatives.

% Arthur Montanari, 06/2022
% If you use this code, please reference:
% A. N. Montanari, L. Freitas, D. Proverbio, J. Gonçalves, 
% “Functional observability and subspace reconstruction in nonlinear
% systems” (2022).

clear all; close all; clc;
addpath([pwd,'/Systems/'])

% The following nonlinear systems were implemented:
% - rossler
% - lorenz63
% - cord
% - neuronHR
nonlinearsys = 'lorenz63'

% Simulation time
dt = 1e-2;
t0 = 0; 
tf = 1010;
trans = 1000;            % transiente time
st = 1e-2;              % sampling time
N = length(trans/dt:round(st/dt):tf/dt);      % data points

%% Nonlinear system
if contains(nonlinearsys, 'rossler')
    %% ==================== ROSSLER SYSTEM ==================== %%
    % System dimension
    n = 3;
    
    % Parameters (spiral a = 0.398, screw a = 0.450)
    a = 0.398; b = 2; c = 4;
    
    % Numerical simulation
    s0 = randn(1,n);
    [taux,saux] = odeRK(@(t,s)rossler(t,s,a,b,c),[t0 dt tf],s0);
    
    % Remove transient, samples data
    t = taux(trans/dt:round(st/dt):end);
    s = saux(trans/dt:round(st/dt):end,:);
    
    % Symbolic functions
    X = sym('x',[1 n],'real')';
    
    % Nonlinear system function
    f = [- X(2) - X(3);
        X(1) + a*X(2);
        b + X(3)*(X(1)-c)];
    
    % Observation function
    h = X(3);
    
    % Estimation function
    g = atan(X(2)/X(1));
    
elseif contains(nonlinearsys, 'lorenz63')
    %% ==================== LORENZ SYSTEM ==================== %%
    % System dimension
    n = 3;
    
    % Parameters
    sigma = 10; rho = 28; beta = 8/3;
    
    % Numerical simulation
    s0 = 0.1*randn(1,n);
    [taux,saux] = odeRK(@(t,s)lorenz(t,s,sigma,rho,beta),[t0 dt tf],s0);
    
    % Remove transient, samples data
    t = taux(trans/dt:round(st/dt):tf/dt);
    s = saux(trans/dt:round(st/dt):tf/dt,:);
    
    % Symbolic functions
    X = sym('x',[1 n],'real')';
    
    % Nonlinear system function
    f = [sigma*(X(2) - X(1));
        rho*X(1) - X(2) - X(1)*X(3);
        X(1)*X(2) - beta*X(3)];
    
    % Observation function
    h = X(1);
    
    % Estimation function
    g = X(2);
    
elseif contains(nonlinearsys, 'cord')
    %% ==================== CORD SYSTEM ==================== %%
    % System dimension
    n = 3;
    
    % Parameters
    a = 0.258; b = 4.033; Fp = 8; G = 1; order = 1;
    
    % Numerical simulation
    s0 = 0.1*randn(1,n);
    [taux,saux] = odeRK(@(t,s)cord(t,s,a,b,Fp,G,order),[t0 dt tf],s0);
    
    % Remove transient, samples data
    t = taux(trans/dt:round(st/dt):tf/dt);
    s = saux(trans/dt:round(st/dt):tf/dt,:);
    
    % Symbolic functions
    X = sym('x',[1 n],'real')';
    
    % Nonlinear system function
    f = [-X(2)^order - X(3)^order - a*X(1) + a*Fp;
        X(1)*X(2) - b*X(1)*X(3) - X(2) + G;
        b*X(1)*X(2) + X(1)*X(3) - X(3)];

    % Observation function
    h = X(2);
    
    % Estimation function
    g = atan(X(3)/X(2));

elseif contains(nonlinearsys, 'neuronHR')
    %% ==================== HINDMARSH-ROSE NEURON ==================== %%
    % System dimension
    n = 3;
    
    % Parameters
    a = 1; b = 3; c = 1; d = 5;
    I = 3.25; r = 0.001; sdyn = 4; xR = -(1+sqrt(5))/2;%-8/5; 
    
    % Numerical simulation
    s0 = 0.3 + (0.7-0.3).*rand(1,3); % s0 = randn(3,1);
    [taux,saux] = odeRK(@(t,s)neuronHR(t,s,a,b,c,d,I,r,sdyn,xR),[t0 dt tf],s0);
    
    % Remove transient, samples data
    t = taux(trans/dt:round(st/dt):tf/dt);
    s = saux(trans/dt:round(st/dt):tf/dt,:);
    
    % Symbolic functions
    X = sym('x',[1 n],'real')';
    
    % Nonlinear system function
    f = [X(2) - a*X(1)^3 + b*X(1)^2 + I - X(3);
        c - d*X(1)^2 - X(2);
        r*(sdyn*(X(1)-xR) - X(3))];
    
    % Observation function 
    h = X(2);
    
    % Estimation function
    g = X(1); 

end

%% Observability matrices

% Computes observability matrix
[Oc, Lieh] = nonlinearobsvmatrix(X,f,h,n-1);
Dg = simplify( jacobian(g,X) );

% Coefficient of (functional) observability
kappaPhi = zeros(N,1);
kappaPsi = zeros(N,1);
Lieh_x = zeros(N,n);
for k = 1:N
    if mod(k,100) == 0; disp(['Counting ',num2str(k),'/',num2str(N)]); end
    
    % Computes numerical values at time k
    xvec = s(k,:)';                   % state vector at time k
    Oc_x = double(subs(Oc,X,xvec));   % obsv. matrix at time k
    Phi_x = inv(Oc_x);                % map Phi at k
    Dg_x = double(subs(Dg,X,xvec));   % gradient of g(x) at time k
    Psi_x = Dg_x*Phi_x;               % map Psi at time k
    
    % Computes Lie derivatives (differential embedding) at time k
    Lieh_x(k,:) = double(subs(Lieh,X,s(k,:)'));
    
    % Computes condition number
    normp = 2;
    kappaPhi(k) = norm(Phi_x,normp);
    kappaPsi(k) = norm(Psi_x,normp);
end

%% Coefficients of observability
disp([' '])
disp(['Nonlinear system: ', nonlinearsys])
disp(['Measurement function:']); h
disp(['Functional:']); g
disp(['Observability matrix:']); Oc
disp(['Coefficient of complete observability (average): ', num2str(mean(kappaPhi))])
disp(['Coefficient of functional observability (average): ', num2str(mean(kappaPsi))])

%% Plots

% Original state space and embedding space
figure(1);
subplot(121); plot3(s(:,1),s(:,2),s(:,3));
xlabel('x_1'); ylabel('x_2'), zlabel('x_3');
subplot(122); plot3(Lieh_x(:,1),Lieh_x(:,2),Lieh_x(:,3));
xlabel('y'); ylabel('ydot'), zlabel('yddot');

% Coefficients of observability
figure(2)
sample = 1; dotsize = 20;
colormap((pink))

subplot(121); scatter3(s(1:sample:end,1),s(1:sample:end,2),s(1:sample:end,3),...
    dotsize,kappaPhi(1:sample:end),'o','filled')
colorbar; set(gca,'ColorScale','log'); 
caxis([min([kappaPhi; kappaPsi]) max([kappaPhi; kappaPsi])]);
xlabel('x_1'); ylabel('x_2'), zlabel('x_3');

subplot(122); scatter3(s(1:sample:end,1),s(1:sample:end,2),s(1:sample:end,3),...
    dotsize,kappaPsi(1:sample:end),'o','filled')
colorbar; set(gca,'ColorScale','log'); 
caxis([min([kappaPhi; kappaPsi]) max([kappaPhi; kappaPsi])]);
xlabel('x_1'); ylabel('x_2'), zlabel('x_3');
