%% Functional observability and early-warning signal of the Epileptor.
% This code evaluates the complete and functional observability of the
% Epileptor model. It also computes the SVDO coefficient for indirect
% quantification of the system's observability from time-series data. The
% SVDO can be used for early-warning signal of seizure-like events in the 
% Epileptor model.

% Arthur Montanari, 06/2022
% If you use this code, please reference:
% A. N. Montanari, L. Freitas, D. Proverbio, J. Gonçalves, 
% “Functional observability and subspace reconstruction in nonlinear
% systems” (2022).

clear all; close all; clc;
addpath('..');

%% Epileptor model simulation

% Simulation time
dt = 1e-2;
t0 = 0;
tf = 10000;
T = length(0:dt:tf);                % number of data points
st = 0.1;                           % sampling time
N = length(0:round(st/dt):tf/dt);   % sampled data points

% Initial conditions
x(:,1) = [0.022;0.91;3.84;-1.11;0.73;0];
taux(1) = t0;

% Parameters
n = 6;                              % system dimension
x0 = -1.6;      y0 = 1;
Irest1 = 3.1;   Irest2 = 0.45;
tau0 = 2857;    tau2 = 10;
gamma = 0.01;

% Noise
nu = [0.01 0.1];                    % stochastic system
% nu = [0 0];                       % deterministic system

% Euler-Maruyama integration
for k = 2:T
    % Deterministic part
    fvec = epileptor(taux(k-1),x(:,k-1),x0,y0,Irest1,Irest2,tau0,tau2,gamma); 
    taux(k) = taux(k-1) + dt;
    
    % First subsystem (x1,y1)
    x(1,k) = x(1,k-1) + fvec(1) * dt + sqrt(nu(1))*sqrt(dt)*randn;
    x(2,k) = x(2,k-1) + fvec(2) * dt + sqrt(nu(1))*sqrt(dt)*randn;
    
    % Slow variable z
    x(3,k) = x(3,k-1) + fvec(3) * dt;
    
    % Second subsystem (x2,y2) (0.1/0.25/0.71/1)
    x(4,k) = x(4,k-1) + fvec(4) * dt + sqrt(nu(2))*sqrt(dt)*randn;
    x(5,k) = x(5,k-1) + fvec(5) * dt + sqrt(nu(2))*sqrt(dt)*randn;
    
    % Dummy variable g
    x(6,k) = x(6,k-1) + fvec(6) * dt;
end
x = x';

% Down-sampling
t = taux(1:round(st/dt):end);
s = x(1:round(st/dt):end,:);

% Measured time-series data
y = [s(:,1) + s(:,4)];

% Functional
z = s(:,3);

%% Plots of time-series data
figure(1)
subplot(121)
plot(t,y(:,1),'b',t,z,'r');
xlabel('t'); ylabel('y(t), z(t)')
subplot(122)
plot3(s(:,1),s(:,4),s(:,3));
xlabel('x_1'); ylabel('x_4'), zlabel('x_3');

%% Computes the observability matrix

% Symbolic functions
X = sym('x',[1 n],'real')';

% Nonlinear system function
f1 = [X(1)^3 - 3 * X(1)^2;
      (X(4) - 0.6 * (X(3) - 4)^2) * X(1)];
f2 = [0;
      6 * (X(4) + 0.25)];
for i = 1:2
    for j = 1:2
        f{i,j} = [X(2) - f1(i) - X(3) + Irest1;
                  y0 - 5 * X(1)^2 - X(2);
                  1/tau0 * (4 * ( X(1) - x0 ) - X(3));
                  - X(5) + X(4) - X(4)^3 + Irest2 + 2 * X(6) - 0.3 * (X(3) - 3.5);
                  (1/tau2) * (- X(5) + f2(j));
                  - gamma * (X(6) - 0.1 * X(1))];
    end
end

% Observation function
h = [X(1) + X(4)];

% Estimation function
g = X(3);

% Observability matrices corresponding to the discontinuous functions f1,f2
for i = 1:2
    for j = 1:2
        [Oc{i,j}, Lieh{i,j}] = nonlinearobsvmatrix(X,f{i,j},h,n-1);
    end
end

% Gradient of the functional g(x)
Dg = simplify( jacobian(g,X) );

%% Coefficient of (functional) observability
kappaPhi = zeros(N,1); 
kappaPsi = zeros(N,1);
for k = 1:N
    if mod(k,100) == 0; disp(['Counting ',num2str(k),'/',num2str(N)]); end
    
    % System trajectory
    xvec = s(k,:)';                      % state vector at time k 
    if xvec(1) < 0;     i = 1; else; i = 2; end
    if xvec(4) < -0.25; j = 1; else; j = 2; end 
    
    % Computes numerical values
    Oc_x = double(subs(Oc{i,j},X,xvec)); % obsv. matrix at time k
    Phi_x = inv(Oc_x);                   % map Phi at k
    Dg_x = double(subs(Dg,X,xvec));      % gradient of g(x) at time k
    Psi_x = Dg_x*Phi_x;                  % map Psi at time k
    
    % Computes Lie derivatives (differential embedding) at time k
    Lieh_x(k,:) = double(subs(Lieh{i,j},X,s(k,:)'));
    
    % Computes condition number
    normp = 2;
    kappaPhi(k) = norm(Phi_x,normp);
    kappaPsi(k) = norm(Psi_x,normp);
end

%% Coefficients of observability
disp([' '])
disp(['Measurement function:']); h
disp(['Functional:']); g
disp(['Coefficient of complete observability (average): ', num2str(mean(kappaPhi))])
disp(['Coefficient of functional observability (average): ', num2str(mean(kappaPsi))])

% Plot
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

%% Time-delay Embedding
tau = 20; 
dE = 5;
emb{1} = y((dE-1)*tau+1:end,1);
for i = 2:dE
    emb{i} = y((dE-i)*tau+1:end-(i-1)*tau,1)
end
y = [];
for i = 1:dE
    y(:,i) = [zeros((dE-1)*tau,1); emb{i}];
end

% Embedding space
trans = 1;
figure(3)
plot3(y(trans:end,1),y(trans:end,2),y(trans:end,3));
xlabel('y(t)'); ylabel('y(t-\tau)'); zlabel('y(t-2\tau)');

%% Early-warning signal

% Moving window
window = 50;          % window size
count = 1;
for i = window+1:window/2:N
    % Data window
    data = y(i-window:i,:);
    data = data./max(y);

    % Early warning signal
    [U,K,V] = svd(data(:,1:dE));
    sigma_dE(count) = K(dE,dE);

    % Time instant
    ewt(count) = t(i);
    count = count + 1;
end

% Plot
figure(4)
xaxis = [0e3 10e3];

subplot(311);
plot(t,y(:,1),'Color',[0 0.4470 0.7410]); hold on;
plot(t,s(:,3),'Color',[0.6350 0.0780 0.1840]);
xlim(xaxis);

subplot(312);
semilogy(t,kappaPsi,'Color',[0 0.4470 0.7410]);
xlim(xaxis);

subplot(313);
semilogy(ewt,sigma_dE,'Color',[0 0.4470 0.7410]);
xlim(xaxis);
