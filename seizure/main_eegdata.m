%% Early-warning signal of seizure using EEG data.
% This code uses the SVDO coefficient as an early-warning signal of seizure
% events using human intracranial EEG data.

% Arthur Montanari, 06/2022
% If you use this code, please reference:
% A. N. Montanari, L. Freitas, D. Proverbio, J. Gonçalves, 
% “Functional observability and subspace reconstruction in nonlinear
% systems” (2022).

% To use this code, first download the public data available at 
% P. Karoly, M. Cook, L. Kuhlmann, D. Freestone, D. Grayden, E. Nurse, 
% A. Lai, D. Payne, W. D’Souza, U. Seneviratne, S. Berkovic, T. O’Brien, 
% B. Litt, D. Himes, K. Leyde, D. Soudry, S. Ahmadizadeh, M. Maturana, 
% K. Dell, Melbourne NeuroVista Seizure Prediction Trial, University of 
% Melbourne. Dataset 10.26188/5b6a999fa2316 (2018)

clear all; close all; clc;

%% Loads EEG data
load Seizure_010.mat

% Parameters
T = length(data);
sf = 400;                   % sampling frequency (Hz)
st = 1/sf;                  % sampling time
t = 0:st:T*st-st;           % time
channel = 2;                % selects channel

figure(1)
subplot(121);
plot(t,data(:,channel));
xlabel('t'); ylabel('EEG')

%% Time-delay embedding
y = data(:,channel);
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
figure(1)
subplot(122);
plot3(y(:,1),y(:,2),y(:,3));
xlabel('y(t)'); ylabel('y(t-\tau)'); zlabel('y(t-2\tau)');

%% Early-warning signal

% Moving window
window = 400;          % window size
count = 1;
for i = 1:window/4:T-window
    % Data window
    datawindow = y(i:i+window-1,:);
    
    % Early warning signal
    [U,K,V] = svd(datawindow(:,1:dE));
    sigma_dE(count) = K(dE,dE);
    
    % Time instant
    sigma_dE(count) = (norm(K,2));         % early-warning signal
    ewt(count) = t(i+window);          % time horizon
    count = count + 1;
end

% Plot EWS
figure(2);
subplot(211); 
plot(t,y(:,1),'Color',[0 0.4470 0.7410]); 
xlim([0 max(t)])
subplot(212); 
plot(ewt,sigma_dE,'Color',[0.8500 0.3250 0.0980]); hold on;
plot([71 71],[0 4000],'--g');
xlim([0 max(t)])
