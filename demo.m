%%% Demo Script for Laplace State Space Filter (LSSF) %%%
%
% Citation:
% J. Neri, P. Depalle and R. Badeau, "Laplace State Space Filter with 
% Exact Inference and Moment Matching," IEEE International 
% Conference on Acoustics, Speech and Signal Processing (ICASSP), 
% pp. 5880-5884, Barcelona, Spain, 2020. 
%
% Author: Julian Neri
% Affil: McGill University
% Date: May 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;
addpath('data','src','utils');

%% Load Data
% Choose from: 'laplace noise', 'outliers', 'noise switch'.
trial = 'outliers';
[y,x_tru,~,~,A,~,C,~] = data_loader(trial);

%% Laplace State Space Filtering

% Data noise scale R (Laplace) and state noise covariance Q (Guassian)
R = .1; 
Q = 1e-3*eye(2);
% Isotropic Gaussian prior over initial state p(x_1).
m0 = zeros(2,1); P0 = eye(2);

% Infer
[mu,V,y_hat,S_hat,ell] = lssf(y,m0,P0,A,Q,C,R);

%% Plot Results

% Filtered Observations
t = 1:length(y_hat);
figure('pos',[0 1000 500 250]); hold on;
scatter(t,y,'filled','DisplayName','Data y'); 
plot_uncertainty(t,y_hat,S_hat);
xlabel('Time (n)'); ylabel('y'); box on; 
legend('Location', 'northoutside','Orientation','horizontal'); 

% Inferred Latent State Sequence
figure('pos',[0 300 500 300]);
subplot(211); hold on; 
plot(t,x_tru(1,:),'linewidth',2,'DisplayName','State x (true)'); 
plot_uncertainty(t,mu(1,:),squeeze(V(1,1,:)));
xlabel('Time (n)'); ylabel('dim. 1'); box on; axis tight;
legend('Location', 'northoutside','Orientation','horizontal'); 
subplot(212); hold on; 
plot(t,x_tru(2,:),'linewidth',2,'DisplayName','State x (true)'); 
plot_uncertainty(t,mu(2,:),squeeze(V(2,2,:)));
xlabel('Time (n)'); ylabel('dim. 2'); box on; axis tight; 
legend('Location', 'northoutside','Orientation','horizontal');

% Log-Marginal Likelihood ln p(y_n | y_{1:n-1})
figure('pos',[0 0 500 200]);
plot(ell,'k','linewidth',1.5); axis tight; xlabel('Time'); 
title('Log-Marginal Likelihood');