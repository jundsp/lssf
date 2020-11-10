% Laplace State Space Filter (LSSF)
% with exact inference and moment matching
% 
% Inputs:
% y: data sequence (1 X N)
% m0: Gaussian initial state prior mean (D X 1)
% P0: Gaussian initial state prior covariance (D X D)
% A: system dynamics matrix (D X D)
% Q: state noise covariance matrix (D X D)
% C: output matrix (1 X D)
% R: scale of the Laplace noise (scalar)
% 
% Outputs:
% mu: Posterior means
% V: Posterior covariances
% y_hat: marginal likelihood mean
% S_hat: marginal likelihood variance
% logml: log-marginal likelihood, log p(y).
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

function [mu, V, y_hat, S_hat, logml] = lssf(y,m0,P0,A,Q,C,R)

    [~,N] = size(y);
    dimx = size(m0,1);
    
    mu = zeros(dimx,N);
    V = zeros(dimx,dimx,N);
    y_hat = zeros(1,N);
    S_hat = zeros(1,N);
    logml = zeros(1,N);

    m = m0;
    P = P0;
    for n = 1:N
        % Update
        [mu(:,n), V(:,:,n), logml(n)] = lssf_update(m,P,y(:,n),C,R);
        % Predict
        [m, P] = lssf_predict(mu(:,n),V(:,:,n),A,Q);
        % Expected marginal likelihood mean and variance of the data
        y_hat(:,n) = C*mu(:,n);
        S_hat(:,n) = C*V(:,:,n)*C' + 2*R^2;
    end
end

