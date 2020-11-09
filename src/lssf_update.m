% Update Step for LSSF
%
% Exact inference (locally). Laplace likelihood, Gaussian prior.
%
% Inputs:
% m: Gaussian predictive distribution mean
% P: Gaussian predictive distribution covariance matrix
% y: data observation
% C: output matrix
% R: scale of the Laplace noise
% 
% Outputs:
% mu: Posterior mean
% V: Posterior covariance
% logml: log-marginal likelihood, log p(y).
%
%
% Citation:
% J. Neri, P. Depalle and R. Badeau, "Laplace State Space Filter with 
% Exact Inference and Moment Matching," IEEE International 
% Conference on Acoustics, Speech and Signal Processing (ICASSP), 
% pp. 5880-5884, Barcelona, Spain, 2020. 
%
% Author: Julian Neri
% Affil: McGill University, Montreal, CA
% Date: May 1, 2020

function [mu, V, logml] = lssf_update(m,P,y,C,R)

    S = C*P*C';
    yt = y - C*m;

    e1 = (S/R - yt)/sqrt(2*S);
    e2 = (S/R + yt)/sqrt(2*S);
    
    % Scaled complementary error function erfcx().
    Phi_m = erfcx(e1);
    Phi_p = erfcx(e2);

    delta_y = (Phi_m-Phi_p)./(Phi_m+Phi_p);
    Delta_y = 1./(Phi_m+Phi_p).*(4*(Phi_m*Phi_p)./(Phi_m+Phi_p)-2*R/sqrt(pi/2*S));
    
    % Check for inf and set to analytic value at limit.
    if isinf(Phi_m)
        delta_y = 1;
        Delta_y = 0;
    end
    if isinf(Phi_p)
        delta_y = -1;
        Delta_y = 0;
    end
    if isnan(delta_y) || isinf(delta_y)
        if Phi_p > Phi_m
            delta_y = -1;
        else
            delta_y = 1;
        end
    end
    if isnan(Delta_y) || isinf(Delta_y)
        Delta_y = 0;
    end
    
    % Gain
    k = P*C'/R;
    
    % Posterior mean and covariance
    mu = m + k*delta_y;
    V  = P + (k*k')*Delta_y;
    
    % Compute the Log Marginal Likelihood
    phi_m = erfc(e1);
    phi_p = erfc(e2)*exp(2*yt/R);
    logml = -log(4*R) + (S-2*R*yt)/(2*R^2) + log(phi_m+phi_p);
 
end
