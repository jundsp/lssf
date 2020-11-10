% Prediction Step for LSSF
%
% Identical to the Kalman filter's update, because state transition and
% previous state's posterior are Gaussian.
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

function [m,P] = lssf_predict(mu,V,A,Q)
    m = A*mu;
    P = A*V*A' + Q;
end