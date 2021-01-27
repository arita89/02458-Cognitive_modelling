clc;
clear all;
close all;
format long;
dash = ('----------------------------------');
%% observer responds according to equal variance 
%signal detection model


d_prime = 1; %<-------USED FOR WHAT? mu_S  only?
lambda = 0.7;

%% Bayesian observer
%--------------Bayes--------------------
% P(A given B)P(B) = P(B given A)* P(A) 

% P(A given B) = posterior
% P (B given A) = Likelihood (not probability!)
% P (A) = prior

%----------Formula-------------------------
%P(A given B) "almost=" P(B given A)* P(A) 
% Posterior "almost=" Likelihood * Prior

% lambda = criterion
% phi = normPDF %Normal probability density function
% norminv is the inverse of normPDF <-----------------

% Likelihood = phi (lambda, mu, sigma)

% y = normpdf(x,mu,sigma) returns the <----- RIGHT ONE
% pdf of the normal distribution 
% with mean mu and standard deviation sigma,
% evaluated at the values in x.

% p = normcdf(x,mu,sigma) returns the cdf of the normal distribution 
% with mean mu and standard deviation sigma,
% evaluated at the values in x.

x = lambda;
mu_N = 0;
sigma_N = 1; %why?

mu_S=1;
sigma_S=1;%why?

L_X_N = normpdf(x,mu_N,sigma_N)
L_X_S = normpdf(x,mu_S,sigma_S)

Posterior_S = L_X_N/(L_X_N+L_X_S)