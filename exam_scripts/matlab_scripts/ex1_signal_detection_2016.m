clc;
clear all;
close all;
dash = ('----------------------------------');
%% SIGNAL DETECTION 

% N = noise (mu_N = 0, sigma_N)
% CR =  correct rejection (true negative)
% FA = False alarm (false positive)
%----------------------
% SN = signal+ noise (mu_SN , sigma_SN)
% or simplified (mu_S , sigma_S)
% H = Hit (true positive)
% M = Miss (false negative)
%----------------------
% Hit + Miss = tot_obs
% Correct Regection + False alarm = total_observations
%----------------------
% P(H)+P(M) = 1
% P(CR) + P(FA) = 1

%% DATA
tot_obs = 60;
H = [41,39,32]
CR = [19,18,28]
%M = [];
%FA = [];

l = length(H) % number of subjects

%% initialization
% H = ones(1,l);
% CR = ones(1,l);
M = ones(1,l); % row vector (1 row, 3 columns)
FA= ones(1,l);
dprime = ones(1,l); %Initialize d_prime
disp (dash)

%% complementary data

for i = 1:l
    %i
    M(i) = tot_obs - H(i);%Miss  (false negative)
    FA(i) = tot_obs - CR(i); % False alarm, (false positive)
    %H(i) = tot_obs - M(i);%Miss  (false negative)
    %CR(i) = tot_obs - FA(i); % False alarm, (false positive)
    %disp (dash)
end

M
FA
H
CR
disp (dash)

%% calculate probabilities

P_H = H/tot_obs
P_FA = FA/tot_obs
P_M = M/tot_obs
P_CR = CR/tot_obs
disp (dash)

%% d_prime - Definition

d_prime = @(Prob_hits,Prob_FalseAlarm) norminv(Prob_hits)-norminv(Prob_FalseAlarm);

% the observer capability on distinguishing two different stimuli
% can be seen as the distance between the mean values of the two normal curves, N and SN
% the approximation of d? is not highly affected by the Criterion. 
% d? is an index of detectability of the input and it is defined as sensitivity


% d' = Z(P_H)-Z(P_FA) = Z_Noise - Z_SignalPlusNoise
% Z = norminv () STANDARDIZE

% y = normpdf(x,mu,sigma) returns the <----- RIGHT ONE
% pdf of the normal distribution 
% with mean mu and standard deviation sigma,
% evaluated at the values in x.

% p = normcdf(x) returns the cumulative distribution function (cdf) 
% of the standard normal distribution,
% evaluated at the values in x.

for i = 1:l
    disp (i)
    dprime(i) = d_prime (P_H(i),P_FA(i));
    disp (dash)
end

dprime

%% d_prime - Interpretation
% the observer capability on distinguishing two different stimuli
% can be seen as the distance between the mean values of the two normal curves, N and SN
% the approximation of d? is not highly affected by the Criterion. 
% d? is an index of detectability of the input and it is defined as sensitivity

% d_prime = 0 or <0 ----> stimuli cannot be distinguished
% the bigger the d_prime the more stimuli can be distinguished