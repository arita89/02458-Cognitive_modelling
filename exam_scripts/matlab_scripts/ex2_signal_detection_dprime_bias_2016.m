clc;
clear all;
close all;
dash = ('----------------------------------');
%% SIGNAL DETECTION 
% cat detect SOUND
tot_obs = 60; %number of tests
l =1; %number of test subjects
%----------------------
% NO sound played (Noise , N)
CR = 60-35
FA = 35
%----------------------
% YES sound played (SignalplusNoise , SN)
H = 57
M = 60-57

%% calculate probabilities

P_H = H/tot_obs
P_FA = FA/tot_obs
P_M = M/tot_obs
P_CR = CR/tot_obs
%disp (dash)

%% d_prime - Definition

d_prime = @(Prob_hits,Prob_FalseAlarm) norminv(Prob_hits)-norminv(Prob_FalseAlarm);

% the observer capability on distinguishing two different stimuli
% can be seen as the distance between the mean values of the two normal curves, N and SN
% the approximation of d? is not highly affected by the Criterion. 
% d? is an index of detectability of the input and it is defined as sensitivity


% d' = Z(P_H)-Z(P_FA) = Z_Noise - Z_SignalPlusNoise
% Z = norminv () STANDARDIZE

%d? and lambda are measured in standard deviations (??)

for i = 1:l
    %disp (i)
    dprime(i) = d_prime (P_H(i),P_FA(i));
    %disp (dash)
end

dprime;

%% criterion or BIAS
% how easy is the test subject to answer "yes i perceive the input"
% lax, lower then dprime/2 =  very easy to say yes
% conservative, higher = more difficult to say yes

% lambda = criterion
% phi = normpfd() <----------------- Gaussian PDF
% normcfd() <----------------- Gaussian cumulative 
% norminv is the inverse of normPDF <-----------------

% criterion = lambda = - (norminv(P_FA))
lambda = -(norminv(P_FA));
%%
disp(dash)
disp(['The sensitivity is ',num2str(dprime),'?'])
disp(['The bias is ',num2str(lambda),'?'])

%% Signal Detectoin Theory

% lambda = criterion
% phi = ? = normpfd() 
% y = normpdf(x,mu,sigma) returns the <----- RIGHT ONE
% pdf of the normal distribution 
% with mean mu and standard deviation sigma,
% evaluated at the values in x.
% norminv is the inverse of normpfd <-----------------

% another function
% p = normcdf(x) returns the cumulative distribution function (cdf) 
% of the standard normal distribution,
% evaluated at the values in x.

% P(CR) = phi(lambda)
% P(FA) = phi (-lambda)= 1-phi(lambda) 
% P(M) = phi ((lambda-mu_SN)/sigma_SN)
% P(H) = 1- phi ((lambda-mu_SN)/sigma_SN) 

%% the psychometric function is a cumulative Gaussian ?((x-x0)/?) 
% x is the sound level in db
% x0 is the 50% threshold  = CRITERION = lambda
% ? is the standard deviation = sound_level/dprime

x = 5; %sound_level
% x = sigma*dprime
sigma = x/dprime;
fiftythreshold = lambda*sigma;
%%
disp(dash)
disp(['The sound level is ',num2str(x),' dB'])
disp(['The standard deviation is ',num2str(sigma),' dB'])
disp(['The 50% of thethreshold is ',num2str(fiftythreshold ),' dB'])