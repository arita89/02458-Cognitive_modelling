clear
clc
format long;  
dash = ('----------------------------------');
%%
x = 5; %intensity of the stimulus
d_prime = 1.1; 
d_prime/2

%% examples of possible critereons 
criterion = [-0.5 0.55 1.5];

%% the psychometric function is a cumulative Gaussian  
% x is the stimulus intensity
% x0 is the 50% threshold  = CRITERION = lambda
% ? is the standard deviation = sound_level/dprime
%x = sigma*dprime
sigma = x/d_prime;
fiftythreshold = criterion*sigma

%%
% sensitivity d? = ?-1 (P_H) - ?-1 (P_FA)
% d_prime= Z_H-Z_FA;

% criterion or bias = ?-1(P_CR)
% P_CR = 1- P_FA


% P_FA = normcdf(Z_FA)
% P_H = normcdf(Z_H)
% P_M = 1-P_H
% P_CR = 1-P_FA
% 
% format long
% 
% disp(' Intensity               P(Hit)             P(Miss)   ')
% disp ([ x' P_H' (1-P_H)'])
%%
disp(dash)
disp(['The sound level is ',num2str(x),' cd/m^3'])
disp(['The standard deviation is ',num2str(sigma),' cd/m^3'])
disp(['The 50% of the threshold is ',num2str(fiftythreshold ),' cd/m^3'])