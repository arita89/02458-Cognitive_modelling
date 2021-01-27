clc;
clear all;
close all;
format long;
dash = ('----------------------------------');
%% SIGNAL DETECTION for UNEQUAL VARIANCE BTW NOISE AND SIGNAL

% when we have rate responses 
% we have different ORDERED n categories
% categories are separated by n-1 criterions (lambda1,lambda2..etc)
% we want to calculate for each criterion: P_FA and P_H
% ?(-1)(P_H) and ?(-1)(P_FA) fit a straight line

n = 6; %number categories, == len(N)== len(s)
n_lambdas = n-1;
t_n= 100; % number of NOISE inputs 

t_sn= 100;% number of SIGNAL+noise inputs
tot_obs = 100;%for simplicity, since its the same for Noise and Signal

N = [2,31,57,8,2,0]; %Noise
N_CR = [2,31,57];
N_FA = [8,2,0]; 
%------------------------------
S = [2,14,27,38,15,4]; %Signal(+noise)
S_M = [2,14,27];
S_H = [38,15,4];
%%
No_responses_N = sum(N(1:3))
No_responses_S = sum(N(4:6))
Yes_responses_N = sum(S(1:3))
Yes_responses_S = sum(S(4:6))

Tot_NO = [No_responses_N,No_responses_S]
Tot_YES = [Yes_responses_N,Yes_responses_S]

%% 

%two_by_two_table(:,:,n_lambdas)= zeros(2, 2, n_lambdas)
for i = 1:n_lambdas
   disp (['criterion ',num2str(i)])
   N_NO_per_lambda =  sum(N(1:i))
   S_NO_per_lambda =  sum(S(1:i))
   N_YES_per_lambda =  sum(N(i:n_lambdas))
   S_YES_per_lambda =  sum(N(i:n_lambdas))
   %two_by_two_table(:,:,i) = [N_NO_per_lambda,N_YES_per_lambda;S_NO_per_lambda,S_YES_per_lambda]
   disp (dash)
end

for i = 1:n_lambdas
    disp (['criterion ',num2str(i)])
   %two_by_two_table(i) 
end

%% cumulatives

Cum_N = cumsum(N)
Cum_S = cumsum(S)
disp(dash)

%% initialization
% P_H = zeros(1,n); %row vectors
% P_FA = zeros(1,n);
%%
for i= 1:n-1 % 1,2,3,4,5,n
    P_FA(i) = (tot_obs-Cum_N(i))/tot_obs;
    P_H(i) = (tot_obs-Cum_S(i))/tot_obs;
end

P_FA 
P_H 
%% standardize
Z_FA = norminv(P_FA)
Z_H = norminv(P_H)


%%	ROC RECECIVER OPERATING PARAMETERS
% Z_FA (x axis IN ROC)
% Z_H (y axis IN ROC)
% mu_S
% sigma_S

%% in Gaussian coordinates ROC characteristics is a straight line
% ?-1(P_H)=(1/?_S)*(?-1(P_FA))+(?_S/?-S)
% Z_H = (1/?_S)*Z_FA +(?_S/?-S)
% y = mx +q

% y = P(1)x +P(2)<---- using polyfit
% m= P(1)= (1/?_S) <---- SLOPE
% Q = P(2)=(?_S/?_S) <---- INTERCEPT

%% linear fit to data
c = 4;
x = Z_FA(1:c);
y = Z_H(1:c);
% 
% % fitlm <---- WRONG
% tbl = table(x,y)
% model = fitlm(tbl,'linear')
% anova(model,'summary')
% 
% A = model.Coefficients.Estimate(2:end)'; %<--------- take the transpose to have row vector
% m = A(1)%slope
% q = model.Coefficients.Estimate(1) % intercept
% s_S = 1/m
% m_S = q*s_S

% polyfit
%p = normcdf(x,mu,sigma) returns the cdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
P = polyfit(x,y,1); 
yfit = P(1)*x+P(2);
sigma_S = 1/P(1)
mean_S = P(2)*sigma_S

% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(P(1)) '*x + ' num2str(P(2))])
 
%% plot in Gaussian coordinates the ROC
% receiver operating caratheristic ROC curve  z_H over z_FA
% z_H = (1/sigma_signal)* z_FA + (mu_signal/sigma_signal)

figure(1)
hold on
x = Z_FA;
y = Z_H;
scatter(x,y,25,'b','o') 

%% d_prime - Definition

d_prime = @(Prob_hits,Prob_FalseAlarm) norminv(Prob_hits)-norminv(Prob_FalseAlarm);

%%
% how to check if data comes from equal variance or unequal variance? 
% estimate d_prime for all the criterions
% if d_prime = constant --> equal variance
% if d_prime = NOT constant --> UNequal variance
disp(dash)
for i= 1:length(P_FA)
    dprime_approx(i) = d_prime(P_H(i),P_FA(i)); 
end

dprime_approx

% D_prime changes a lot with the criteriom == unequal variance