clc, clear all, close all

%% Exercise 3.3 Bayesian Multisensory Integration
%import data

x_load_a  = importdata('xA.txt');
x_load_v  = importdata('xV.txt');
x_load_av = importdata('xAV.txt');

%%
test_n = 2;
training_data = [x_load_a x_load_v x_load_av];
n_sample = size(training_data,1);

%% Strong fusion models - Function defintion

sigma_av = @(s_a,s_v) (s_a^2 * s_v^2)/ (s_a^2 + s_v^2);
w_av = @(s_a,s_v) (s_v^2)/ (s_a^2 + s_v^2);
mu_av = @(mu_a,mu_v,s_a,s_v) w_av(s_a,s_v)*mu_a + (1-w_av(s_a,s_v))*mu_v;

log_likelihood_solo = @(x,mu,st) normlike([mu,st],x);
log_likelihood_sum = @(mu_a,mu_v,st_a,st_v) log_likelihood_solo(training_data(:,1),mu_a,st_a) + log_likelihood_solo(training_data(:,2),mu_v,st_v) + log_likelihood_solo(training_data(:,3),mu_av(mu_a,mu_v,st_a,st_v),sigma_av(st_a,st_v)); 

x_hat_0 = [0,0,1,1];
x_hat_av = fminunc(@(x) log_likelihood_sum(x(1),x(2),x(3),x(4)),x_hat_0)
m_av_hat = mu_av(x_hat_av(1), x_hat_av(2), x_hat_av(3), x_hat_av(4))
sigma_av_hat = sigma_av(x_hat_av(3), x_hat_av(4))

%%

