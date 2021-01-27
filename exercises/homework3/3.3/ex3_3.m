clc, clear all, close all

%% Exercise 3.3 Bayesian Multisensory Integration
%import data first part of exercise

x_a  = importdata('xA.txt');
x_v  = importdata('xV.txt');
x_av = importdata('xAV.txt');

%import data second part of exercise
ex2 = 0; %<------------- SWITCH
if ex2 == 1
    x_a  = importdata('xAlapse.txt');
    x_v  = importdata('xVlapse.txt');
    x_av = importdata('xAVlapse.txt');
end
%%
test_n = 2;
training_data = [x_a x_v x_av];
n_sample = size(training_data,1);

%% Strong fusion models - Function defintion

sigma_av = @(s_a,s_v) (s_a^2 * s_v^2)/ (s_a^2 + s_v^2);
w_av = @(s_a,s_v) (s_v^2)/ (s_a^2 + s_v^2);
mu_av = @(mu_a,mu_v,s_a,s_v) w_av(s_a,s_v)*mu_a + (1-w_av(s_a,s_v))*mu_v;

log_likelihood_solo = @(x,mu,st) normlike([mu,st],x);
log_likelihood_sum = @(mu_a,mu_v,st_a,st_v) log_likelihood_solo(training_data(:,1),mu_a,st_a) + log_likelihood_solo(training_data(:,2),mu_v,st_v) + log_likelihood_solo(training_data(:,3),mu_av(mu_a,mu_v,st_a,st_v),sigma_av(st_a,st_v)); 

x_hat_0 = [0,0,1,1];
x_hat_av = fminunc(@(x) log_likelihood_sum(x(1),x(2),x(3),x(4)),x_hat_0);
m_av_hat = mu_av(x_hat_av(1), x_hat_av(2), x_hat_av(3), x_hat_av(4));
sigma_av_hat = sigma_av(x_hat_av(3), x_hat_av(4));

%%
m_a    = x_hat_av(1);
m_v    = x_hat_av(2);
s_a    = x_hat_av(3);
s_v    = x_hat_av(4);
m_av   = m_av_hat;
s_av   = sigma_av_hat;


real_m_a    = mean(x_a);
real_s_a     = std(x_a);
real_m_v    = mean(x_v);
real_s_v    = std(x_v);
real_m_av   = mean(x_av);
real_s_av   = std(x_av);

Variables = {'m_a';'s_a';'m_v';'s_v';'m_av';'s_av'};
Real = [real_m_a;real_s_a;real_m_v;real_s_v;real_m_av;real_s_av];
Calculated = [m_a;s_a;m_v;s_v;m_av;s_av];
T = table (Variables,Real, Calculated)

%%
plot_n =1;
x = -20:1:+30;

figure(plot_n)
fig_idx = 1;
subplot(1,3,fig_idx);
histogram(x_a,'Normalization','probability','displaystyle','stairs','binmethod','integers')
hold on;

y = normpdf(x,m_a,s_a);
plot(x,y,'r','linewidth',2);
title('x_a')

fig_idx = fig_idx + 1;
subplot(1,3,fig_idx);
histogram(x_v,'Normalization','probability','displaystyle','stairs','binmethod','integers')
hold on;

y = normpdf(x,m_v,s_v);
plot(x,y,'r','linewidth',2);
title('x_v')

fig_idx = fig_idx + 1;
subplot(1,3,fig_idx);
histogram(x_av,'Normalization','probability','displaystyle','stairs','binmethod','integers')
hold on;

y = normpdf(x,m_av,s_av);
plot(x,y,'r','linewidth',2);
title('x_a_v')

sgtitle('Normalized Histograms and pdfs')



