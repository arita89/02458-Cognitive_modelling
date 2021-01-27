clc
close all
clear all

%% Bayesian Likelihood

%% Data
% J
mu_J_R = 600;
sigma_J_R = 50; 

mu_J_NR = 500;
sigma_J_NR = 50; 

%M
mu_M_R = 580;
sigma_M_R = 20; 

mu_M_NR = 520;
sigma_M_NR = 20; 

%C
mu_C_R = 400;
sigma_C_R = 100; 

mu_C_NR = 550;
sigma_C_NR = 100; 

%--------------------------------------------------------------------------
% Priors
P_J = 0.10;
P_J_R = 0.15;
P_J_NR = 1-P_J_R;

P_M = 0.50;
P_M_R = 0.80;
P_M_NR = 1-P_M_R;

P_C = 1-P_J-P_M;
P_C_R = 0.10;
P_C_NR = 1-P_C_R;

%% Ex 1P(J_R| 540<X<550)
% posterior probability oF J_R given Interval

%P(J_R| 540<X<550) =  [P(540<X<550 | J_R) * (P_J_R) ] / P(540<X<550 | J)

%----------------------------------------------------------------------------
% numerator [P(540<X<550 | J_R) * (P_J_R) ]
% P(540<X<550 | J_R) = interval_likelihood(interval,mu,std)
% (P_J_R)
interval = [540,550];
numerator = interval_likelihood(interval,mu_J_R,sigma_J_R)*P_J_R;
%----------------------------------------------------------------------------
% denominator P(540<X<550) 
% = P(540<X<550 | J_R)*(P_J_R) +P(540<X<550 | J_NR)*(P_J_NR)
denominator = interval_likelihood(interval,mu_J_R,sigma_J_R)*P_J_R+ interval_likelihood(interval,mu_J_NR,sigma_J_NR)*P_J_NR;
%----------------------------------------------------------------------------
%%solution ex1
P_ex1= numerator/denominator

%% Ex 2 P(random_fruit_R| 540<X<550)

% =[(J_R| 540<X<550)*P_J+(M_R| 540<X<550)*P_M+(C_R| 540<X<550)*P_C] / P(540<X<550) 

% numerator = same as previous for each fruit
%? = ?1 ? P(Juju) + ?2 ? P(mongos) + ?3? P(chakavas)
%?1 = P(540-550| Juju,ripe) ? P(Juju,ripe),
%?2 = P(540-550| mongos,ripe) ? P(mongos,ripe),
%?3 = P(540-550| chakavas,ripe) ? P(chakavas,ripe),

% interval, mu_X_R, sigma_X_R, mu_X_NR, sigma_X_NR,P_X_R,P_X_NR 
pp_J_R = Num(interval, mu_J_R, sigma_J_R, mu_J_NR, sigma_J_NR,P_J_R,P_J_NR)
pp_M_R = Num(interval, mu_M_R, sigma_M_R, mu_M_NR, sigma_M_NR,P_M_R,P_M_NR)
pp_C_R = Num(interval, mu_C_R, sigma_C_R, mu_C_NR, sigma_C_NR,P_C_R,P_C_NR)
numerator2= pp_J_R*P_J+pp_M_R*P_M+pp_C_R*P_C;

% denominator P(540<X<550) 
% = sum(over X) of [P(540<X<550 | X_R)*(P_X_R) +P(540<X<550 |
% X_NR)*(P_X_NR)]*P_X

pp_I_J = Den(interval, mu_J_R, sigma_J_R, mu_J_NR, sigma_J_NR,P_J_R,P_J_NR)
pp_I_M = Den(interval, mu_M_R, sigma_M_R, mu_M_NR, sigma_M_NR,P_M_R,P_M_NR)
pp_I_C = Den(interval, mu_C_R, sigma_C_R, mu_C_NR, sigma_C_NR,P_C_R,P_C_NR)

denominator2 =pp_I_J*P_J+pp_I_M*P_M+pp_I_C*P_C;

%----------------------------------------------------------------------------
%%solution ex2
P_ex2= numerator2/denominator2

%% Ex 3 


% Labels
% 1 : Juju-Ripe
% 2 : Juju - Not Ripe
% 3 : Mongo - Ripe
% 4 : Mongo - Not Ripe
% 5 : Chakava - Ripe
% 6 : Chakava - Not Ripe
trials = 1000;
p_fruits = [P_J * P_J_R P_J * P_J_NR P_M*P_M_R P_M*P_M_NR P_C*P_C_R P_C*P_C_NR];
%p_fruits = [Prior_juju * Prior_juju_ripe, Prior_mongo*Prior_mongo_ripe, Prior_chakava*Prior_chakava_ripe, Prior_juju * Prior_juju_not_ripe, Prior_mongo*Prior_mongo_not_ripe, Prior_chakava*Prior_chakava_not_ripe];

pd = makedist('Multinomial','Probabilities',p_fruits);
fruits = random(pd, trials, 1);

% Normal Distributions:
%total_mu  = [mu_ripe, mu_not_ripe];
%total_std = [std_ripe, std_not_ripe];
matrix_normal = [mu_J_R, sigma_J_R; mu_J_NR, sigma_J_NR; mu_M_R, sigma_M_R; mu_M_NR, sigma_M_NR; mu_C_R, sigma_C_R; mu_C_NR, sigma_C_NR ];
%matrix_normal = [mu_juju_ripe, std_juju_ripe; mu_mongo_ripe, std_mongo_ripe; mu_chakava_ripe, std_chakava_ripe; mu_juju_not_ripe, std_juju_not_ripe; mu_mongo_not_ripe, std_mongo_not_ripe; mu_chakava_not_ripe, std_chakava_not_ripe ];


counter = 0;
for i = 1:trials
    
    mu    = matrix_normal(fruits(i),1);
    sigma = matrix_normal(fruits(i),2);
    
    r = normrnd (mu , sigma);
    interval_ripe_2 = [r-5,r+5];
    [response, probability] = is_ripe(interval_ripe_2, prior, prior_ripe, mu_ripe, std_ripe,prior_not_ripe, mu_not_ripe, std_not_ripe);
    
    if response == rem(fruits(i),2)
        counter = counter + 1;
    end
end

disp ('The accuracy of the monkey is :');
disp(counter/1000);


%% Functions
function x= likelihood(b,mu,std) %likelihood up to value X<b
    x = normcdf(b,mu,std);
end

function x= interval_likelihood(interval,mu,std) %likelihood of interval  a<X<b
    dummy = normcdf(interval,mu,std);
    x = dummy(2)-dummy(1); %L(X<b)-L(X<a)
end

function x = posterior_probability(interval, mu_X_R, sigma_X_R, mu_X_NR, sigma_X_NR,P_X_R,P_X_NR)
    N = interval_likelihood(interval,mu_X_R,sigma_X_R)*P_X_R;
    D = interval_likelihood(interval,mu_X_R,sigma_X_R)*P_X_R+ interval_likelihood(interval,mu_X_NR,sigma_X_NR)*P_X_NR;
    x = N/D;
end

function x = Num(interval, mu_X_R, sigma_X_R, mu_X_NR, sigma_X_NR,P_X_R,P_X_NR)
    x = interval_likelihood(interval,mu_X_R,sigma_X_R)*P_X_R;
end

%P(540<X<550 | X_R)*(P_X_R) +P(540<X<550 | X_NR)*(P_X_NR)
function x = Den(interval, mu_X_R, sigma_X_R, mu_X_NR, sigma_X_NR,P_X_R,P_X_NR)
    x = interval_likelihood(interval,mu_X_R,sigma_X_R)*P_X_R+ interval_likelihood(interval,mu_X_NR,sigma_X_NR)*P_X_NR;
end

function [response, probability] = is_ripe(interval, prior_exist, prior_ripe, mu_ripe, std_ripe,prior_not_ripe, mu_not_ripe, std_not_ripe) 


t = 0;
b = 0;
    for i= 1:length(mu_ripe)
        ripe_like = normcdf(interval,mu_ripe(i),std_ripe(i));
        ripe_like_int = ripe_like(2)-ripe_like(1);
        not_ripe_like = normcdf(interval,mu_not_ripe(i),std_not_ripe(i));
        not_ripe_like_int = not_ripe_like(2)-not_ripe_like(1);
        
        dummy = prior_ripe(i) * ripe_like_int;
        t = t + dummy * prior_exist(i); 
        b = b + prior_exist(i) * ( dummy + prior_not_ripe(i)*not_ripe_like_int);

    end
probability = t/b;
response = 0;
    if probability > 0.5
        response = 1;
    end
end