clc
close all
clear all

%% Prior probabilities
% juju-fruits
Prior_juju          = 0.10;
Prior_juju_ripe     = 0.15;
Prior_juju_not_ripe = 1-Prior_juju_ripe;
mu_juju_ripe        = 600;
std_juju_ripe       = 50;
mu_juju_not_ripe    = 500;
std_juju_not_ripe   = 50;

% mongo berries
Prior_mongo          = 0.50;
Prior_mongo_ripe     = 0.80;
Prior_mongo_not_ripe = 1-Prior_mongo_ripe;
mu_mongo_ripe        = 580;
std_mongo_ripe       = 20;
mu_mongo_not_ripe    = 520;
std_mongo_not_ripe   = 20;

% chakava
Prior_chakava          = 0.40;
Prior_chakava_ripe     = 0.10;
Prior_chakava_not_ripe = 1-Prior_chakava_ripe;
mu_chakava_ripe        = 400;
std_chakava_ripe       = 100;
mu_chakava_not_ripe    = 550;
std_chakava_not_ripe   = 100;


mu_ripe      = [mu_juju_ripe, mu_mongo_ripe, mu_chakava_ripe];
mu_not_ripe  = [mu_juju_not_ripe, mu_mongo_not_ripe, mu_chakava_not_ripe];
std_ripe     = [std_juju_ripe, std_mongo_ripe, std_chakava_ripe];
std_not_ripe = [std_juju_not_ripe, std_mongo_not_ripe, std_chakava_not_ripe];

prior          = [Prior_juju , Prior_mongo , Prior_chakava ];
prior_ripe     = [Prior_juju_ripe, Prior_mongo_ripe, Prior_chakava_ripe];
prior_not_ripe = [Prior_juju_not_ripe, Prior_mongo_not_ripe, Prior_chakava_not_ripe];
%x=0:800;

like = @likelihood; 

%% Question 1
interval_ripe_1 = [540,550];
likelihood_juju_ripe     = like(interval_ripe_1,mu_juju_ripe,std_juju_ripe);
likelihood_juju_not_ripe = like(interval_ripe_1,mu_juju_not_ripe,std_juju_not_ripe);


P_mystery_juju  = Prior_juju_ripe*likelihood_juju_ripe/( Prior_juju_ripe*likelihood_juju_ripe+ Prior_juju_not_ripe*likelihood_juju_not_ripe);
disp ('The probability of the juju fruit being ripe is :');
disp(P_mystery_juju);



%% Question 2


[response, probability] = is_ripe(interval_ripe_1, prior, prior_ripe, mu_ripe, std_ripe,prior_not_ripe, mu_not_ripe, std_not_ripe);
disp ('The probability of the random fruit being ripe is :');
disp(probability);


%% Question 3

% Labels
% 1 : Juju-Ripe
% 2 : Juju - Not Ripe
% 3 : Mongo - Ripe
% 4 : Mongo - Not Ripe
% 5 : Chakava - Ripe
% 6 : Chakava - Not Ripe
trials = 1000;
p_fruits = [Prior_juju * Prior_juju_ripe Prior_juju * Prior_juju_not_ripe Prior_mongo*Prior_mongo_ripe Prior_mongo*Prior_mongo_not_ripe Prior_chakava*Prior_chakava_ripe Prior_chakava*Prior_chakava_not_ripe];
%p_fruits = [Prior_juju * Prior_juju_ripe, Prior_mongo*Prior_mongo_ripe, Prior_chakava*Prior_chakava_ripe, Prior_juju * Prior_juju_not_ripe, Prior_mongo*Prior_mongo_not_ripe, Prior_chakava*Prior_chakava_not_ripe];

pd = makedist('Multinomial','Probabilities',p_fruits);
fruits = random(pd, trials, 1);

% Normal Distributions:
%total_mu  = [mu_ripe, mu_not_ripe];
%total_std = [std_ripe, std_not_ripe];
matrix_normal = [mu_juju_ripe, std_juju_ripe; mu_juju_not_ripe, std_juju_not_ripe;mu_mongo_ripe, std_mongo_ripe; mu_mongo_not_ripe, std_mongo_not_ripe; mu_chakava_ripe, std_chakava_ripe; mu_chakava_not_ripe, std_chakava_not_ripe ];
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

function x= likelihood(interval,mu,std) 

dummy = normcdf(interval,mu,std);
x = dummy(2)-dummy(1);

end