close all
clear all

%% Prior probabilities
P_juju_ripe = @(x) (0.15 * normpdf(x, 600, 50));
P_juju_not_ripe =@(x) (0.85 * normpdf(x, 500, 50));
P_mongo_ripe = @(x) (0.80 * normpdf(x, 580, 20));
P_mongo_not_ripe =@(x) (0.20 * normpdf(x, 520, 20));
P_chakava_ripe = @(x) (0.10 * normpdf(x, 400, 100));
P_chakava_not_ripe =@(x) (0.90 * normpdf(x, 550, 100));

%x=0:800;
%% Question 1
P_mystery_juju  = (integral(P_juju_ripe,540,550)) / (integral(P_juju_ripe,540,550) + integral(P_juju_not_ripe,540,550));
disp ('The probability of the juju fruit being ripe is :');
disp(P_mystery_juju);

% plot(P_juju_ripe(x)); hold on;
% plot(P_juju_not_ripe(x));
% legend('juju ripe', 'juju not ripe');

%% Question 2
[response, probability] = is_ripe(540, 550);
disp ('The probability of the random fruit being ripe is :');
disp(probability);

% figure();
% plot(0.10 * P_juju_ripe(x)); hold on;
% plot(0.10 * P_juju_not_ripe(x)); hold on;
% plot(0.50 * P_mongo_ripe(x)); hold on;
% plot(0.50 * P_mongo_not_ripe(x)); hold on;
% plot(0.40 * P_chakava_ripe(x)); hold on; 
% plot(0.40 * P_chakava_not_ripe(x));
% legend('juju ripe', 'juju not ripe', 'mongo ripe', 'mongo not ripe', 'chakava ripe', 'chakava not ripe');

%% Question 3

% Labels
% 1 : Juju-Ripe
% 2 : Juju - Not Ripe
% 3 : Mongo - Ripe
% 4 : Mongo - Not Ripe
% 5 : Chakava - Ripe
% 6 : Chakava - Not Ripe

p_fruits = [0.10*0.15 0.10*0.85 0.50*0.80 0.50*0.20 0.40*0.10 0.40*0.9];
pd = makedist('Multinomial','Probabilities',p_fruits);
fruits = random(pd, 1000, 1);

% Normal Distributions:
matrix_normal = [600, 50; 500, 50; 580, 20; 520, 20; 400, 100; 550, 100];
counter = 0;
for i = 1:1000
    mu = matrix_normal (fruits(i), 1);
    sigma = matrix_normal (fruits(i), 2);
    
    r = normrnd (mu , sigma);
    [response, probability] = is_ripe(r - 5, r + 5);
    
    if response == rem(fruits(i),2)
        counter = counter + 1;
    end
end

disp ('The accuracy of the monkey is :');
disp(counter/1000);

%% Functions
function [response, probability] = is_ripe(a, b) 

    P_juju_ripe = @(x) (0.15 * normpdf(x, 600, 50));
    P_juju_not_ripe =@(x) (0.85 * normpdf(x, 500, 50));
    P_mongo_ripe = @(x) (0.80 * normpdf(x, 580, 20));
    P_mongo_not_ripe =@(x) (0.20 * normpdf(x, 520, 20));
    P_chakava_ripe = @(x) (0.10 * normpdf(x, 400, 100));
    P_chakava_not_ripe =@(x) (0.90 * normpdf(x, 550, 100));

    t1 = (0.10 * integral(P_juju_ripe, a, b)); 
    t2 = (0.50 * integral(P_mongo_ripe, a, b));
    t3 = (0.40 * integral(P_chakava_ripe, a, b));

    b1 = 0.10 * (integral(P_juju_ripe, a, b) + integral(P_juju_not_ripe, a, b));
    b2 = 0.50 * (integral(P_mongo_ripe, a, b) + integral(P_mongo_not_ripe, a, b));
    b3 = 0.40 * (integral(P_chakava_ripe, a, b) + integral(P_chakava_not_ripe, a, b));

    probability = (t1+t2+t3)/(b1+b2+b3);
    response = 0;
    if probability > 0.5
        response = 1;
    end
end