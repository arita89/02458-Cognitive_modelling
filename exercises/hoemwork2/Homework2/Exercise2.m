close all
clear all

%% probabilities of fruits being ripe
P_juju_ripe = @(x) (0.15 * normpdf(x, 600, 50));
P_juju_not_ripe =@(x) (0.85 * normpdf(x, 500, 50));
P_mongo_ripe = @(x) (0.80 * normpdf(x, 580, 20));
P_mongo_not_ripe =@(x) (0.20 * normpdf(x, 520, 20));
P_chakava_ripe = @(x) (0.10 * normpdf(x, 400, 100));
P_chakava_not_ripe =@(x) (0.90 * normpdf(x, 550, 100));

x=0:800;
%% question 1
P_mystery_juju  = (integral(P_juju_ripe,540,550)*0.15)/(integral(P_juju_ripe,540,550)+integral(P_juju_not_ripe,540,550));
disp ('The probability of the juju fruit being ripe is :');
disp(P_mystery_juju);


plot(P_juju_ripe(x)); hold on;
plot(P_juju_not_ripe(x));
legend('juju ripe', 'juju not ripe');

%% question 2
t1 = (integral(P_juju_ripe,540,550)*0.15*0.10); 
t2 = (integral(P_mongo_ripe,540,550)*0.80*0.50);
t3 = (integral(P_chakava_ripe,540,550)*0.10*0.40);

b1 = 0.10*(integral(P_juju_ripe,540,550)+integral(P_juju_not_ripe,540,550));
b2 = 0.50*(integral(P_mongo_ripe,540,550)+integral(P_mongo_not_ripe,540,550));
b3 = 0.40*(integral(P_chakava_ripe,540,550)+integral(P_chakava_not_ripe,540,550));

P_random_fruit = (t1+t2+t3)/(b1+b2+b3);

disp ('The probability of the random fruit being ripe is :');
disp(P_random_fruit);

figure();
plot(0.10 * P_juju_ripe(x)); hold on;
plot(0.10 * P_juju_not_ripe(x)); hold on;
plot(0.50 * P_mongo_ripe(x)); hold on;
plot(0.50 * P_mongo_not_ripe(x)); hold on;
plot(0.40 * P_chakava_ripe(x)); hold on; 
plot(0.40 * P_chakava_not_ripe(x));
legend('juju ripe', 'juju not ripe', 'mongo ripe', 'mongo not ripe', 'chakava ripe', 'chakava not ripe');
