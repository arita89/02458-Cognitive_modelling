clc;
clear all;
close all;
dash = ('----------------------------------');
%% SIGNAL DETECTION 
% 1 partecipant
% confidence rate of her responses: high,medium, low
rate = ['H','M','L'];
tot_obs = 100; %number of tests for Noise and Sound
l =3; %number of test subjects/type of answers in this case
%----------------------
% NO sound played (Noise , N)

CR = [2,31,57]
FA = [0,2,8]
%----------------------
% YES sound played (SignalplusNoise , SN)
H = [38,15,4]
M = [27,14,2] %<------careful! invert!


%% calculate probabilities

P_H = H/tot_obs
P_FA = FA/tot_obs
P_M = M/tot_obs
P_CR = CR/tot_obs
%disp (dash)

%% d_prime - Definition
dprime = ones(1,l);
d_prime = @(Prob_hits,Prob_FalseAlarm) norminv(Prob_hits)-norminv(Prob_FalseAlarm);

for i = 1:l
    %disp (i)
    dprime(i) = d_prime (P_H(i),P_FA(i));
    %disp (dash)
end

dprime

%% criterion 
% criterion = lambda = - (norminv(P_FA))
lambda = -(norminv(P_FA));

%%
disp(dash)
for  i = 1:l
    disp ([rate(i),' confidence'])
    disp(['The sensitivity is ',num2str(dprime(i))])
    disp(['The bias is ',num2str(lambda(i))])
    disp(dash)
end 
