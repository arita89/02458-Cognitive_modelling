%% P2 - SDT - N and SN with unequal variance
clear
clc
format short;  
dash = ('----------------------------------');
%%
tot_obs = 50; 
n = 3; % NUMBER CLASSES

N = [15,6,29]; %Noise [yes, no, maybe]
N_yes = 29; %FA
N_no = 15;
N_maybe = 6;
%------------------------------
S = [27,4,19]; %Signal(+noise)
S_yes = 19; %H
S_no = 27;
S_maybe = 4;

%% cumulatives

Cum_N = cumsum(N)
Cum_S = cumsum(S)

%% probabilities


P_FA = [(6+29)/50, (29)/50]
P_H = [(4+19)/50,(19)/50]


%% standardize
Z_FA = norminv(P_FA)
Z_H = norminv(P_H)

d_prime= Z_H-Z_FA;

%% plot in Gaussian coordinates the ROC
% receiver operating caratheristic ROC curve  z_P(hit) over z_P(False alarm)
% z_H = (1/sigma_signal)* z_FA + (mu_signal/sigma_signal)

figure(1)
hold on
x = Z_FA
y = Z_H
scatter(x,y,25,'b','o') 

figure(2)
hold on
x = Z_FA;
y = Z_H;
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);;%<-- linear fit
slope = P(1);
% if slope of ROC curve is <1 the standard deviation of Signal-noise distribution is
% greater then the standad deviation of the Noise
hold on;
plot(x,yfit,'b-'); 
xlabel('Z_F_A')
ylabel('Z_H')
legend( 'data fit','model')
title ('Receiver Operating Characteristics-gaussian Coordinates')
grid on

figure(3) %ROC in P_H and P_FA
hold on
x = P_FA;
y = P_H;
scatter(x,y,25,'b','o') 
hold on;
plot(x,y,'b-');  
xlabel('p_F_A')
ylabel('P_H')
legend('data', 'isosensitivity contour')
title ('Receiver Operating Characteristics')
grid on

%%retrieve mu and sigma
% z_H = (1/sigma_signal)* z_FA + (mu_signal/sigma_signal)
% yfit = P(1)*x+P(2);
% P(1)=(1/sigma_signal)
% P(2)=(mu_signal/sigma_signal)
%mu = mean(N);
%sigma = std(N);

sigma_SN_approx= 1/P(1)
mu_SN_approx = P(2)*sigma_SN_approx
