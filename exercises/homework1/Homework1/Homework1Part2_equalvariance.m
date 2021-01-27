%%Homework1-part2
clear
clc
format short;  

mu =0;
sigma =1;
t = 1000;
t_n= t/2; % how much of the input is noise N and how much is signal plus noise sn
t_sn= t-t_n;
C = [-0.5 0.5 1.5];%<-- CRITERION (ALSO C or Beta)to modify to change the observer attitude! 
d= ones(size(C));
N= normrnd(0,1,[t_n 1]);%noise
%mu = mean(N);
%sigma = std(N);
SN = normrnd (1,1,[t_sn 1]);%signal plus noise
disp('The Criterions are')
disp (C)
%%
FA= zeros(size(C));% False Alarm
CR= zeros(size(C));% Correct Rejection
H=zeros(size(C));% Hit
M=zeros(size(C));% Miss
P_FA= zeros(size(C));%  p(yes/N)
P_CR= zeros(size(C));%  p(no/N)
P_H=zeros(size(C));% p(Hit)
P_M=zeros(size(C));% p(Miss)

%p= zeros(1,4);
for j =1:3
%%Noise
 for i=1:(t_n)
    if N(i) <= C(j)
        CR(j)=CR(j)+1;
    else FA(j)=FA(j)+1;
    end
 end

%%Signal and Noise
 for i=1:(t_sn)
    if SN(i) > C(j) 
        H(j)=H(j)+1;%--> p(yes/SN)
    else M(j)=M(j)+1;% --> p(no/SN)
    end
 end
 
P_FA(j)= FA(j)/t_n;%  p(yes/N)
P_CR(j)= CR(j)/t_n;%  p(no/N)
P_H(j)= H(j)/t_sn;% p(Hit)
P_M(j)= M(j)/t_sn;% p(Miss)
end

p = [P_H P_FA P_CR P_M];
%P_H_z = zscore(P_H);
%P_FA_z = zscore(P_FA); 
%P_CR_z =zscore(P_CR); 
%P_M_z = zscore(P_M); % standardization of probabilities 

T1=table(C',H',M',FA',CR');
T1.Properties.VariableNames = {'Criterions','H','M','FA','CR'};
T1

T2=table(C',P_H',P_M',P_FA',P_CR');
T2.Properties.VariableNames = {'Criterions','P_H','P_M','P_FA','P_CR'};
T2


%calculate d_prime for each experiment
%% ANSWER  for a particular separation of Noise and Signal+Noise distributions the value of d' is CONSTANT for all possible criterion positions
for j =1:3
 z_FA(j) = norminv(P_FA(j));
 z_N(j) = norminv(1.0-P_FA(j));% z_score(1.0-p(FA))-> z_N
 z_H(j)= norminv(P_H(j));
 z_SN(j)= norminv(1.0-P_H(j));% z_score(1.0-p(H))-> z_SN
 d_approx(j) = z_N(j)-z_SN(j);% dprime = z_N-z_SN = z(pHit)-z(pFA)= z_H-z_FA 
end

OUT = ['Number of trials:',num2str(t)];
disp(OUT)
T3=table(C',d',d_approx');
T3.Properties.VariableNames = {'Criterions','d_prime','d_approx' };
T3

%%plot in Gaussian coordinates the ROC
% receiver operating caratheristic ROC curve  z_P(hit) over z_P(False alarm)
% z_H = (1/sigma_signal)* z_FA + (mu_signal/sigma_signal)

figure(1)
x = z_FA;
y = z_H;
scatter(x,y,25,'b','o') 
title ('Data points') 

figure(3)
hold on
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
%yfit =x+P(1);%<-- fit
hold on;
plot(x,yfit,'b-'); 
for j =1:3
m(j)= x(j)+ d_approx(j); % z_H = z_FA+ d_approx <--- mymodel
end
plot(x,m,'r-') 
xlabel('Z_F_A')
ylabel('Z_H')
legend('data fit','model')
title ('Receiver Operating Characteristics-gaussian Coordinates')
grid on

figure(2) %ROC in P_H and P_FA
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

%%new estimate of d' with the fit
%d' = norminv(P_Hit)+C
for j =1:3
d_eval1(j) = z_H(j)+C(j);
end
%or
%yfit = P(1)*x+P(2);
%P(2) = yfit-P(1)^*x;
%for x=0.1
%d' = yfit-P(1)*0.1
%d_eval2 = yfit-P(1)*0.1

OUT = ['Number of trials:',num2str(t)];
disp(OUT)
T4=table(C',d',d_approx',d_eval1');
T4.Properties.VariableNames = {'Criterions','d_prime','d_approx', 'd_eval_ROC' };
T4

