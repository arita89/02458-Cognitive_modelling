%%Homework1-part3
clear
clc
format short;  

%mu =0;
%sigma =1;
%ex = ['Ex1) moderate criterion' 'Ex2) conservative criterion' 'Ex3) lax criterion'];
IL = 1:5; %five intensiy levels
x = [0.4 0.9 1.2 1.7 2.3];%<-- CRITERION (ALSO C or Beta)to modify to change the observer attitude! 
N = length(x);

n =50
H = [1 6 13  32 49]; %P(hit)  = 'yes' rate
P_H = H/n;

IL = 1:5;
disp(' Int Lev   Intensity     P(Hit)   P(Miss)   ')
disp ([IL' x' P_H' (1-P_H)'])


figure(1)
hold on
scatter(x,P_H,25,'b','o')
xlabel('stimulus intensity')
ylabel('probablity of detection ')
legend('sperimental value of hits')
grid on

fcn = @(b,x) normcdf(x, b(1), b(2)); 
MSE = @(b) norm(P_H - fcn(b,x));   % Norm Residual
B = fminsearch(MSE, [0; 10]);   % Estimate Parameters

% hp, the parameters are the ones that are used in stevens law?
% phi = B(2)*x.^(B(1));

Xplot = linspace(min(x), max(x));
figure(2)
plot(x, P_H, 'bo')
hold on
plot(Xplot, fcn(B,Xplot))
plot(x, phi, 'r')
xlabel('stimulus intensity')
ylabel('probablity of detection ')
%legend('sperimental value of hits', 'data fit','stevens model of perception of stimulus')
legend('sperimental value of hits', 'data fit')
hold off
grid

fcn(1,B)
fcn(2,B)