%% ex 5
clear all;
clc;
close all;
addpath(genpath('FastICA_25'));
%% Import sound data
filname = {'Guitar1.wav', 'Guitar2.wav'}

%%
n_sample = [2,100001];
S = [];
S_dst = [];
S_dummy = [];
Fs = []
j = 1;
for i = 1:2:4
    [S_dummy(:,i:i+1),Fs] = audioread(string(filname(j)), n_sample);
    j= j+1;
end

% linearly combine the two sound track, such that: 
%                   $$x_1(k) = a_{11}\cdot s_1(k) + a_{12}\cdot s_2(k) $$
%                   $$x_2(k) = a_{21}\cdot s_1(k) + a_{22}\cdot s_2(k) $$
A_dist_s = [18 20; -30 +10];
%A_dist_s = [18 20; 18 +20];
S_dst = S_dummy(:,1:2:4);
S(1,:) = A_dist_s(1,1)*S_dst(:,1) + A_dist_s(1,2)*S_dst(:,2) ;
S(2,:) = A_dist_s(2,1)*S_dst(:,1) + A_dist_s(2,2)*S_dst(:,2) ;



%%
%sound(S_sound(:,1))
%sound(S_sound(:,2))
%%

% 1.3 plot the linear combination against each other in a 2D plot
plot_n = 1;
figure(plot_n)

subplot(2,2,1),plot(S_dst(:,1))
xlabel(['k samples']), ylabel(['sound'])
legend('d_1(k)'), title('A'), grid on

subplot(2,2,2), plot(S_dst(:,2))
xlabel(['k samples']), ylabel(['sound'])
legend('d_2(k)'), title('B'), grid on

subplot(2,1,2), plot(S(2,:)), hold on, plot(S(1,:))
xlabel(['k samples']), ylabel(['linear combinations'])
legend('s_1(k)','s_2(k)'), title('C'), grid on
plot_n = plot_n + 1;

figure(plot_n)
plot(S(1,:),S(2,:),'.')
xlabel(['s_1(k) = a_{11}\cdot d_1(k) + a_{12}\cdot d_2(k)'])
ylabel(['s_2(k) = a_{21}\cdot d_1(k) + a_{22}\cdot d_2(k)'])
grid on

plot_n = plot_n + 1;
%%

[W_ica, X_ica, A_ica] = fastica(S);

W_ica_test = A_ica*S;

S_ica       = X_ica(:,:)*W_ica(:,:); % S_ica = X_ica*A_ica*S
S_ica_1     = X_ica(:,1)*W_ica(1,:);
S_ica_2     = X_ica(:,2)*W_ica(2,:);



%%
sound(S_ica_2)
%%
sound(S_dst(:,2))

%% ex 3-4
cov_s     = round(cov(S'),2);
cov_w_ica = round(cov(W_ica'),2);
cov_ica   = round(cov(S_ica'),2);
cov_ica_1 = round(cov(S_ica_1'),2);
cov_ica_2 = round(cov(S_ica_2'),2);

%%
plot_n = 1;
figure(plot_n)
% original data
subplot(2,2,1),plot(S(1,:), S(2,:),'.')
title(splitlines(['A',newline,'Original data',newline,'\rm Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
xlabel(['x_1(k)']), ylabel(['x_2(k)']) 
grid on, axis equal, axis square

% score
subplot(2,2,2), plot(W_ica(1,:),W_ica(2,:),'.')
title(splitlines(['C',newline,'ICA',newline,'\rm score',newline,'Var_1=',num2str(cov_w_ica(1,1)),'; Var_2=',num2str(cov_w_ica(2,2)),'; Cov=',num2str(cov_w_ica(1,2))]  ) )
xlabel(['x_1(k)']), ylabel(['x_2(k)']) 
grid on, axis equal, axis square


% ICA reconstruction
subplot(2,1,2), plot(S_ica(1,:),S_ica(2,:),'.b', S_ica_1(1,:),S_ica_1(2,:),'.m',S_ica_2(1,:),S_ica_2(2,:),'.g')
xlabel(['s^{\prime}_1(k)']), ylabel(['s^{\prime}_2(k)']),legend('data', 'IC_1', 'IC_2')
title(splitlines(['E',newline,'ICA',newline,'\rm Reconstructed data',newline,'Var_1=',num2str(cov_ica(1,1)),'; Var_2=',num2str(cov_ica(2,2)),'; Cov=',num2str(cov_ica(1,2))]  ) )
grid on, axis equal, axis square



%% ex 5.a 
RMSE_reco = sqrt(sum(((S_ica'-S_dst).^2)/size(S,2)))

%% ex.5.b
% Reconstructed sound
sound(S_ica)
%% Original sound
sound(S_dst)
