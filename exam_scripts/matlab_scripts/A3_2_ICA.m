clc;
clear all;
close all;

addpath(genpath('FastICA_25'));

%% MODELING COGNITION - ICA EXERCISE
%
%  exercise 1
%  PCA analysis on a bivariate distribution with excess kurtosis

% 1.1 create two lapalcian distributions s_1(k) and s_2(k)
mu_s = 0;
sigma_s = 0.5;
dim_s = [2,10000];
dist_s = randlpl(mu_s,sigma_s,dim_s(1),dim_s(2));
A_dist_s = [10 20; -2 +10]; %matrix of weights

% 1.2 linearly combine the two distributions, such that: 
%                   $$x_1(k) = a_{11}\cdot s_1(k) + a_{12}\cdot s_2(k) $$
%                   $$x_2(k) = a_{21}\cdot s_1(k) + a_{22}\cdot s_2(k) $$
s_1 = A_dist_s(1,1)*dist_s(1,:) + A_dist_s(1,2)*dist_s(2,:); 
s_2 = A_dist_s(2,1)*dist_s(1,:) + A_dist_s(2,2)*dist_s(2,:);  

%%
% 1.3 plot the linear combination against each other in a 2D plot
plot_n = 1;
figure(plot_n)

subplot(3,1,1), plot(dist_s(1,:)), hold on, plot(dist_s(2,:))
xlabel(['k samples']), ylabel(['bivariate distributions'])
legend('d_1(k)','d_2(k)'), title('A'), grid on

subplot(3,1,2), plot(s_1), hold on, plot(s_2)
xlabel(['k samples']), ylabel(['linear combinations'])
legend('s_1(k)','s_2(k)'), title('B'), grid on

subplot(3,1,3); %same thing as below
scatter(s_1,s_2,'filled')
xlabel(['s_1(k)']), ylabel(['s_2(k)'])
title('s_1 against s_2'), grid on
plot_n = plot_n + 1;

figure(plot_n)
plot(s_1,s_2,'.')
xlabel(['s_1(k) = a_{11}\cdot d_1(k) + a_{12}\cdot d_2(k)'])
ylabel(['s_2(k) = a_{21}\cdot d_1(k) + a_{22}\cdot d_2(k)'])
grid on

plot_n = plot_n + 1;

%% ex 2 - PCA on the joint distribution

S = [s_1; s_2]; %column 
display('Kurtosis excess')
kurtosis(S);

% Principal Components analysis
[X,W,latent] = pca(S'); %take transpose to have raw vector

% check how many eigenvalues are needed to have 95% variance
i = 0;
sum_variance = 0;
var_threshold = sum(latent)*0.95;
while sum_variance <= var_threshold
   i = i+1;
   sum_variance = sum_variance + latent(i);
end
display('Min PCA necessary')
n_pca = i

figure(plot_n), bar(latent)
title( 'Principal components variance' )
xlabel('PC'), ylabel('Variance'),grid on
plot_n = plot_n+1;

%% reconstruction of S
S_pca = X*W';
S_pca_1 = X(:,1)*W(:,1)';
S_pca_2 = X(:,2)*W(:,2)';


cov_s     = round(cov(S'),2);
cov_W     = round(cov(W),2);
cov_s_pca = round(cov(S_pca'),2);

figure(plot_n)

% Original data
subplot(2,2,1),plot(S(1,:), S(2,:),'.')
title(splitlines(['A - Original data',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
xlabel(['s_1(k)']), ylabel(['s_2(k)']) 
grid on,axis equal,axis square

% PCA score
subplot(2,2,2),plot(W(:,1), W(:,2),'.')
title(splitlines(['B - Scores',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
xlabel(['s_1(k)']), ylabel(['s_2(k)']) 
grid on,axis equal,axis square

% Reconstructed data with PCA
subplot(2,1,2), plot(S_pca(1,:), S_pca(2,:),'.b'), hold on, plot(S_pca_1(1,:),S_pca_1(2,:),'.r'), hold on, plot(S_pca_2(1,:),S_pca_2(2,:),'.c')
title(splitlines(['C - Reconstructed data',newline,'Var_1=',num2str(cov_s_pca(1,1)),'; Var_2=',num2str(cov_s_pca(2,2)),'; Cov=',num2str(cov_s_pca(1,2))]  ) )
xlabel(['s^{\prime}_1(k)']),ylabel(['s^{\prime}_2(k)']),legend('data', 'PC_1', 'PC_2','Location','southeast')
grid on, axis equal,axis square
sgtitle('Principal components analysis')
plot_n = plot_n+1;


figure(plot_n)
figure_title = ['A','B','C','D'];
legend_title = ['var';'cov';'cov';'var'];
subp_id = 1;
max_pca = 2;
for pca_idx_row = 1:1:max_pca
    for pca_idx_col = 1:1:max_pca
        subplot(max_pca,max_pca,subp_id),plot(S_pca(pca_idx_row,:), S_pca(pca_idx_col,:),'.')
        title(figure_title(subp_id)), xlabel('PC_'+string(pca_idx_row)),ylabel('PC_'+string(pca_idx_col))
        legend([legend_title(subp_id,:),' ',num2str(cov_s_pca(pca_idx_row,pca_idx_col))])
        grid on, axis equal
        subp_id = subp_id + 1;
    end
end
sgtitle('Principal components analysis')

plot_n = plot_n+1;

%% ex 3 Independent components analysis - ICA

[W_ica, X_ica, A_ica] = fastica(S);

W_ica_test = A_ica*S;

S_ica       = X_ica(:,:)*W_ica(:,:); % S_ica = X_ica*A_ica*S
S_ica_1     = X_ica(:,1)*W_ica(1,:);
S_ica_2     = X_ica(:,2)*W_ica(2,:);

cov_ica     = round(cov(W_ica'),2);
cov_s_ica   = round(cov(S_ica'),2);

%% ex 3-4 plot score, and then compare it with the original data

figure(plot_n)
% original data
subplot(2,3,1),plot(S(1,:), S(2,:),'.')
title(splitlines(['A',newline,'Original data',newline,'\rm Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
xlabel(['x_1(k)']), ylabel(['x_2(k)']) 
grid on, axis equal, axis square

% score
subplot(2,3,2),plot(W(:,1), W(:,2),'.')
title(splitlines(['B',newline,'PCA',newline,'\rm score',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
xlabel(['x_1(k)']), ylabel(['x_2(k)']), 
grid on, axis equal, axis square


subplot(2,3,3), plot(W_ica(1,:),W_ica(2,:),'.')
title(splitlines(['C',newline,'ICA',newline,'\rm score',newline,'Var_1=',num2str(cov_ica(1,1)),'; Var_2=',num2str(cov_ica(2,2)),'; Cov=',num2str(cov_ica(1,2))]  ) )
xlabel(['x_1(k)']), ylabel(['x_2(k)']) 
grid on, axis equal, axis square

% PCA reconstruction
subplot(2,2,3), plot(S_pca(1,:), S_pca(2,:),'.b', S_pca_1(1,:),S_pca_1(2,:),'.r', S_pca_2(1,:),S_pca_2(2,:),'.c')
title(splitlines(['D',newline,'PCA',newline,'\rm Reconstructed data',newline,'Var_1=',num2str(cov_s_pca(1,1)),'; Var_2=',num2str(cov_s_pca(2,2)),'; Cov=',num2str(cov_s_pca(1,2))]  ) )
xlabel(['s^{\prime}_1(k)']), ylabel(['s^{\prime}_2(k)']), legend('data', 'PC_1', 'PC_2')
grid on, axis equal, axis square

% ICA reconstruction
subplot(2,2,4), plot(S_ica(1,:),S_ica(2,:),'.b', S_ica_1(1,:),S_ica_1(2,:),'.m',S_ica_2(1,:),S_ica_2(2,:),'.g')
xlabel(['s^{\prime}_1(k)']), ylabel(['s^{\prime}_2(k)']),legend('data', 'IC_1', 'IC_2')
title(splitlines(['E',newline,'ICA',newline,'\rm Reconstructed data',newline,'Var_1=',num2str(cov_s_ica(1,1)),'; Var_2=',num2str(cov_s_ica(2,2)),'; Cov=',num2str(cov_s_ica(1,2))]  ) )
grid on, axis equal, axis square

%%

plot_n = plot_n+1;
figure(plot_n)
subp_id = 1;
max_ica = 2;
for pca_idx_row = 1:1:max_ica
    for pca_idx_col = 1:1:max_ica
        subplot(max_ica,max_ica,subp_id),plot(S_ica(pca_idx_row,:), S_ica(pca_idx_col,:),'.')
        xlabel('IC_'+string(pca_idx_row)),ylabel('IC_'+string(pca_idx_col))
        legend([legend_title(subp_id,:),' ',num2str(cov_s_ica(pca_idx_row,pca_idx_col))])
        grid on, axis equal
        subp_id = subp_id + 1;
    end
end
sgtitle('Independent components analysis')
plot_n = plot_n+1;

%% functions

function x = randlpl(mu, b, m, n)
    x = mu - b*sign(rand(m,n)-0.5).*log(1-2*abs(rand(m,n)-0.5));
end  