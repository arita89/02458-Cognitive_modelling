clc;
clear all;
close all;

%% ex 1
mu_1 = 0;
sigma_1 = 0.5;
dim = [2,10000];
dist_1 = randlpl(mu_1,sigma_1,dim(1),dim(2));
x_comb_1 = dist_1(1,:)+ 2*dist_1(2,:); %c1 = x1+2*x2
x_comb_2 = dist_1(2,:) - dist_1(1,:); %c2 =  x1 -x2
plot(x_comb_1,x_comb_2,'.')
xlabel(['x_1 = d_1 + 2 \cdot d_2'])
ylabel(['x_2 = d_1 - d_2'])
grid on

%% ex 2 - PCA on the joint distribution

S = [x_comb_1' x_comb_2'];
% Principal Components
[X,W,latent] = pca(S);

% check how many eigenvalues are needed to have 95% variance
i = 0;
sum_variance = 0;
var_threshold = sum(latent)*0.95;
while sum_variance <= var_threshold
   i = i+1;
   sum_variance = sum_variance + latent(i);
end
n_pca = i;

% reconstruction of S, except for the scaling, in the pc-space
S_pca = X*W';


plot_n = 1;
figure(plot_n), bar(latent)
title( 'pca variance' )
grid on
plot_n = plot_n+1;

cov_s = round(cov(S ),2);
cov_W = round(cov(W ),2);

figure(plot_n)
subplot(1,2,1),plot(S(:,1), S(:,2),'.')
title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
grid on
subplot(1,2,2),plot(W(:,1), W(:,2),'.')
xlabel('PCA 1'),ylabel('PCA 2')
title(splitlines(['The score',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
grid on
axis equal

plot_n = plot_n+1;


figure(plot_n)
subp_id = 1;
max_pca = 2;
for pca_idx_row = 1:1:max_pca
    for pca_idx_col = 1:1:max_pca
        subplot(max_pca,max_pca,subp_id),plot(W(:,pca_idx_row), W(:,pca_idx_col),'.')
        xlabel('PCA '+string(pca_idx_row)),ylabel('PCA '+string(pca_idx_col))
        legend(num2str(cov_W(pca_idx_row,pca_idx_col)))
        grid on
        axis equal
        subp_id = subp_id + 1;
    end
end

plot_n = plot_n+1;

%% ex 3
[ica_vect_sound, A_ica, W_ica] = fastica(S');
S_ica_w = W_ica*ica_vect_sound; % remember to argument about this, not fully reconstructed 

S_ica = A_ica*ica_vect_sound; % = A_ica*S*W_ica

%% ex 3-4
cov_ica = round(cov(ica_vect_sound'),2);
cov_s_ica = round(cov(S_ica'),2);
cov_s_pca = round(cov(S_pca'),2);
figure(plot_n)
subplot(2,4,1),plot(S(:,1), S(:,2),'.')
title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
grid on
axis equal
subplot(2,4,2),plot(W(:,1), W(:,2),'.')
xlabel('PCA 1'),ylabel('PCA 2')
title(splitlines(['The score',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
grid on
axis equal

subplot(2,4,3), plot(ica_vect_sound(1,:),ica_vect_sound(2,:),'.')
xlabel('ICA 1'),ylabel('ICA 2')
title(splitlines(['The score',newline,'Var_1=',num2str(cov_ica(1,1)),'; Var_2=',num2str(cov_ica(2,2)),'; Cov=',num2str(cov_ica(1,2))]  ) )
grid on
axis equal

subplot(2,2,3), plot(S_pca(1,:),S_pca(2,:),'.')
xlabel('reconstructed 1'),ylabel('reconstructed 2')
title(splitlines(['Reconstructed pca data',newline,'Var_1=',num2str(cov_s_pca(1,1)),'; Var_2=',num2str(cov_s_pca(2,2)),'; Cov12=',num2str(cov_s_pca(1,2))]  ) )
grid on
axis equal
subplot(2,2,4), plot(S_ica(1,:),S_ica(2,:),'.b', S_ica_w(1,:),S_ica_w(2,:),'.r')
xlabel('reconstructed 1'),ylabel('reconstructed 2')
title(splitlines(['Reconstructed ica data',newline,'Var_1=',num2str(cov_s_ica(1,1)),'; Var_2=',num2str(cov_s_ica(2,2)),'; Cov12=',num2str(cov_s_ica(1,2))]  ) )
grid on
axis equal

%% ex 5
clear all;
clc;
close all;
filname = {'Guitar1.wav', 'Guitar2.wav'}
n_sample = [2,100000];
S_sound = [];
S_dummy = [];
%%
j = 1;
y = [];
for i = 1:2:4
    
    S_dummy(:,i:i+1) = audioread(string(filname(j)), n_sample);
    j= j+1;
end

%S_sound(:,1) = 3*S_dummy(:,1) - S_dummy(:,2) + 2*S_dummy(:,3) + 5*S_dummy(:,4);
%S_sound(:,2) = 24*S_dummy(:,1) + 3*S_dummy(:,2) - 7*S_dummy(:,3) - 2*S_dummy(:,4);
S_sound(:,1) = 3*S_dummy(:,1) + 2*S_dummy(:,3) ;
S_sound(:,2) = 24*S_dummy(:,1) - 7*S_dummy(:,3)  ;

%%
%sound(S_sound(:,1))
%sound(S_sound(:,2))
%%

[ica_vect_sound, A_ica_sound, W_ica_sound] = fastica(S_sound');

S_ica_sound = A_ica_sound(:,1)*ica_vect_sound(1,:); % = A_ica*S*W_ica ----------> this to extract only respect one ica
S_ica_1_sound = A_ica_sound(:,1)*ica_vect_sound(1,:); % = A_ica*S*W_ica ----------> this to extract only respect one ica
S_ica_2_sound = A_ica_sound(:,2)*ica_vect_sound(2,:);


%sound1 = W_ica_sound(:,1) * S_sound;


%S_sounnd1 = A_ica_sound(:,1)*sound1;
%% ex 3-4
cov_ica = round(cov(score_sound'),2);

%%
plot_n = 1;
cov_s_sound = round(cov(S_sound ),2);

cov_ica_sound = round(cov(ica_vect_sound'),2);
cov_s_ica_sound = round(cov(S_ica_sound'),2);
%cov_s_pca = round(cov(S_pca'),2);
figure(plot_n)
subplot(1,3,1),plot(S_sound(:,1), S_sound(:,2),'.')
title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s_sound(1,1)),'; Var_2=',num2str(cov_s_sound(2,2)),'; Cov=',num2str(cov_s_sound(1,2))]  ) )
grid on
axis equal

subplot(1,3,2), plot(S_ica_sound(1,:),S_ica_sound(2,:),'.')
xlabel('ICA 1'),ylabel('ICA 2')
title(splitlines(['The score',newline,'Var_1=',num2str(cov_ica_sound(1,1)),'; Var_2=',num2str(cov_ica_sound(2,2)),'; Cov=',num2str(cov_ica_sound(1,2))]  ) )
grid on
axis equal

subplot(1,3,3), plot(S_ica_sound(1,:),S_ica_sound(2,:),'.b', S_ica_w_sound(1,:),S_ica_w_sound(2,:),'.r')
xlabel('reconstructed 1'),ylabel('reconstructed 2')
%% ex 5.a
plot_n = plot_n + 1
figure(plot_n)
subplot(3,2,1),plot(S_sound(:,1))
legend('S_sound 1')
grid on
%axis equal
subplot(3,2,2),plot(S_sound(:,2))
legend('S_sound 2')
legend('S_sound ica 1')
grid on
%axis equal
subplot(3,2,3),plot(ica_vect_sound(1,:))
hold on
plot(S_dummy(:,1))
grid on
%axis equal
subplot(3,2,4),plot(S_ica_sound(2,:))
legend('S_sound ica 2')
grid on
%axis equal
%% ex.5.b
sound(S_sound(:,1))
sound(ica_vect_sound(1,:))

%% ex.5.b
sound(S_sound(:,2))
sound(S_ica_sound(2,:))
