%%Exam 2019 P4- Faces and PCA

clear all;
clc;
close all;
dash = ('----------------------------------');
%% 0-Local paths
myFolder = " C:/Users/arian/OneDrive/Documents/1semestre/02458 Cognitive Modelling/exam_scripts/faces";

a = readtable("C:/Users/arian/OneDrive/Documents/1semestre/02458 Cognitive Modelling/exam_scripts/smile_intensity.txt");

%% 1- Load images
filePattern = fullfile(myFolder, '*.jpg');
jpegFiles = dir(filePattern);

%% from I to S
% S transformed matrix
tot_faces = 400; %length(jpegFiles)
I_dim_r = 260;
I_dim_c = 360; %resolution 260x360


S = []; %S [400*93600]
for k = 1:length(jpegFiles) %400
    baseFileName = jpegFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim_r, I_dim_c]); 
    S(k,:) = reshape(double(I),[1,I_dim_r*I_dim_c]);
end

%% PCA 
% Principal Components
S_mean = mean(S,1);

[X,W,latent,tsqu,explained] = pca(S-S_mean);

% check how many eigenvalues are needed to have 95% variance
i = 0;
sum_variance = 0;
var_threshold = 90;
while sum_variance <= var_threshold
   i = i+1;
   sum_variance = sum_variance + explained(i);
end
n_pca = i;

%
sum_pc_6c = sum(explained(1:6));
sum_pc_90 = sum(explained(1:n_pca));

plot_n = 1;
% Plot PCA 
figure(plot_n), 
subplot(1,3,1)
bar(explained)
title( 'pca - variance' )
grid on

n = n_pca;
subplot(1,3,2)
x = [1:n]; y = explained(1:n);
ylabel('Var %')
bar(x,y)
title([num2str(n),' pca (',num2str(sum_pc_90),'% variance)']) 
grid on

n = 6;
subplot(1,3,3)
x = [1:n]; y = explained(1:n);
ylabel('Var %')
bar(x,y)
title([num2str(n),' pca (',num2str(sum_pc_6c),'% variance)']) 
for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
grid on

plot_n = plot_n+1;

%% draw first 6 components as images
% gives size error that i hadnt time to debug

% n_pca = 6;
% I_pca   = W(:,1:n_pca)*X(:,1:n_pca)';
% I_pcas = {};
% for pca_idx = 1:1:6
%     I_pcas{pca_idx} = W(:,pca_idx)*X(:,pca_idx)';
% end
% 
% I_reconstructed = [];
% I_r = {};
% i = 1; % iterate of the I_pca matrix row 
% 
% for j = 1:I_dim_r
%     for k = 1:I_dim_c
%         I_reconstructed(j,k) =  reshape(I_pca(i,:), I_dim_r,I_dim_c) ;
%         i = i+1;
%     end
% end
% 
% 
% for pca_idx = 1:1:6
%     i = 1; 
%     for j = 1:I_dim_r
%         for k = 1:I_dim_c
%             I_r{pca_idx}(j,k) =  reshape(I_pcas{pca_idx}(i,:), I_dim_r,I_dim_c) ;
%             i = i+1;
%         end
%     end
% end
% 
% 
% cov_s = round(cov(S ));
% cov_W = round(cov(W ));
% 
% figure(plot_n)
% subplot(1,2,1),plot(S(:,1), S(:,2),'.')
% title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
% grid on
% axis equal
% subplot(1,2,2),plot(W(:,1), W(:,2),'.')
% xlabel('PCA 1'),ylabel('PCA 2')
% title(splitlines(['The score',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
% grid on
% axis equal
% sgtitle(string(fig_title(figure_t)))
% 
% plot_n = plot_n+1;
% figure(plot_n)
% 
% for subp_id = 1:1:6
%     subplot(3,2,subp_id),imshow(I_r{subp_id}(1:I_dim_r,1:I_dim_c,[]))
%     title(['pca ',num2str(subp_id)])
%     grid on,axis equal
% end
% 
% 
% plot_n = plot_n+1;

%% Linear regression model f
n_sel = 20; % first 20 components

W_sel = W (:,1:n_sel); % scores
X_sel = X (:,1:n_sel); % coeff

y = a.smile_intensity; % target vector, all smiles_intensities

model = fitlm(W_sel,y)
anova(model,'summary')

A = model.Coefficients.Estimate(2:end)'; %<--------- take the transpose to have row vector
b = model.Coefficients.Estimate(1); % intercept, b

%% Features Selection
c = 10;
opts = statset('Display','iter');
fun= @(xtrain, ytrain, xtest, ytest) sum((predict(fitlm(xtrain, ytrain),xtest) - ytest).^2);
[fs,history] = sequentialfs(fun,W_sel,y,'cv',c,'options',opts);

%% 8- Selected coefficients (X_sel) and Selected scores (W_sel)
W_best = W_sel(:,fs(1,:) == 1);
X_best = X_sel(:,fs(1,:) == 1);

S_pca = W_best*X_best'+S_mean;
e = S-S_pca;

index_fs = [];
dummy_idx = 1;
for i = 1:1:size(fs,2)
    if fs(1,i) == 1
        index_fs(dummy_idx) = i;
        dummy_idx = dummy_idx +1;
    end
end

%% generate new faces
model = fitlm(W_best,y)
anova(model,'summary')

A = model.Coefficients.Estimate(2:end)'; %<--------- take the transpose to have row vector
b = model.Coefficients.Estimate(1); % intercept, b


plot_n = plot_n+ 1;    
figure(plot_n)      
fig_idx = 1; 
for_step = 0.5;
dim_sub_plt = size(-0.5:for_step:1.5,2);
for y_0 = -1:for_step:1           %<------- selected score we want to create a face for (eg very smiley = +1) 
    alpha = (y_0-b)./(norm(A)^2);  
    x = alpha * A;               
    S_new = (x*X_best')+b + S_mean ;  
    I_new= reshape(uint8(S_new),[I_dim_r,I_dim_c]);   
    subplot(1,dim_sub_plt,fig_idx), imshow(I_new, [])  
    if y_0 == -1
        title('Not Smiling')  
    elseif y_0 == 0
        title('Neutral')  
    
    elseif y_0 == 1
        title('Smiling')  
    else
        title([num2str(y_0*100),'%']) 
    end
    fig_idx = fig_idx + 1; 
end
