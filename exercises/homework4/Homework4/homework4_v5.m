clear all;
clc;
close all;

%% 0-Local paths, update with you own, dont delete the others ;) 
%results_file_path = "C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/results.txt";
results_file_path = "/home/silvia-neurorobotics/Documents/cognitive_modeling/week11/Stimulus_presentation_script/results.txt";

%results_file_path ="C:/Users/arita/Desktop/homework_4/results.txt"
%myFolder = 'C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/pictures';
%myFolder = "C:/Users/arita/Desktop/homework_4/pictures";
myFolder = '/home/silvia-neurorobotics/Documents/cognitive_modeling/week11/Stimulus_presentation_script/pictures';

%% 1-Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["subject_id", "trial_number", "stimulus_filename", "answer", "reaction_time"];
opts.VariableTypes = ["double", "double", "string", "categorical", "double"];
opts = setvaropts(opts, 3, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [3, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
results = readtable(results_file_path, opts);
summary(results)

% Clear temporary variables
clear opts

%% 2-Clean the data (remove data below 200ms and above 2 sec)
toDelete = results.reaction_time < 200*10^(-3);
results(toDelete,:) = [];
toDelete = results.reaction_time >2;
results(toDelete,:) = [];

%% 3-From binary to continuos

% for each subject normalize the reaction time 
% quality of the answer from not very smiling to very smiling
user_id = [1,2,111089,20];
zero_size_results = zeros(size(results,1),1);
results.VarName6 = zero_size_results;
for i = 1:1:size(user_id,2)
    dummy_user_id = results.subject_id == user_id(i);
    results(dummy_user_id,:).VarName6 =  1-normalize(results(dummy_user_id,:).reaction_time, 'range', [0 1]);
end

toChngsign = results.answer == 'Not smiling';
results(toChngsign,6).VarName6 = - results(toChngsign,:).VarName6;



%% 4- from I to S
I_dim = 300;

filePattern = fullfile(myFolder, '*.jpg');
jpegFiles = dir(filePattern);
% S transformed matrix, 
S = [];
for k = 1:length(jpegFiles)
    baseFileName = jpegFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]); 
    S(k,:) = reshape(double(I),[1,I_dim*I_dim]);
    
end
%% 5-PCA 
% Principal Components
S_mean = mean(S,1);

[X,W,latent,tsqu,explained] = pca(S-S_mean);

% check how many eigenvalues are needed to have 95% variance
i = 0;
sum_variance = 0;
var_threshold = 95;
while sum_variance <= var_threshold
   i = i+1;
   sum_variance = sum_variance + explained(i);
end
n_pca = i;

%
sum_pc6 = sum(explained(1:6));
sum_pc118 = sum(explained(1:n_pca));

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
title([num2str(n),' pca (',num2str(sum_pc118),'% variance)']) 
grid on

n = 6;
subplot(1,3,3)
x = [1:n]; y = explained(1:n);
ylabel('Var %')
bar(x,y)
title([num2str(n),' pca (',num2str(sum_pc6),'% variance)']) 
grid on

plot_n = plot_n+1;

%% 5-b A LOOK AT IMAGES RECONSTRUCTED WITH PC = 118
% sub scores with that explains the 95%
W_sub = W(:,1:n_pca);
X_sub = X(:,1:n_pca);

S_pca = W_sub*X_sub' + S_mean;

num_images =4; % how many examples you want to plot
n = n_pca;
for k = 1:num_images %max is the length(jpegFiles) = 507
          
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);   
    
    %reconstruced image
    I_pca = reshape(uint8(S_pca(k,:)),[I_dim,I_dim]);
    
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title([num2str(n),' pca (',num2str(sum_pc118),'% variance)']) 
    grid on
end 

%%
PCA_reconstruct = @(W,X,i,j,m) W(:,i )*X(:,j)' + m;

%% Plot pca 

n_features = n_pca;
plot_n = 1
figure(plot_n) 
supbplot_n = 1;
step_for = 1;
sub_plt_dim = 8;
for n_c = 1:step_for:n_pca

    S_reconstruct = PCA_reconstruct(W_sub, X_sub,n_c,n_c,0);
 
    S_reconstruct_2 = reshape(double(S_reconstruct(k,:)),[ I_dim,I_dim]);
        
    if n_c == n_pca/2
        supbplot_n = 1
        plot_n = plot_n+ 1;
        figure(plot_n) 
    end
    subplot(sub_plt_dim,sub_plt_dim,supbplot_n) 
    imshow(S_reconstruct_2,[])   
    title(['pca ',num2str(n_c) ]) 
    supbplot_n = supbplot_n +1; 
end

%% 6- target vector y
statarray = grpstats(results,'stimulus_filename','mean','DataVars',{'VarName6'});
y   = statarray.mean_VarName6; %VarName6 is smiling/not smiling

%% 7- SELECT BEST FEATURES
c = 10;
opts = statset('Display','iter');
fun= @(xtrain, ytrain, xtest, ytest) sum((predict(fitlm(xtrain, ytrain),xtest) - ytest).^2);%<----------------this works, faster
[fs,history] = sequentialfs(fun,W_sub,y,'cv',c,'options',opts);

%% 8- Selected coefficients (X_sel) and Selected scores (W_sel)
W_sel = W_sub(:,fs(1,:) == 1);
X_sel = X_sub(:,fs(1,:) == 1);

S_pca = W_sel*X_sel'+S_mean;
e = S-S_pca;

index_fs = [];
dummy_idx = 1;
for i = 1:1:size(fs,2)
    if fs(1,i) == 1
        index_fs(dummy_idx) = i;
        dummy_idx = dummy_idx +1;
    end
end
%% Plot selected pca 
n_features = size(X_sel,2);
plot_n = 1
figure(plot_n) 
supbplot_n = 1;
step_for = 1;
sub_plt_dim = round(sqrt(size(1:step_for:n_features,2)));
for n_c = 1:step_for:n_features
    S_reconstruct = PCA_reconstruct(W_sel,X_sel,n_c,n_c,0);
    S_reconstruct_2 = reshape(double(S_reconstruct(k,:)),[ I_dim,I_dim]);
    subplot(sub_plt_dim,sub_plt_dim,supbplot_n) 
    imshow(S_reconstruct_2,[])  
    title(['pca ',num2str(index_fs(n_c)) ])
    supbplot_n = supbplot_n +1; 
end

%% 8-b A LOOK AT IMAGES RECONSTRUCTED WITH PC Selected

num_images =4; % how many examples you want to plot

for k = 1:num_images %max is the length(jpegFiles) = 507
          
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);   
    
    %reconstruced image
    I_pca = reshape(uint8(S_pca(k,:)),[I_dim,I_dim]);
    
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title([num2str(n_features),' selected components ']) 
    grid on
end 

%% 9- Linear regression model
model = fitlm(W_sel,y)
anova(model,'summary')

A = model.Coefficients.Estimate(2:end)'; %<--------- take the transpose to have row vector
b = model.Coefficients.Estimate(1); % intercept, b


%% 10- Lets get a face with a certain level of smileyness (y_0)

fig_idx = 1; 
plot_n = 1;    
figure(plot_n)      
for_step = 0.5;
dim_sub_plt = size(-1:for_step:1,2);
for y_0 = -1:for_step:1           %<------- selected score we want to create a face for (eg very smiley = +1) 
    alpha = (y_0-b)./(norm(A)^2); %<------- formulas from last slides 
    x = alpha * A;                %<------- formulas from last slides 
    S_new = (x*X_sel')+b + S_mean ;   %<------- formulas from Paolo
    I_new= reshape(uint8(S_new),[I_dim,I_dim]);   
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
