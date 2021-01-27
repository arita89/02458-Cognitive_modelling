clear all;
clc;
close all;

%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /home/silvia-neurorobotics/Documents/cognitive_modeling/week11/Stimulus_presentation_script/results.txt
%
% Auto-generated by MATLAB on 14-Nov-2019 10:36:43

%% Setup the Import Options
%results_file_path = "C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/results.txt";
results_file_path ="C:/Users/arita/Desktop/homework_4/results.txt"
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


%% Clean the data (remove data below 200ms and above 2 sec)
toDelete = results.reaction_time < 200*10^(-3);
results(toDelete,:) = [];
toDelete = results.reaction_time >2;
results(toDelete,:) = [];

%% From binary to continuos

% for each subject normalize the reaction time 
% quality of the answer from not very smiling to very smiling
user_id = [1,2];
zero_size_results = zeros(size(results,1),1);
results.VarName6 = zero_size_results;
for i = 1:1:size(user_id,2)
    dummy_user_id = results.subject_id == user_id(i);
    m = mean(results(dummy_user_id,:).reaction_time);
    results(dummy_user_id,:).VarName6 =  1-normalize(results(dummy_user_id,:).reaction_time, 'range', [0 1]);
end

toChngsign = results.answer == 'Not smiling';
results(toChngsign,6).VarName6 = - results(toChngsign,:).VarName6;



%%
I_dim = 300;

plot_n = 1;
I_dim = 300;
size_patch = 10;
S_dim = I_dim/size_patch;
%myFolder = 'C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/pictures';
myFolder = "C:/Users/arita/Desktop/homework_4/pictures";
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.jpg');
jpegFiles = dir(filePattern);
length(jpegFiles)
% S transformed matrix, 
S = [];
for k = 1:length(jpegFiles)
    baseFileName = jpegFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);

   
    S(k,:) = reshape(double(I),[1,I_dim*I_dim]);
    
end
%%
% Principal Components
[X,W,latent,tsqu,explained] = pca(S);

% check how many eigenvalues are needed to have 95% variance
i = 0;
sum_variance = 0;
var_threshold = 95;
while sum_variance <= var_threshold
   i = i+1;
   sum_variance = sum_variance + explained(i);
end
n_pca = i;
% sub scores with that explains the 95%
W_sub = W(:,1:n_pca);
%%
figure(plot_n), 
subplot(1,2,1),bar(latent)
title( 'pca variance' )
grid on
subplot(1,2,2),bar(latent(1:6))
title( 'zoom pca variance' )
%sgtitle(string(fig_title('figure_t')))
grid on
plot_n = plot_n+1;

%%
statarray = grpstats(results,'stimulus_filename','mean','DataVars',{'VarName6'});
y   = statarray.mean_VarName6; %VarName6 is smiling/not smiling
%% SELECT BEST FEATURES
%rmse = @(X,y) sqrt(immse(X,y));

%err = immse(X,Y) calculates the mean-squared error (MSE) between the arrays X and Y. 
%X and Y can be arrays of any dimension,but must be of the same size and class.


%%Ari attempt 1
%y = y_train;
%L = @(X,y) loss(fitrensemble(X,y),X,y);
%inmodel = sequentialfs(L, W_sub,y);


%%Split W_sub in train and test set , dummy set to test criterion
% Cross varidation (train: 70%, test: 30%)
%cv = cvpartition(size( W_sub,1),'HoldOut',0.3)
%idx = cv.test;
% Separate to training and test data
%X_Train_dummy = X(~idx,:);
%X_Test_dummy  = X(idx,:);
%y_Train_dummy = y(~idx,:);
%y_Test_dummy  = y(idx,:);

%immse(y_Test_dummy, predict(fitrensemble(X_Train_dummy,y_Train_dummy), X_Test_dummy))%<----------------this works
%https://se.mathworks.com/help/stats/fitrensemble.html

%%Ari attempt 2
%criterion = @(XTRAIN,ytrain,XTEST,ytest)  immse(y_Test, predict(fitrensemble(X_Train,y_Train), X_Test));
%inmodel = sequentialfs(criterion, W_sub,y);

%%Ari attempt 3
%c = 10;
%opts = statset('Display','iter');
%fun =@(X_Train,y_Train,X_Test,y_Test)loss(fitrensemble(X_Train,y_Train),X_Test,y_Test);%<----------------this works
%[fs,history] = sequentialfs(fun,W_sub,y,'cv',c,'options',opts)

%%Ari attempt 4
%c = 10;
%opts = statset('Display','iter');
%fun = @(X_Train,y_Train,X_Test,y_Test)immse(y_Test, predict(fitrensemble(X_Train,y_Train), X_Test));
%[fs,history] = sequentialfs(fun,W_sub,y,'cv',c,'options',opts)

%%Ari attempt 5
c = 10;
opts = statset('Display','iter');
fun= @(xtrain, ytrain, xtest, ytest) sum((predict(fitlm(xtrain, ytrain),xtest) - ytest).^2);;%<----------------this works, faster
[fs,history] = sequentialfs(fun,W_sub,y,'cv',c,'options',opts)


%% Selected coefficients (X_sel) and Selecte scores (W_sel)

W_sel = W_sub(:,fs(1,:) == 1);
X_sel = X(:,fs(1,:) == 1);

S_pca = W_sel*X_sel';
e = S-S_pca;
% S = W*X' +b
% S = S_pca +b;
%% Linear regression model
model = fitlm(W_sel,y)
anova(model,'summary')

A = model.Coefficients(2:end,1) % vector of the coefficients of predictors, A 
b = model.Coefficients (1,1) % intercept, b

%%
y_0 = -0.99; % selected score we want to create a face for (eg very smiley = +1)
S_mean = mean(S,1);
alpha = (y_0-b.Estimate)./(norm(A.Estimate,2)); 
x = alpha * A.Estimate;  
S_new = (x'*X_sel')+ S_mean;

I_mean= reshape(uint8(S_mean),[I_dim,I_dim]);  
I_new= reshape(uint8(S_new),[I_dim,I_dim]);   
          
% PLOTS 
% plot original picture and PC reconstructed  
plot_n = plot_n +1;    
figure(plot_n)      
fig_idx = 1; 
subplot(1,2,fig_idx), imshow(I_mean) 
title('mean')          
fig_idx = fig_idx + 1; 
subplot(1,2,fig_idx), imshow(I_new,[])  
title('new happy face')  
%imwrite(I_pca,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 
%ypred = predict(model,'newface')
%%

for k = 1:num_images%length(jpegFiles) 
    I_pca = reshape(uint8(S_new(k,:)),[I_dim,I_dim]);       
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);           
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title(['scaled ',num2str(n_pca(c)),'pca'])  
    imwrite(I_pca,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 
end 
%% 
S_pca = W_sel*X_sel';
 
%reconstruction of S, except for the scaling, in the pc-space 

plot_n = 2; 
num_images =4; 

for k = 1:num_images%length(jpegFiles) 
    I_pca = reshape(uint8(S_pca(k,:)),[I_dim,I_dim]);       
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);           
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title(['scaled ',num2str(n_pca(c)),' pca (95% variance)'])  
    imwrite(I_pca,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 
end 



%%
% we look at the first tot images reconstructed with 118 components 
%n_pca = 118  
%n_pca = 6 
%n_pca = n  
for c = 1: length(n_pca) 
    n_pca = 1
    S_n_pca = W(:,1:n_pca)*X(:,1:n_pca)'; 
end 
for k = 1:num_images%length(jpegFiles) 
    I_pca = reshape(uint8(S_n_pca(k,:)),[I_dim,I_dim]);       
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);           
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title(['scaled ',num2str(n_pca(c)),'pca'])  
    imwrite(I_pca,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 
end 

%%
%%
for c = 1: length(n_pca) 
    n_pca = 1
    S_mean = b; 
end 
for k = 1:num_images%length(jpegFiles) 
    I_pca = reshape(uint8(S_mean(k,:)),[I_dim,I_dim]);       
    %original image 
    baseFileName = jpegFiles(k).name; 
    fullFileName = fullfile(myFolder, baseFileName); 
    I = imresize(rgb2gray(imread(fullFileName)),[I_dim, I_dim]);           
    % PLOTS 
    % plot original picture and PC reconstructed  
    plot_n = plot_n +1;    
    figure(plot_n)      
    fig_idx = 1; 
    subplot(1,2,fig_idx), imshow(I) 
    title('original')          
    fig_idx = fig_idx + 1; 
    subplot(1,2,fig_idx), imshow(I_pca,[])  
    title(['scaled ',num2str(n_pca(c)),'pca'])  
    imwrite(I_pca,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 
end 

 




























