clear all;
clc;
close all;

%% 0-Local paths, update with you own, dont delete the others ;) 
results_file_path = "C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/results.txt";
%results_file_path ="C:/Users/arita/Desktop/homework_4/results.txt"
myFolder = 'C:/Users/arita/Documents/private/Data Engineering/1semestre/02458 Cognitive Modelling/homework4/pictures';
%myFolder = "C:/Users/arita/Desktop/homework_4/pictures";

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


plot_n = 0;

%% 2-Clean the data (remove data below 200ms and above 2 sec)
toDelete = results.reaction_time < 200*10^(-3);
results(toDelete,:) = [];
toDelete = results.reaction_time >2;
results(toDelete,:) = [];

%% 3-From binary to continuos

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



%% 4- from I to S
I_dim = 300;

plot_n = plot_n +1;
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
%% 5-PCA 
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
%
sum_pc6 = sum(explained(1:6));
sum_pc118 = sum(explained(1:118));

% Plot PCA 
figure(plot_n), 
subplot(1,3,1)
bar(explained)
title( 'pca - variance' )
grid on

n = 118;
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
W_118 = W(:,1:118);
X_118 = X(:,1:118);
S_pca = W_118*X_118';

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
    title([num2str(n),' pca (',num2str(sum_pc118),'% variance)']) 
    grid on
end 

%% 6- target vector y
statarray = grpstats(results,'stimulus_filename','mean','DataVars',{'VarName6'});
y   = statarray.mean_VarName6; %VarName6 is smiling/not smiling

%% 7- SELECT BEST FEATURES
c = 10;
opts = statset('Display','iter');
fun= @(xtrain, ytrain, xtest, ytest) sum((predict(fitlm(xtrain, ytrain),xtest) - ytest).^2);;%<----------------this works, faster
[fs,history] = sequentialfs(fun,W_sub,y,'cv',c,'options',opts)

%% 8- Selected coefficients (X_sel) and Selected scores (W_sel)
W_sel = W_sub(:,fs(1,:) == 1);
X_sel = X(:,fs(1,:) == 1);

S_pca = W_sel*X_sel';
e = S-S_pca;

%% 8-b A LOOK AT IMAGES RECONSTRUCTED WITH PC Selected
S_pca = W_sel*X_sel';
s = size(W_sel)
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
    title([num2str(s(2)),' selected components ']) 
    grid on
end 

%% 9- Linear regression model
model = fitlm(W_sel,y)
anova(model,'summary')

A = model.Coefficients(2:end,1); % vector of the coefficients of predictors, A 
A = (A.Estimate)'; %<--------- take the transpose to have row vector
b = model.Coefficients(1,1); % intercept, b
b = (b.Estimate)'; %<-----------take the transpose cause i am silly and i transpose numbers


%% 10- Lets get a face with a certain level of smileyness (y_0)

y_0 = 0.99; % selected score we want to create a face for (eg very smiley = +1)
S_mean = mean(S,1);
alpha = (y_0-b)./(norm(A,2)); %<------- formulas from last slides 
x = alpha * A;                %<------- formulas from last slides 
S_new = (x*X_sel')+ S_mean;   %<------- formulas from Paolo

I_mean= reshape(uint8(S_mean),[I_dim,I_dim]);  
I_new= reshape(uint8(S_new),[I_dim,I_dim]);   
          
% PLOTS 
% plot general mean and new image 
plot_n = plot_n +1;    
figure(plot_n)      
fig_idx = 1; 
subplot(1,2,fig_idx), imshow(I_mean) 
title('mean face')          
fig_idx = fig_idx + 1; 
subplot(1,2,fig_idx), imshow(I_new, [])  
title('new happy face')  

%imwrite(I_new,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file 


 




























