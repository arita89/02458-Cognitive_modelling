clc;
clear all;
close all;
%% Import image
fig_title = {'woodBW.jpg'};
%fig_title = {'MonaLisaBW.jpg'};
I_dim = 1000;
size_patch = 10;
S_dim = I_dim/size_patch;

plot_n = 1;
n_raws = 2;
n_columns = 3;
fig_idx =0;

for figure_t = 1:1:size(fig_title,1)
    I = imread(string(fig_title(figure_t)));
    
    % if color image convert it to greyscale
    %if size(I,3)== 3
    %    I = rgb2gray(I);
    %end
    % make sure that all the pic have same size
    I = imresize(I,[I_dim I_dim]) ;
    
    % S transformed matrix, 
    % each row is a 10x10 patch converted into 1x100 vector
    S = [];
    k = 1;
    for i = 1:10:size(I,1)
        for j = 1:10:size(I,2)
            S(k,:) = reshape(double(I(i:i+size_patch-1,j:j+size_patch-1)),1,100 );
            k = k+1;
        end
    end

    % Principal Components
    [X,W,latent,tsquared, explained, mu] = pca(S);

    % check how many eigenvalues are needed to have 95% variance
    i = 0;
    sum_variance = 0;
    var_threshold = sum(latent)*0.95;
    while sum_variance <= var_threshold
       i = i+1;
       sum_variance = sum_variance + latent(i)
    end
    n_pca = i;
    
    %sum(latent(1:21))/sum(latent)
    %sum(explained(1:21))/sum(explained)

    
    % reconstruction of S, except for the scaling, in the pc-space
    I_pca   = W(:,1:n_pca)*X(:,1:n_pca)';
    I_pcas = {};
    for pca_idx = 1:1:6
        I_pcas{pca_idx} = W(:,pca_idx)*X(:,pca_idx)';
    end
    
    I_reconstructed = [];
    I_r = {};
    i = 1; % iterate of the I_pca matrix row 

    for j = 1:size_patch:I_dim
        for k = 1:size_patch:I_dim
            I_reconstructed(j:j+size_patch-1,k:k+size_patch-1) =  reshape(I_pca(i,:), size_patch,size_patch) ;
       
            i = i+1;
        end
    end
    
    
    for pca_idx = 1:1:6
        i = 1; 
        for j = 1:size_patch:I_dim
            for k = 1:size_patch:I_dim
                I_r{pca_idx}(j:j+size_patch-1,k:k+size_patch-1) =  reshape(I_pcas{pca_idx}(i,:), size_patch,size_patch) ;
                i = i+1;
            end
        end
    end
    
    
    cov_s = round(cov(S ));
    cov_W = round(cov(W ));
    
    figure(plot_n)
    subplot(1,2,1),plot(S(:,1), S(:,2),'.')
    title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
    grid on
    axis equal
    subplot(1,2,2),plot(W(:,1), W(:,2),'.')
    xlabel('PCA 1'),ylabel('PCA 2')
    title(splitlines(['The score',newline,'Var_1=',num2str(cov_W(1,1)),'; Var_2=',num2str(cov_W(2,2)),'; Cov=',num2str(cov_W(1,2))]  ) )
    grid on
    axis equal
    sgtitle(string(fig_title(figure_t)))
    
    plot_n = plot_n+1;
    figure(plot_n)
    
    for subp_id = 1:1:6
        subplot(3,2,subp_id),imshow(I_r{subp_id}(1:size_patch,1:size_patch),[])
        title(['pca ',num2str(subp_id)])
        grid on,axis equal
    end
    
    
    plot_n = plot_n+1;
end

%-------------------------------------