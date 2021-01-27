clc;
clear all;
close all;
%% Import image
fig_title = {'MonaLisaBW.jpg';'beach.jpg';'mountain.jpg';'Autumn.jpg'};
plot_n = 1;
I_dim = 1000;
size_patch = 10;
S_dim = I_dim/size_patch;

for figure_t = 1:1:size(fig_title,1)
    I = imread(string(fig_title(figure_t)));
    
    % if color image convert it to greyscale
    if size(I,3)== 3
        I = rgb2gray(I);
    end
    % make sure that all the pic have same size
    I = imresize(I,[I_dim I_dim]) ;
    
    % plot original picture
    fig_idx = 1;
    figure(plot_n), subplot(2,2,fig_idx), imshow(I)
    title('original')
    plot_n =plot_n + 1;
    fig_idx = fig_idx + 1;

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
    I_pca = W(:,1:n_pca)*X(:,1:n_pca)';

    %%
    I_project = [];
    i = 1; % iterate of the I_pca matrix row 

    for j = 1:size_patch:I_dim
        for k = 1:size_patch:I_dim
            I_project(j:j+size_patch-1,k:k+size_patch-1) =  reshape(I_pca(i,:), size_patch,size_patch) ;
            i = i+1;
        end
    end
    
    I_project = uint8(I_project);
    
    % plot
    subplot(2,2,fig_idx), imagesc(I_project)
    title('pca scene')
    fig_idx = fig_idx + 1;
    subplot(2,2,fig_idx), imshow(I_project)
    title(['not scaled ',num2str(n_pca),'pca (95% variance)'])
    fig_idx = fig_idx + 1;
    subplot(2,2,fig_idx), imshow(I_project,[])
    title(['scaled ',num2str(n_pca),' pca (95% variance)'])
    sgtitle(string(fig_title(figure_t)))
    fig_idx = fig_idx + 1;

    figure(plot_n), 
    subplot(1,2,1),bar(latent)
    title( 'pca variance' )
    grid on
    subplot(1,2,2),bar(latent(1:6))
    title( 'zoom pca variance' )
    sgtitle(string(fig_title(figure_t)))
    grid on
    plot_n = plot_n+1;
    
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
    subp_id = 1;
    for pca_idx_row = 1:1:6
        for pca_idx_col = 1:1:6
            subplot(6,6,subp_id),plot(W(:,pca_idx_row), W(:,pca_idx_col),'.')
            %title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
            xlabel('PCA '+string(pca_idx_row)),ylabel('PCA '+string(pca_idx_col))
            legend(num2str(cov_W(pca_idx_row,pca_idx_col)))
            grid on
            axis equal
            subp_id = subp_id + 1;
        end
    end
    sgtitle(string(fig_title(figure_t)))
    
    plot_n = plot_n+1;
end

%%
%close all