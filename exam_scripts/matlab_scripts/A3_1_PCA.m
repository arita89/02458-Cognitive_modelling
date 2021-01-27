clc;
clear all;
close all;
%% Import image
fig_title = {'beach.jpg'};
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
    
    %sum(latent(1:21))/sum(latent)
    %sum(explained(1:21))/sum(explained)

    % reconstruction of S, except for the scaling, in the pc-space
    n_pca = [100,i,6];
    plot_n = 1;
    for c = 1: length(n_pca)
        c;
        I_pca = W(:,1:n_pca(c))*X(:,1:n_pca(c))';
        I_project = [];
        i = 1; % iterate of the I_pca matrix row 

        for j = 1:size_patch:I_dim
            for k = 1:size_patch:I_dim
                I_project(j:j+size_patch-1,k:k+size_patch-1) =  reshape(I_pca(i,:), size_patch,size_patch) ;
                i = i+1;
            end
        end
    
        %I_project = uint8(I_project); #?
        %I_project = uint16(I_project);
        
        % PLOTS
        % plot original picture and PC reconstructed without scaling
        fig_idx = 1;
        plot_n = 1 +(c-1);
        figure(plot_n)
        subplot(1,2,fig_idx), imshow(I)
        title('original')
        
        fig_idx = fig_idx + 1;
        subplot(1,2,fig_idx)
        I_project = imresize(I_project,[I_dim I_dim]) ;
        imagesc(I_project);
        colormap('gray');
        title('pca scene- imagesc')
        
        sgtitle(['Reconstruction of Image with',num2str(n_pca(c)),' principal components'])
        imwrite(I_project,sprintf('%d.jpg',plot_n))% save a dedicated .jpg file
        
        %increment plot_n of one and re-set fig_idx
        plot_n =plot_n + 1 +(c-1);
        fig_idx = 1;
        
        %plot original and PC reconstructed with scaling
        figure(plot_n)
        subplot(1,2,fig_idx), imshow(I_project,[]) 
        colormap('gray');
        title(['scaled ',num2str(n_pca(c)),' pca (95% variance)']) 
        %sgtitle(string(fig_title(figure_t)))
        imwrite(I_project,sprintf('%d.jpg',plot_n))
        
        fig_idx = fig_idx + 1;
        subplot(1,2,fig_idx), imshow(I_project)
        title(['not scaled ',num2str(n_pca(c)),'pca (95% variance)'])
   
    end
    %I_project = uint16(I_project);
    %imwrite(I_project,'mod.jpg')
%------------------------------------------------------------------
%n = 100pc
    plot_n =plot_n + 1;
    figure(plot_n)
    
    fig_idx = 1;
    subplot(2,1,fig_idx)
    n = 100;
    x = [1:n]; y = explained(1:n);
    ylabel('Var %')
    bar(x,y)
    title(['% Variance over',num2str(n),' principal components'])
    grid on
    
    fig_idx = fig_idx + 1;
    subplot(2,1,fig_idx)
    n = 100;
    x = [1:n]; y = latent(1:n);
    ylabel('Var')
    bar(x,y)
    title(['% Variance over',num2str(n),' principal components'])
    grid on

    sgtitle(string(fig_title(figure_t)))
 %--------------------------------------------------------   
 % n= 21pc   
    plot_n =plot_n + 1;
    figure(plot_n)
    fig_idx = 1;
    subplot(2,1,fig_idx)
    n = 21;
    x = [1:n]; y = explained(1:n);
    bar(x,y)
    title(['Variance over ',num2str(n),' principal components'])
    ylabel('Var %')
    for i1=1:numel(y)
        text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
    end  
    grid on   
 
 
    fig_idx = fig_idx + 1;
    subplot(2,1,fig_idx)
    a = cumsum(latent)/sum(latent);
    b = cumsum(explained);
    x = [1:21]; y = b(1:21);
    bar(x,y)
    title(['Cumulative variance over ',num2str(n),' principal components'])
    xlabel('First 21 PCs')
    ylabel('Var %')
    for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
    end
    sgtitle(string(fig_title(figure_t)))
%-----------------------------------------------------------------
%n =6pc
    plot_n =plot_n + 1;
    figure(plot_n)
    fig_idx = 1;
    subplot(2,1,fig_idx)
    n = 6;
    x = [1:n]; y = explained(1:n);
    bar(x,y)
    
    title(['Variance over ',num2str(n),' principal components'])
    ylabel('Var %')
    name = {'pc1';'pc2';'pc3';'pc4';'pc5';'pc6'}; 
    set(gca,'xticklabel',name)
    for i1=1:numel(y)
        text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
    end  
    grid on
    
    fig_idx = fig_idx + 1;
    subplot(2,1,fig_idx)
    a = cumsum(latent)/sum(latent);
    b = cumsum(explained);
    x = [1:n]; y = b(1:n);
    bar(x,y)
    title(['Cumulative variance over ',num2str(n),' principal components'])
    ylabel('Var %')
    name = {'pc1';'pc2';'pc3';'pc4';'pc5';'pc6'}; 
    set(gca,'xticklabel',name)
    sgtitle(string(fig_title(figure_t)))
    for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
    end
    grid on
 %--------------------------------------------------------   
 
    cov_s = round(cov(S ));
    cov_W = round(cov(W ));
    plot_n = plot_n+1;
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
    
%%--------------------------------------------------------   
    %plot_n = plot_n+1;   
    %figure(plot_n)
    %subp_id = 1;
    %for pca_idx_row = 1:1:6
        %for pca_idx_col = 1:1:6
            %subplot(6,6,subp_id),plot(W(:,pca_idx_row), W(:,pca_idx_col),'.')
            %title(splitlines(['Original data covariance',newline,'Var_1=',num2str(cov_s(1,1)),'; Var_2=',num2str(cov_s(2,2)),'; Cov=',num2str(cov_s(1,2))]  ) )
            %xlabel('PCA '+string(pca_idx_row)),ylabel('PCA '+string(pca_idx_col))
            %legend(num2str(cov_W(pca_idx_row,pca_idx_col)))
            %grid on
            %axis equal
            %subp_id = subp_id + 1;
        %end
    %end
    %sgtitle(string(fig_title(figure_t)))

end

%%
%close all