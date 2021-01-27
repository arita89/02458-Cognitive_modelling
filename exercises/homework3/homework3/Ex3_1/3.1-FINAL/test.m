i = 15;   
n_pca = [100,i,6];
plot_n =  1  
    for c = 1: length(n_pca)
        plot_n = plot_n+1;
        A = [c,n_pca(c),plot_n];
        disp(A)
    end