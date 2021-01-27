%J = imread('Autumn.jpg');
J = imread('MonaLisaBW.jpg');
%J = rgb2gray(J);
%figure, imagesc(J)
size (J)
I = imcrop(J,[1 1 799 1139]); %resizing image to shape 400x700
%[rows, columns, numberOfColorChannels] = size(I)

figure (1)
%subplot(2,1,1),imshow(RGB)
subplot(2,1,1),imshow(J)
title('original')
subplot(2,1,2),imshow(I)
title('grey scale and cropped')


imSz = size(I)
patchSz = [10 10];
xIdxs = [1:patchSz(2):imSz(2) imSz(2)+1]; % +10, 1:401
yIdxs = [1:patchSz(1):imSz(1) imSz(1)+1]; % +10, 1:701
r = length(yIdxs)-1 %40
c = length(xIdxs)-1 %70
patches = cell(r,c);
S = zeros(r*c,patchSz(1)*patchSz(2)); %initialize empty matrix (2800 rows, 100 columns)
k = 1; %counter for matrix lines

for i = 1:r
    Isub = I(yIdxs(i):yIdxs(i+1)-1,:);
    for j = 1:c
        patches{i,j} = Isub(:,xIdxs(j):xIdxs(j+1)-1); %creating patches 
        sz = [1,numel(patches{i,j})]; %raw vector = array [1,100] or 
        %sz = [numel(patches{i,j}),1]; %column vector = array [100,1] ? 
        vectors{i,j} = reshape(patches{i,j},sz);
        S(k,:) = vectors{i,j}; %filling initialized matrix with vectors
        k = k+1;
    end        
end
%k
%patches
%vectors
%celldisp(patches(40,70))
%celldisp(vectors(40,70))


%P_test = patches{2,3}; check how a patch gets vectorized
%V_test =vectors{2,3}; 
 
[coeff, score, latent, tsquared, explained, mu] =  pca(S); % X= coefficients; W = scores
[U,T,V] = svd(S);
covarianceMatrix = cov(S);
[E,D] = eig(covarianceMatrix);
%eigval = latent;
%corr(W) --> orthogonality check
%%verifying that XW.', is a reconstruction of S, except for the scaling.
%size(score * coeff');
score = score(:,1:6); %select this two lines to see what happens with only the first 6 components
coeff = coeff(:,1:6);
S_hat = score*coeff'; %X*W'
S2_hat = U*T*V';

%figure (10)
%imshow(S2_hat)


%S_hat = reshape(S_hat, r,c);
%S_hat = X*W' + repmat(mu1,100,28);
vectors_hat=cell(1,100);
patches_hat=cell(r,c);
I_reconstruct = zeros(size(I));
I2_reconstruct = zeros(size(I));
vectors2_hat=cell(1,100);
patches2_hat=cell(r,c);
%k = 2801;
t = k;
for i = r:-1:1 %40 to 1
    for j = c:-1:1  %1 to 70
        if t > 0
            vectors_hat = S_hat(t-1,:);
            vectors2_hat = S2_hat(t-1,:);
            patches_hat{i,j} = reshape(vectors_hat, [10,10]);
            %patches_hat{i,j} = reshape(vectors_hat, [1,10]);
            patches2_hat{i,j} = reshape(vectors2_hat, [10,10]);%re-creating patches
            t = t-1;
        end
    end
end
I_reconstruct = cell2mat(patches_hat);
I2_reconstruct = cell2mat(patches2_hat);

figure (2)
subplot(1,3,1),imagesc(I)
title('original')
subplot(1,3,2),imagesc(I2_reconstruct)
title('SVD reconstructed') 
subplot(1,3,3),imagesc(I_reconstruct)
title('PCA reconstructed') %which does NOT contain the same amount of info...
%%
%variance explained by principal components
%two ways to get it, same stuff b = a*100
a = cumsum(latent)/sum(latent); %percentual
b = cumsum(explained); %decimals

figure(3)
Y = explained;
bar(Y)
%text([1:length(Y)],Y', num2str(Y','%0.2f'),'HorizontalAlignment','center','VerticalAlignment','bottom')
box off


figure (4)
x = [1:50]; y = a(1:50); 
bar(x,y)
title('Cumulative variance explained')
xlabel('principal components 1-50')
ylabel('Cum Var')
for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end    

figure (5)
x = [1:6]; y = explained(1:6);
bar(x,y)
title('6 principal components')
ylabel('Var')
name = {'pc1';'pc2';'pc3';'pc4';'pc5';'pc6'}; 
set(gca,'xticklabel',name)
for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
    
