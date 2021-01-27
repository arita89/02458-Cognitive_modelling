%%
X = coeff(1:6,:);
W = score(1:6, :);
eigval = latent;
%corr(W) --> orthogonality check
%%verifying that XW.', is a reconstruction of S, except for the scaling.
size(W * X')
S2_hat = W*X';


%S_hat = reshape(S_hat, r,c);
%S_hat = X*W' + repmat(mu1,100,28);
v2_hat=cell(1,100);
p2_hat=cell(r,c);
I2_reconstruct = zeros(size(I));

t2 = 6;
for i = 6:-1:1 %40 to 1
    for j = c:-1:1  %70 to 1
        if t > 0
            v2_hat = S2_hat(t2-1,:);
            p2_hat{i,j} = reshape(v2_hat, [10,10]); %re-creating patches
            t2 = t2-1;
        end
    end
end
I2_reconstruct = cell2mat(p2_hat);

figure (6)
subplot(2,1,1),imshow(I)
title('original')
subplot(2,1,2),imshow(I_reconstruct)
title('PCA2 reconstructed with 6 PC') %which does NOT contain the same amount of info...
