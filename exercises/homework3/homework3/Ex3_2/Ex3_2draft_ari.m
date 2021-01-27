%% EX1
%Laplacian distribution 
% random variable has a Laplace(mu, b) distribution:
% with mean mu and standard deviation sigma.
% mu : mean or location parameter
% b = scale parameter 
% [m, n] : the dimension of y.

%randlpl(mu, b, m, n)
%10000 samples from a Laplacian
% [1, 10000]
mu = 10
b = 3
m = 1
n = 10000

d1 = randlpl(mu, b, m, n);
d2 = randlpl(mu, b, m, n);

lc1 = 2*d1-d2;
lc2 = d1+0.5*d2;

figure(1)
%plot (lc1,lc2)
scatter (lc1,lc2,'filled')
xlabel ('lc1 = 2*d1-d2')
ylabel ('lc2 = d1+0.5*d2')
title ('1) Plotting linear combination of distributions against each other')

%% EX 2
C = cov (lc1, lc2)
[coeff,latent,explained] = pcacov(C)
figure (2)
biplot(coeff)
