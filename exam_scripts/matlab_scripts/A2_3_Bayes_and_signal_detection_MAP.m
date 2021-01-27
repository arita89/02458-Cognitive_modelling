%% assignment 2 part 3 signal detection and max a posteriori rule SEE ex4_bayesian_observer2016
clc;
clear all;
close all;

prior_signal = [0.5,0.95,0.15];
d = 1.5;
c=[];
for i=1:1:size(prior_signal,2)
a = log(1-prior_signal(i))-log(prior_signal(i));
c(i) = (a*d^2)/(2*d);
end
display('criterion:')
display(c)
