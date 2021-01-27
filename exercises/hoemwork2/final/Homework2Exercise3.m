prob_signal = 0.5;
% prob_signal = 0.95;
% prob_signal = 0.15;

x0 = 0;
criterion = @(x) ((normpdf(x - 1.5)/normpdf(x) - prob_signal/(1 - prob_signal))^2);

fminunc(criterion ,x0)