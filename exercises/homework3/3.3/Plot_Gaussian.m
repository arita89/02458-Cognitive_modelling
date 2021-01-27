function p = plotnormal(x,m,s)
a = min(x); b = max(x);
d = a + (b-a) * rand(1, 500);
f = gauss_distribution(d, m, s);
plot(x,f,'.')
%grid on
%title('Bell Curve')
xlabel('Randomly produced numbers')
ylabel('Gauss Distribution') 

end