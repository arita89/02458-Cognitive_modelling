function p = plotnormal(x,m,s)
a = min(x); b = max(x);
d = a + (b-a) * rand(1, 500);
f = normpdf(x, m, s);
plot(x,f)
end