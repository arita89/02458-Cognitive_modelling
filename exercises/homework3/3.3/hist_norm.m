function p = hist_norm(x,m,s)
n1= normpdf (x,m,s);
histfit(n1,100,'Normal')
%grid on
%title('Bell Curve')
end