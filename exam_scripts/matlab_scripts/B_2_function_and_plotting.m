%% function holder
fun1= @(a,b) a.^b %esponenziale
fun2= @(a,b) log(a)./b

%% Plot+ Title, with two subplots and title

plot_n = 0; %initialization
fig_idx = 0;
n_raws = 3;
n_columns = 1;

plot_n =plot_n + 1;
figure(plot_n)
    
fig_idx = fig_idx + 1; 
subplot(n_raws,n_columns,fig_idx)
n = 100; % max len
x = [1:1:n]; %x axis 
y = fun1(2,1:n); %
ylabel('label y')
scatter(x,y,25,'b','o')
%bar(x,y) % uncomment for histogram
title(['first title with a number ',num2str(n),' in the title'])
grid on

fig_idx = fig_idx + 1; 
subplot(n_raws,n_columns,fig_idx)
n = -100;
t = [n:1:0]; 
xlabel('label t')
y1 = fun1(2,n:0);
y2 = fun2(2,n:0); 
plot (t,y1,'b-',t,y2,'r-'); %plot both on the same graph
legend('fun1','fun2')
title('plain title')
grid on

fig_idx = fig_idx + 1; 
subplot(n_raws,n_columns,fig_idx)
n = 100; % max len
start =1;
step = 10;
x = [start:step:n]; %x axis 
y = cumsum(1:step:n); %
xlabel('label x')
ylabel('label y')
bar(x,y) % uncomment for histogram
title(['histogram with a number ',num2str(step),' just to...'])

for idx=1:numel(y)%count element in array
text(x(idx),y(idx),num2str(y(idx),'%0.2f'),...
           'HorizontalAlignment','center',...
           'VerticalAlignment','bottom')
end
grid on

sgtitle('overall title')