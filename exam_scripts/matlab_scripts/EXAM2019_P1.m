%% P1  strong fusion model
clear
clc
format short;  
dash = ('----------------------------------');
%%
%xA = load('xAlapse.txt');
x_v = load('xV.txt');
x_av = load('xAV.txt');

%% QUESTION
% xa = ? 
% mu_a = ?
% sigma_a = ?
% how does the observer reacts when there is only the audio stimulus
% no visual stimulus
%%
real_m_v    = mean(x_v)
real_s_v    = std(x_v)
real_m_av   = mean(x_av) % will be 3 values
real_s_av   = std(x_av) % will be 3 values
%% sigma_a and w

b = (real_s_av).^2 % will be 3 values
a = (real_s_v).^2
for i = 1:3
    s_a(i) = ((b(i).*a)/(a-b(i))).^(0.5);
    w(i) = (real_s_v).^2/(s_a(i).^2+(real_s_v).^2);
    m_a(i) = (real_m_av(i)- (1-w(i))*real_m_v)/w(i);
end

for i = 1:3
    disp (['audio from location number ',num2str(i)])
    disp (['standard deviation: ',num2str(s_a(i))])
    disp (['mean: ',num2str(m_a(i))])
    disp (dash)
end



