%%Homework1-part1
clear
clc
format short;  

% y_f = k*log(t./T)  % Perceived Brightness with Fechner law
% build y_f = A*x where 
% x(1) = k; x(3) = T,threshold
% y_f = x(1)*log(t)-x(1)*log(x(3))<--- x(1)*log(x(3)=x(2)
% y_f = x(1)*log(t)-x(2))
% A = [a1' a2']
%a1 = ones(1,max(t)).*log(t);
%a2 = ones(1,max(t)).*(-1);

t = 1:100;  
y = 10*(t).^0.33;% Perceived Brightness with Stevens law
A = [[ones(1,max(t)).*log(t)]' [ones(1,max(t)).*(-1)]'];
x = pinv(A)*y'; % <-- returns x(1)=k and x(2)
x(3) = exp(x(2)/x(1));
y_f = x(1)*log(t./x(3));% <-- y_f = k*log(t./T);

ye = 0.00015*(t).^3.3;%Perceived electric shock with Stevens law
B = [[ones(1,max(t)).*log(t)]' [ones(1,max(t)).*(-1)]'];
z = pinv(B)*ye'; % <-- returns z(1)=k and z(2)
z(3) = exp(z(2)/z(1));
ye_f = z(1)*log(t./z(3)); % <-- y_f = k*log(t./T);

%plot 
figure(1)
plot (t, y,'b-',t,y_f,'r-');
legend('Steven','Fechner')
xlabel('Luminance')
ylabel('Perceived brightness')
grid on

figure(2)
plot (t,ye,'b-',t,ye_f,'r-');
legend('Steven','Fechner')
xlabel('Electroshock magnitude')
ylabel('Perceived Electroshock')
grid on









