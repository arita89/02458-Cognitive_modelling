xA = importdata('xA.txt');
xV = importdata('xV.txt');
xAV = importdata('xAV.txt');

mA=mean(xA);
sA=std(xA);

mV=mean(xV);
sV=std(xV);

%mAV=mean(xAV);
%sAV=std(xAV);

%w = ((s2^2)/((s1^2)+(s2^2)))

%fusion = @(xAV,m1,m2,s1,s2) ((mAV-(((s2^2)/((s1^2)+(s2^2)))*m1+(1-((s2^2)/((s1^2)+(s2^2))))*m2)^2)/((s1^2)*(s2^2))/((s1^2)+(s2^2)))
fusion = @(xAV,n,m1,m2,s1,s2) (-n/2)*log(2*phi)+(-n/2)*log(((s1^2)*(s2^2))/((s1^2)+(s2^2)))-(1/(2*(((s1^2)*(s2^2))/((s1^2)+(s2^2)))))*sum((xAV-((((s2^2)/((s1^2)+(s2^2)))*m1+(1-((s2^2)/((s1^2)+(s2^2))))*m2))).^2)

x0 = [xAV,50,1,1,1,1];
x = fminsearch(fusion,x0)