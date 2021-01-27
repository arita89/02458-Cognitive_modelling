
%% FOR cicle

for x = -1:0.5:1 %start:step:end             
    
    if x == -1
        disp('first option')  
    elseif x == 1
        disp('last option')   
    else
        disp([num2str(x*100),'%']) 
    end
 
end