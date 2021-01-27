fileID = 'C:\Users\arian\OneDrive\Documents\1semestre\02458 Cognitive Modelling\homework4\Stimulus presentation script\results';
formatSpec = '*.txt';
data = readtable(fileID, 'HeaderLines',1);
[M,N]= size(data)

for i = 1 : N
   for j = 1 : M
      val = A{i,j};
      if (ischar(val))
         val = str2num(val);
      end
      M(i,j) = val;
   end
end

time = data(:,5);
time = table2array(time)
[n_entries,a] = size(time);
mean(time)
max (time)
min (time)

data(time < 0.20) = [];
 
