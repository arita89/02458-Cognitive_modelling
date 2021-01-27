I = imread('leafs.jpg');
figure (1), imagesc(I)
%X = reshape(I,size(I,1)*size(I,2),3);
%figure (2), imagesc(X)
imSz = size(I)
patchSz = [10 10];
xIdxs = [1:patchSz(2):imSz(2) imSz(2)+1];
yIdxs = [1:patchSz(1):imSz(1) imSz(1)+1];
patches = cell(length(yIdxs)-1,length(xIdxs)-1);
for i = 1:length(yIdxs)-1
    Isub = I(yIdxs(i):yIdxs(i+1)-1,:);
    for j = 1:length(xIdxs)-1
        patches{i,j} = Isub(:,xIdxs(j):xIdxs(j+1)-1);
    end
end
patches;
numel(patches{2,62})
figure, imagesc(patches{2,3})
