function [out] =chaoperm(im,pr,pc,num,forward)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[rows,cols] = size(im);
mat = zeros([rows,cols,num+1]);
mat(:,:,1) = im(:,:);
for loc=2:num+1
if(strcmp(forward,'forward'))
for i=1:rows
for j=1:cols
mat(pr(i,j),pc(i,j),loc) =mat(i,j,loc-1);
end
end

elseif(strcmp(forward,'backward'))
for i=1:rows
for j=1:cols
mat(i,j,loc) = mat(pr(i,j), pc(i,j),loc-1);
end
end
end
end
out = mat(:,:,num+1);

end
%*****************************************************

