function [pr,pc] =chaomat(n)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%*****************************************************
%
I=sum(n);
k=size(n,2);
for i=1:k
N(i+1)=1;
for j=1:i
N(i+1)=N(i+1)+n(j);
end
end
N(1)=1;
%N(0)=1;
for cb=1:k
for rb=1:n(cb)
rbstartcol(rb)=mod((rb-1)*I,n(cb));
rbendcol(rb)=mod((rb*I-1),n(cb));
rbstartrow(rb)=fix(((rb-1)*I)/n(cb));
rbendrow(rb)=fix((rb*I-1)/n(cb));
mincol(rb)=min([rbendcol(rb)+1,rbstartcol(rb)]);
maxcol(rb)=max([rbendcol(rb),rbstartcol(rb)-1]);
end
for i=1:I
for j=N(cb):N(cb+1)-1
newindex(i,j-N(cb)+1)=(i-1)*n(cb)+(n(cb)-j+N(cb)-1);
newindexmod(i,j-N(cb)+1)=mod(newindex(i,j-N(cb)+1),n(cb));
newindexquotient(i,j-N(cb)+1)=fix(newindex(i,j-N(cb)+1)/n(cb));
rowblockindex(i,j-N(cb)+1)=fix(newindex(i,j-N(cb)+1)/I)+1;
end
end
for i=1:I
for j=1:n(cb)
for rb=1:n(cb)
if rowblockindex(i,j)==rb;
if newindexmod(i,j)>maxcol(rb)
col=rbendrow(rb)-newindexquotient(i,j)+(n(cb)-1-newindexmod(i,j))*(rbendrow(rb)-rbstartrow(rb));

elseif newindexmod(i,j)>=mincol(rb) && newindexmod(i,j)<=maxcol(rb)
if rbstartcol(rb)>rbendcol(rb)
c=0;
d=-1;
else
c=1;
d=1;
end
col=(rbendrow(rb)- rbstartrow(rb))*(n(cb)-1-maxcol(rb))+(rbendrow(rb)- newindexquotient(i,j)+c)+(maxcol(rb)-newindexmod(i,j))*(rbendrow(rb)-rbstartrow(rb)+d);
else %if newindexmod(i,j)<=mincol(rb)
col=I-mincol(rb)*(rbendrow(rb)-rbstartrow(rb))+(rbendrow(rb)-newindexquotient(i,j)+1)+(mincol(rb)-1-newindexmod(i,j))*(rbendrow(rb)-rbstartrow(rb));
end
row=1+I-N(cb+1)+rowblockindex(i,j);
end
end
pr(i,j+N(cb)-1)=row;
pc(i,j+N(cb)-1)=col;
end
end
end



