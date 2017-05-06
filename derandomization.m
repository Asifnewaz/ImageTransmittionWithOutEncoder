function [f] = derandomization(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[M,N]=size(x) ;
g=im2col(x, [M,N], [M,N], 'distinct');
h=dec2bin(double(g));
[M1,N1]=size(h);
z=zeros (M1,N1) ;
for i=1:M1
for jj=1:N1
z(i,jj)=str2num(h(i,jj));
end;
end;
h_rec = zeros(size(z));
for j=1:8
HH=z(:,j);
FF=reshape(HH,M,N);
n = [10,5,12,5,10,8,14,10,5,12,5,10,8,14,10,5,12,5,10,8,14,10,5,12,5,10,8,14];
%n=[4,1,5,6];
%n = [10,5,12,5,10,8,14,10,5,12,5,10,8,14];
% n=[5 2 1 4 3 1];
%n=[4 2 4 4 2];
[pr,pc] = chaomat(n);
pim = chaoperm(FF,pr,pc,3,'backward');
h_rec(:,j)=reshape(pim,M*N,1);
end
for i=1:M1
for jj=1:N1
zn(i,jj)=num2str(h_rec(i,jj));
end;
end;
hn=bin2dec(zn);
f=col2im(hn, [M,N], [M,N], 'distinct');
end
