function [Tt] = jake(v,N)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% v=120;
% N=6;
fc=4096000;
c=3*10^8;
v=v/3.6;
f=2*10^9; % carrier frequency
Tc=1/fc;%chip duration
nu=1;%1 samples per chip
fs=nu*fc; % sampling frequency
Ts=1/fs;
wM=2*pi*f*v/c;
%%%%%%%%%%%%%%%%%%
tet=rand(1,200)*2*pi;
t=0:128*Tc:5000*128*Tc;
No=16;
%%%%%%%%%%%%%%%%%%
A=hadamard(No);
%%%%%%%%%%%%%%%%%%
for k=1:N
Ttg=zeros(1,length(t));
tet=rand(1,200)*2*pi;
for l=1:No
B_n=pi*l/(No);
An=2*pi*(l-.5)/(4*No);
wn=wM*cos(An);
Ttg=A(k,l)*sqrt(2/No)*(cos(B_n)+j*sin(B_n)).*cos(wn*t+ tet(l))+Ttg;
end
Tt(k,:)=Ttg;

end
Tt=Tt(:,1:5000);
for z1=1:N
Tt(z1,:)=(Tt(z1,:)-mean(Tt(z1,:)))/std(Tt(z1,:));
end


