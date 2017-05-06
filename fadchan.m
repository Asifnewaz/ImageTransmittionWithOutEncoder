function [paths] = fadchan(SP)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

N=8;
OR = 4;
M = 256;
Dop_res = 0.1;
res_accu = 20;
% P = [ 0 -5 -10 ];
% K = [ 0.1 0 0 ];
% tau = [ 0.0 5 10];
% Dop = [ 2 1.5 2.5];
% Dop = [ 0 0 0];
P = [ 0 -5 -10 ];
K = [ 1 0 0 ];
tau = [ 0.0 0.5 1.0 ];
Dop = [ 0.4 0.4 0.4 ];
ant_corr = 0.4;
Fnorm = -1.5113;
P = 10.^(P/10); % calculate linear power
s2 = P./(K+1); % calculate variance
m2 = P.*(K./(K+1)); % calculate constant power
m = sqrt(m2);
L = length(P); % number of taps
paths_r = sqrt(1/2)*(randn(L,N) + j*randn(L,N)).*((sqrt(s2))' * ones(1,N));
paths_c = m' * ones(1,N);
for p = 1:L
D = Dop(p) / max(Dop) / 2; % normalize to highest Doppler
f0 = [0:M*D]/(M*D); % frequency vector
PSD = 0.785*f0.^4 - 1.72*f0.^2 + 1.0; % PSD approximation


filt = [ PSD(1:end-1) zeros(1,M-2*M*D) PSD(end:-1:2)]; % S(f)
filt = sqrt(filt); % from S(f) to |H(f)|
filt = ifftshift(ifft(filt)); % get impulse response
filt = real(filt); % want a real-valued filter
filt = filt / sqrt(sum(filt.^2)); % normalize filter
path = fftfilt(filt, [ paths_r(p,:) zeros(1,M) ]);
paths_r(p,:) = path(1+M/2:end-M/2);
end;
paths = paths_r + paths_c;
paths = paths * 10^(Fnorm/20); % multiply all coefficients with F
end

