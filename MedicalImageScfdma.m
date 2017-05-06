%*****************************************************
%Image Transmission over SC-FDMA
%*****************************************************
clear all
tic % tic starts a stopwatch timer to measure performance
%======= Choose simulation Parameters
SP.FFTsize = 256;
SP.inputBlockSize = 32;
SP.CPsize = 20;
%SP.subband = 15;
SP.subband = 0;  
SP.SNR = 10:1:15;

for sss=1: length(SP.SNR)
%======= Choose channel type

%======= SUI3 Channel
SP.paths= fadchan(SP);

%=====%======= Choose Equalization Type

% SP.equalizerType ='ZERO';
SP.equalizerType ='MMSE';
numSymbols = SP.FFTsize;
Q = numSymbols/SP.inputBlockSize;

%%%%%%%%%========== Channel Generation============

%%%%%%%%%============ SUI3 channel============
%SUI3_1=[SP.paths(1,1) 0 0 SP.paths(2,1) 0 SP.paths(3,1)];
%SP.channel = SUI3_1/norm(SUI3_1);
SP.channel = 1;
%%%%%%%%%========================================
H_channel = fft(SP.channel,SP.FFTsize);

im=imread('cameraman.tif');%image reading

xx=randomization(im);%image randomization

%******************Data Generation********************
f=zeros(256,256);
f=xx;
[M,N]=size(f);
g=im2col(f, [M,N],'distinct');%image to column converter
h=dec2bin(double(g));%pixel value to binary conversion...every value replaced by 8 bits string
[M1,N1]=size(h) ;
z=zeros (M1,N1) ;
clear i j
for i=1:M1
for j=1:N1
z(i,j)= str2num(h(i,j)); %string to number conversion
end;
end;
[M2,N2] = size(z) ;
zz = reshape(z,M2*N2, 1);%parallel data reshaping date to vector


% ********** Dividing the image into blocks***********

nloops = ceil((M2*N2)/SP.inputBlockSize );%number of image %blocks by approximation
new_data = nloops*SP.inputBlockSize ;%new vector proportional to %block size
nzeros = new_data - (M2*N2);%number of zeros to be added to %old data vector
input_data = [zz;zeros(nzeros,1)];%construction of new data %vector
input_data2 = reshape(input_data ,SP.inputBlockSize ,nloops); %reshape the new data to matrix of block size rows to number of %blocks columns save input_data2
%************** transmission ON SC-FDMA ***************
demodata1 = zeros(SP.inputBlockSize ,nloops);% this matrix to %store received data vector
clear jj
for jj =1: nloops % loop for columns
b1= input_data2(:,jj)';%every block size
%%%%%%%%%%%%%%% QPSK Modulation %%%%%%%%%%%%%%%%
tmp = b1;
tmp = tmp*2 - 1;
inputSymbols = (tmp(1,:) + i*tmp(1,:))/sqrt(2);
%%%%%%%%%%%% SC-FDMA Modulation %%%%%%%%%%%%%
inputSymbols_freq = fft(inputSymbols);
inputSamples_ifdma = zeros(1,numSymbols);
inputSamples_lfdma = zeros(1,numSymbols);
%%%%%%%%%%%% Subcarriers Mapping %%%%%%%%%%%%%
inputSamples_ifdma(1+SP.subband:Q:numSymbols) = inputSymbols_freq;
inputSamples_lfdma([1:SP.inputBlockSize]+SP.inputBlockSize*SP.subband) = inputSymbols_freq;
inputSamples_ifdma = ifft(inputSamples_ifdma);
inputSamples_lfdma = ifft(inputSamples_lfdma);

%%%%%%%%%%%%% Add Cyclic Prefix %%%%%%%%%%%%%
TxSamples_ifdma = [inputSamples_ifdma(numSymbols-SP.CPsize+1:numSymbols) inputSamples_ifdma];
TxSamples_lfdma = [inputSamples_lfdma(numSymbols-SP.CPsize+1:numSymbols) inputSamples_lfdma];
%%%%%%%%%%%%% Wireless channel %%%%%%%%%%%%%%
RxSamples_ifdma = filter(SP.channel, 1,TxSamples_ifdma); % Multipath Channel
RxSamples_lfdma = filter(SP.channel, 1,TxSamples_lfdma); % Multipath Channel
%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%
tmp = randn(2, numSymbols+SP.CPsize);
complexNoise = (tmp(1,:) + i*tmp(2,:))/sqrt(2);
noisePower = 10^(-SP.SNR(sss)/10);
%%%%%%%%%%%%%%% Received signal%%%%%%%%%%%%%%
RxSamples_ifdma = RxSamples_ifdma + sqrt(noisePower/Q)*complexNoise;
RxSamples_lfdma = RxSamples_lfdma + sqrt(noisePower/Q)*complexNoise;
%%%%%%%%%%%% Remove Cyclic Prefix%%%%%%%%%%%%
RxSamples_ifdma = RxSamples_ifdma(SP.CPsize+1:numSymbols+SP.CPsize);
RxSamples_lfdma = RxSamples_lfdma(SP.CPsize+1:numSymbols+SP.CPsize);
%%%%%%%%%%%%% SC-FDMA demodulation%%%%%%%%%%%
Y_ifdma = fft(RxSamples_ifdma, SP.FFTsize);
Y_lfdma = fft(RxSamples_lfdma, SP.FFTsize);
%%%%%%%%%%%% subcarriers demapping%%%%%%%%%%%
Y_ifdma = Y_ifdma(1+SP.subband:Q:numSymbols);
Y_lfdma = Y_lfdma([1:SP.inputBlockSize]+SP.inputBlockSize*SP.subband);
%%%%%%%%%%%%%%%%% Equalization %%%%%%%%%%%%%%%%%
H_eff = H_channel(1+SP.subband:Q:numSymbols);
if SP.equalizerType == 'ZERO'
Y_ifdma = Y_ifdma./H_eff;
elseif SP.equalizerType == 'MMSE'
C = conj(H_eff)./(conj(H_eff).*H_eff + 10^(-SP.SNR(sss)/10));
Y_ifdma = Y_ifdma.*C;
end
H_eff = H_channel([1:SP.inputBlockSize]+SP.inputBlockSize*SP.subband);
if SP.equalizerType == 'ZERO'
Y_lfdma = Y_lfdma./H_eff;
elseif SP.equalizerType == 'MMSE'
C = conj(H_eff)./(conj(H_eff).*H_eff + 10^(-SP.SNR(sss)/10));
Y_lfdma = Y_lfdma.*C;
end
EstSymbols_ifdma = ifft(Y_ifdma);
EstSymbols_lfdma = ifft(Y_lfdma);
%%%%%%%%%%%%%%%% demodulation%%%%%%%%%%%%%%%%
EstSymbols_ifdma = sign(real(EstSymbols_ifdma));
EstSymbols_ifdma =(EstSymbols_ifdma+1)/2;
EstSymbols_lfdma = sign(real(EstSymbols_lfdma));
EstSymbols_lfdma = (EstSymbols_lfdma+1)/2;
demodata1_ifdma(:,jj) = EstSymbols_ifdma(:); % the output of scfdma columns%storing of received image data
demodata1_lfdma(:,jj) = EstSymbols_lfdma(:);
% the output of scfdma columns%storing of received image data
end
    
%***************** Received image ***************
[M3,N3] = size(demodata1_ifdma);
% demodata2 = demodata1(:);
yy1_ifdma = reshape (demodata1_ifdma,M3,N3);
%reshaping the matrix to vector
yy1_lfdma = reshape (demodata1_lfdma,M3,N3);
%reshaping the matrix to vector
received_image_ifdma = yy1_ifdma(1:M2*N2);%taking the original data
received_image_lfdma = yy1_lfdma(1:M2*N2);%taking the original data
%************** Regeneration of image **************
zz1_ifdma=reshape(received_image_ifdma,M2* N2,1);
%reshaping to M2*N2 vector
zz1_lfdma=reshape(received_image_lfdma,M2* N2,1);
%reshaping to M2*N2 vector
yy_ifdma = reshape(zz1_ifdma,M2, N2);
yy_lfdma = reshape( zz1_lfdma,M2, N2);
clear i j
for i=1:M1
for j=1:N1
zn_ifdma(i,j)=num2str(yy_ifdma(i,j));
zn_lfdma(i,j)=num2str(yy_lfdma(i,j));
end;
end;
hn_ifdma=bin2dec(zn_ifdma);
hn_lfdma=bin2dec(zn_lfdma);
gn_ifdma=col2im(hn_ifdma, [M,N], [M,N], 'distinct');
gn_lfdma=col2im(hn_lfdma, [M,N], [M,N], 'distinct');
y1_ifdma=derandomization(gn_ifdma);
y1_lfdma=derandomization(gn_lfdma);
y1_ifdma=y1_ifdma/255;% Why this is divide unknown???
y1_lfdma= y1_lfdma/255;
% ***************** The output results****************
figure (1)
imshow(im)
figure (2)
imshow(y1_ifdma)
figure (3)
imshow(y1_lfdma)

MSE1_ifdma=sum(sum((double(im) - y1_ifdma).^2))/(M*N);
PSNR_ifdma(sss)=10*log10((256*256)/MSE1_ifdma);

MSE1_lfdma=sum(sum((double(im) - y1_lfdma).^2))/(M*N);
PSNR_lfdma(sss)=10*log10((256*256)/MSE1_lfdma);
% MSE1_ifdma=sum(sum((double(im)/255-y1_ifdma).^2))/prod(size(im));
% %PSNR_ifdma(sss)=10*log(1/MSE1_ifdma)/log(10);
% PSNR_ifdma=10*log(1/MSE1_ifdma)/log(10);
% % figure (3)
% % imshow(y1_lfdma)
% MSE1_lfdma=sum(sum((double(im)/255-y1_lfdma).^2))/prod(size(im));
% %PSNR_lfdma(sss)=10*log(1/MSE1_lfdma)/log(10);
% PSNR_lfdma=10*log(1/MSE1_lfdma)/log(10);
%y1_ifdma_DFT(:,:,sss)=y1_ifdma;
%y1_lfdma_DFT(:,:,sss)=y1_lfdma;
end

save PSNR_ifdma; 
save PSNR_lfdma;


%save y1_ifdma_DFT; 
%save y1_lfdma_DFT;

%%%%%%%%% Plot the Results %%
figure(44)
plot(SP.SNR,PSNR_ifdma,'rx-',SP.SNR,PSNR_lfdma,'mx-');
legend('DFT-IFDMA','DFT-LFDMA')
xlabel('SNR (dB)'); ylabel('PSNR(dB)');
axis([0 20 0 60])
grid on
toc
%PSNR(e)=10*log(1/MSE1)/log(10);% Peak signal-to-noise ratio 10*log10((max possible pixel value of the image)^2/MSE)
%end

%****************** End of file ********************
%*****************************************************