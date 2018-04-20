%/*************************************************************************
%    > Copyright (C) 2017 All rights reserved.
%    > File Name: ACO_OFDM.m
%    > Author: Chauency
%    > Mail: chauencychan@gmail.com
%    > Created Time: 2017-08-20
%   > Description: Matlab code for ACO-OFDM,
% ************************************************************************/
clc;
clear;
maxSNR = 15;
minSNR = 0;
acrSNR = 1;
lenghtOfIFFT=128;%
numOfSymbolInAnOFDMSymbol = lenghtOfIFFT/4;%For ACO-OFDM, it is lenghtOfIFFT/4 due to the Hermitian
%structure and only the odd subcarriers are used for transmission, it
%is equal to lenghtOfIFFT, otherwise.
bitPerSymbol=4;
numOfOFDMSymbol=5000;
bitLength=  numOfSymbolInAnOFDMSymbol * numOfOFDMSymbol * bitPerSymbol;
bitSource = round(rand(1,bitLength));
bit2Decimal=bi2de(reshape(bitSource',bitPerSymbol,bitLength/bitPerSymbol).','left-msb');
decimalForModulation=reshape(bit2Decimal,numOfSymbolInAnOFDMSymbol,numOfOFDMSymbol);
symbol = qammod(decimalForModulation, 2^bitPerSymbol);
% Hermitian
N=lenghtOfIFFT/2;
symIFFT = zeros(lenghtOfIFFT,numOfOFDMSymbol);
for i=1:numOfSymbolInAnOFDMSymbol
    symIFFT(2*i,:)=symbol(i,:);
end
symIFFT(N+1,:) = symIFFT(N+1,:);
symIFFT(N+2:N*2,:) = flip(conj(symIFFT(2:N,:)),1);%  Hermitian
IFFTOutput = ifft(symIFFT);% Number of IFFT is 256
IFFTOutput((IFFTOutput<0))=0;% ACO

si=1;
s1=zeros(1,maxSNR - minSNR + 1);
fftSymbolForDemod = zeros(numOfSymbolInAnOFDMSymbol, numOfOFDMSymbol);
for snr = minSNR:acrSNR:maxSNR
    receiveSymbol=awgn(IFFTOutput,2*snr,'measured');
    fftSymbol=fft(receiveSymbol);
    for i=1:numOfSymbolInAnOFDMSymbol
        fftSymbolForDemod(i,:)=2*fftSymbol(2*i,:);
    end
    symbolHat = qamdemod(fftSymbolForDemod,2^bitPerSymbol);
    bitHat = de2bi(symbolHat,'left-msb'); 
    bitHat=reshape(bitHat.',1,[]);
    [number_of_errors,bit_error_rate]=biterr(bitSource,bitHat);
    s1(si)=bit_error_rate;
    si=si+1;
end
snr=minSNR:acrSNR:maxSNR;
figure(1);
semilogy(snr,s1,'s-');
xlabel('SNR')
ylabel('BER')
grid on
