%%%%%%%%%%%
%This code is for Fig 1 of "Comparison of Asymmetrically Clipped Optical OFDM and DC-Biased Optical OFDM in AWGN"
%IEEE COMMUNICATIONS LETTERS, VOL. 12, NO. 5, MAY 2008
%%%%%%%%%%%
clc;
close all
%hold on
clear all
maxSNR = 35;
minSNR = 0;
acrSNR = 2;
lenghtOfIFFT=256;%
numOfSymbolInAnOFDMSymbol = lenghtOfIFFT/4;%For ACO-OFDM, it is lenghtOfIFFT/4 due to the Hermitian
%structure and only the odd subcarriers are used for transmission, it
%is equal to lenghtOfIFFT, otherwise.
bitPerSymbol= 10;
numOfOFDMSymbol=5000;
bitLength=  numOfSymbolInAnOFDMSymbol * numOfOFDMSymbol * bitPerSymbol;
bitSource = round(rand(1,bitLength));
bit2Decimal=bi2de(reshape(bitSource',bitPerSymbol,bitLength/bitPerSymbol).','left-msb');
decimalForModulation=reshape(bit2Decimal,numOfSymbolInAnOFDMSymbol,numOfOFDMSymbol);
symbol = qammod(decimalForModulation, 2^bitPerSymbol,'UnitAveragePower', true);
% Hermitian
N=lenghtOfIFFT/2;
symIFFT = zeros(lenghtOfIFFT,numOfOFDMSymbol);
for i=1:numOfSymbolInAnOFDMSymbol
    symIFFT(2*i,:)=symbol(i,:);
end
symIFFT(N+1,:) = symIFFT(N+1,:);
symIFFT(N+2:N*2,:) = flip(conj(symIFFT(2:N,:)),1);%  Hermitian
IFFTOutput = ifft(symIFFT);% Number of IFFT
max(max(IFFTOutput));
size(IFFTOutput);
IFFTOutput((IFFTOutput<0))=0;% ACO
si=1;
s1=zeros(size(minSNR:acrSNR:maxSNR));
fftSymbolForDemod = zeros(numOfSymbolInAnOFDMSymbol, numOfOFDMSymbol);
for snr = minSNR:acrSNR:maxSNR
   I = sum( sum(IFFTOutput.^2) ) / (size(IFFTOutput,1) * size(IFFTOutput,2))/bitPerSymbol;
   %I = sum( sum(IFFTOutput.^2) )/ bitLength;
   Iref = 1/(4*sqrt(pi));
   Eele = I * bitPerSymbol;
   Eopt = sum( sum(IFFTOutput) ) / (size(IFFTOutput,1) * size(IFFTOutput,2));
  % I1 = sum( sum(IFFTOutput) ) / (size(IFFTOutput,1) * size(IFFTOutput,2));
    var  = 10^(-0.1*snr) * I;
    sigma = sqrt(2*var);
    %receiveSymbol=awgn((IFFTOutput),snr,'measured');
    receiveSymbol = sigma * randn(size(IFFTOutput)) + IFFTOutput;
    %receiveSymbol=IFFTOutput;
    fftSymbol=fft(receiveSymbol);
    for i=1:numOfSymbolInAnOFDMSymbol
        fftSymbolForDemod(i,:)=2*fftSymbol(2*i,:);
  
    end
    symbolHat = qamdemod(fftSymbolForDemod,2^bitPerSymbol,'UnitAveragePower', true);
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
axis([0,35,10e-5,1])