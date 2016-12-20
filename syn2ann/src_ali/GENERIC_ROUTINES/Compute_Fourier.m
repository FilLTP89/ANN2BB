function [FT,FAS,f]= Compute_Fourier(signal,dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to compute the Fourier Transform via FFT 
% INPUT:
% signal = time history for which the Fourier Transform is computed
% dt = time interval of input time history
% OUTPUT:
% FT = Fourier Transform (complex-valued) up to the Nyquist frequency
% FAS = Absolute value of Fourier Transform = Fourier Amplitude Spectrum 
%       (real-valued) up to the Nyquist frequency
% f = frequency vector up to the Nyquist frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npun = length(signal);
Nfft = 2^nextpow2(npun); % no. of points (power of 2) on which FFT is performed
df = 1/Nfft/dt; % frequency interval 
%durata = (Nfft-1).*dt; %Ali
%df = 1/durata; %ali
fNy = 1/2/dt; % Nyquist frequency 
f = [0:df:fNy];

FT = dt*fft(signal,Nfft); % application of FFT 
FT = FT(1:length(f)); % Fourier Transform
FAS = abs(FT); % Fourier Amplitude Spectrum 
return
