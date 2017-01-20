% freq_integrate.m
% performs integration in the frequency domain by applying an fft, dividing
% by 2*pi*f, and then ifft back to time domain
% usage:
% out = freq_integrate(data,samples_per_second)
% designed to integrate from acceleration to velocity or velocity to
% displacement
%
% see companion freq_differentiate to differentiate from displacement to
% velocity or velocity to acceleration

% advantage of this method is that frequency domain is more stable than
% time domain and thus you don't get large offsets. It may be suggested to
% do a high pass filter after just to remove long period offset due to the
% lack of integration constant

function out = freq_integrate(data,sps)
    b=zeros(1,length(data));
    a=fft(data);
    for j=1:length(data)
        omega = 2*pi*(sps/2/length(data)) * j;
        b(j) = imag(a(j)) / omega + real(a(j)) / (omega * i);
    end
    out = ifft(b,'symmetric');
    return