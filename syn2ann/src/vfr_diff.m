% freq_differentiate.m
% performs differentiation in the frequency domain by applying an fft,
% multiplying by 2*pi*f, and then ifft back to time domain
%
% usage:
% out = freq_differentiate(ths,samples_per_second)
% designed to differentiate from displacement to velocity or velocity to
% acceleration
%
% see companion freq_integrate to integrate from acceleration to
% velocity or velocity to displacement

% advantage of this method is that frequency domain is more stable than
% time domain and thus you don't get large offsets. It may be suggested to
% do fss high pass filter after just to remove long period offset due to the
% lack of integration constant

function [varargout] = vfr_diff(varargin)
    ths = varargin{1}(:);
    sps = varargin{2};
    ntm = numel(ths);
    
    
    nfr = 2^nextpow2(ntm);
    fss = fft(ths,nfr);
    fss = fftshift(fss);
    dfr = sps/nfr;
    vfr = dfr*(-nfr/2:nfr/2-1)';
    omg = 2*pi*vfr;
    keyboard
    thd = sqrt(-1)*omg.*fss;
    thd = real(ifft(ifftshift(thd)));
    
    varargout{1} = thd(1:ntm,1);
    return
end
