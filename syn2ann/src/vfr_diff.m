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
    ths = [zeros(ntm,1);ths;zeros(ntm,1)];
    ntm = numel(ths);
    nfr = 2^nextpow2(ntm);
    fss = fft(ths,nfr);
    dfr = sps/nfr;
    vfr = dfr*(0:nfr-1)';
    
%     omega = zeros(nfr,1);
    % frequency domain differentation
%     omega(1:nfr/2+1,1)   = 2*pi*vfr(1:nfr/2+1,1);
%     omega(nfr/2+2:nfr,1) = -2*pi*vfr(nfr/2:-1:2,1);
    omega = 2*pi*vfr;
    thd = sqrt(-1)*omega.*fss;
    thd = real(ifft(thd,'symmetric'));
    
    varargout{1} = thd(1:ntm,1);
    return
end
