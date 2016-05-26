%% *Fourier tranform*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _super_fft_: function to evaluate:
%  # Fourier Spectrum (Amplitude and Phase),
%  # Power Spectral Density
%  # Main Period of input signal
%% INPUT:
% * _dt (sampling time step)_
% * _thr (input signal)_
% * _K0_smooth (flag for Konno-Ohmachi smoothing (K0_smooth=1) or not (K0_smooth=0))_
% * _nfr (number of points for FFT computation (optional))_
%% OUTPUT:
% * _vfr (frequency column vector)_
% * _amp (FS amplitude column vector)_
% * _ang (FS phase column vector)_
% * _fsr_ (Fourier Complex Spectrum)
% * _psd (Power Spectral Density column vector)_
% * _nTm (Main Natural Period evaluated between fcut_low and fcut_high)_
%% N.B.:
% nfr should be odd number to be able to transpose the
% conjugate part [size(fsr)=nfr (odd)]
function [varargout]=super_fft(varargin)
    
    %% SET-UP
    %%
    % _time step_
    dtm = varargin{1};
    %%
    % _acceleration_
    thr=varargin{2}(:);
    %%
    % _Konno-Ohmachi Smoothing flag_
    kos = varargin{3};
    
    %% FREQUENCY VECTOR_
    % _nfr definition_
    
    if nargin>3
        nfr=varargin{4};
    else
        nfr=2^nextpow2(numel(thr))+1;
    end
    %%
    % _frequency vector_
    dfr=1/dtm/(nfr-1);     
    vfr=(0:nfr-1)*dfr;    
    vfr = vfr(:);
    
    %% FOURIER RESPONSE
    
    %%
    % _fft complete response_
    % fsr=fft(thr,nfr)./numel(thr);
    fsr = fft(thr,nfr).*dtm;
    
    %%
    % _amplitude of fft response_
    amp = abs(fsr);
    %%
    % _Konno-Ohmachi smoothing_
    if kos
        amp=smooth_KO(amp,40);
    end
    amp = amp(:);
    
    %%
    % _phase angle of fft (to be unwrapped)_
    ang = unwrap(angle(fsr));
    
    %% POWER SPECTRAL DENSITY
    psd=(fsr).*(conj(fsr))/nfr;
    
    %% MAIN PERIOD
    idx0 = find(vfr>=.1,1,'first');
    idx1 = find(vfr<=40,1,'last');
    nTm1 = sum(amp(idx0:idx1).^2./vfr(idx0:idx1));
    nTm2 = sum(amp(idx0:idx1).^2);
    nTm = nTm1/nTm2;
    
    %% OUTPUT
    varargout{1}=vfr(:);
    varargout{2}=amp(:);
    varargout{3}=ang(:);
    varargout{4}=fsr(:);
    varargout{5}=psd(:);
    varargout{6}=nTm;
    return
end
