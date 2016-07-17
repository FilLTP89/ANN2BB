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
% REFERENCES:
% Frequency domain representation returned as a vector, matrix, or 
% multidimensional array. If X is of type single, then fft natively computes 
% in single precision, and Y is also of type single. Otherwise, Y is 
% returned as a type double.
% The size of Y is as follows:
% 
%     For Y = fft(X) or Y = fft(X,[],dim), the size of Y is equal to the size of X.
% 
%     For Y = fft(X,n,dim), the value of size(Y,dim) is equal to n, 
%     while the size of all other dimensions remains as in X.
% 
% If X is real, then Y is conjugate symmetric, and the number of unique 
% points in Y is ceil((n+1)/2).

function [varargout]=super_fft(varargin)
    %% SET-UP
    % _time step_
    dtm = varargin{1};
    % _acceleration_
    thr=varargin{2}(:);
    % _Konno-Ohmachi Smoothing flag_
    kos = varargin{3};
    % _output selection_
    out_sel = 1:6;
    if nargin>3
        out_sel = varargin{4};
    end
    % _frequency vector_
    nfr=2^nextpow2(numel(thr));
    if any(out_sel==1)
        dfr = 1/dtm/(nfr-1);
        vfr = dfr*(0:nfr-1);
        varargout{out_sel==1} = vfr(:);
    end
    %
    % _fft complete response_
    %
    fsr = fft(thr,nfr).*dtm;
    if any(out_sel==4)
        varargout{out_sel==4} = fsr(:);
    end
    %
    % _fft amplitude_
    %
    if any(out_sel==2) || logical(kos)
        amp = abs(fsr);
        % _Konno-Ohmachi smoothing_
        if logical(kos)
            amp=smooth_KO(amp,40);
        end
        varargout{out_sel==2} = amp(:);
    end
    %
    % _phase angle of fft_
    %
    if any(out_sel==3)
        ang = unwrap(angle(fsr));
        varargout{out_sel==3} = ang(:);
    end
    %
    % _power spectral density_
    if any(out_sel==5)
        psd = (fsr).*(conj(fsr))/nfr;
        varargout{out_sel==5} = psd(:);
    end
    %
    % _main period_
    %
    if any(out_sel==6)
        idx0 = find(vfr>=.1,1,'first');
        idx1 = find(vfr<=40,1,'last');
        nTm1 = sum(amp(idx0:idx1).^2./vfr(idx0:idx1));
        nTm2 = sum(amp(idx0:idx1).^2);
        nTm = nTm1/nTm2;
        varargout{out_sel==6} = nTm;
    end
    return
end
