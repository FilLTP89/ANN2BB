%% *Displacement Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _bpf_thd_speed_: function to differentiate displacement record to get velocigram
% and accelerogram by applying butterworth filter and detrending.
%% INPUT:
% * _dtm (sampling time step)_
% * _thd (input accelerogram)_
% * _lfr (corner frequency)_
% * _hfr (cut-off frequency)_
%% OUTPUT:
% * _tha (band-pass filtered acceleration time-history column vector)_
% * _thv (velocity time-history column vector)_
% * _thd (displacement time-history column vector)_
function [varargout] = bpf_thd_speed(varargin)
    %% *SETUP*
    % _time-step_
    dtm = varargin{1};
    % _accelerogram_
    thd = varargin{2}(:);
    % _default corner frequency (high-pass filter)_
    lfr =.15;
    % _default cutoff frequency (low-pass filter)_
    hfr = 25;
    % _default butterworth order_
    bfo = 3;
    % custom corner frequency (high-pass filter)_
    if nargin>=3
        lfr=varargin{3};
    end
    % custom cutoff frequency (low-pass filter)_
    if nargin>=4
        hfr=varargin{4};
    end
    % Nyquist frequency
    fNy = 0.5/dtm;
    if ~isempty(hfr)
        if hfr>fNy
            hfr = fNy;
        end
    end
    
    %% *BUTTERWORTH FILTER*
    %
    
%     [bfb,bfa,~] = create_butter_filter(bfo,lfr,hfr,fNy);
    %% *PROCESSING*
%     thd = filtfilt(bfb,bfa,thd);
    
    %% *BACK TO ACCELERATION*
    thv = avd_diff(dtm,thd);
    % #checked with ali/maria
    thv = tukeywin(numel(thv),0.10).*thv;
    tha = avd_diff(dtm,thv);
    
    %% *OUTPUT*
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end
