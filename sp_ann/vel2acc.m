%% *Velocigram Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Politecnico di Milano - DICA
% Copyright 2016_
%% NOTES
% _vel2acc_: function to differentiate velocigram to get accelerogram
% by applying butterworth filter and detrending.
%% INPUT:
% * _dt (sampling time step)_
% * _thv (input velocigram)_
% * _lfr (corner frequency)_
% * _hfr (cut-off frequency)_
%% OUTPUT:
% * _tha (band-pass filtered acceleration time-history column
% vector)_
% * _thv (velocity time-history column vector)_
% * _thd (displacement time-history column vector)_
function [varargout] = vel2acc(varargin)
    %% SET-UP
    % _time-step_
    dtm = varargin{1};
    % _accelerogram_
    thv = varargin{2}(:);
    % _default corner frequency (high-pass filter)_
    lfr =.01;
    % _default cutoff frequency (low-pass filter)_
    hfr = 25;
    % _default butterworth order_
    bfo = 4;
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
    %% BUTTERWORTH FILTER
    [bfb,bfa,flag] = create_butter_filter(bfo,lfr,hfr,fNy);
    %% PROCESSING
    if flag
        %% PROCESSING VELOCITY
        % _pad definition_
        ntm = numel(thv);
        npd = ceil(40/dtm);
        ntm_pad = ntm + 2*npd;
        thv_pad = zeros(ntm_pad,1);
        % _padding_
        thv_pad(:) = padarray(thv,npd,'both');
        % _base-line correction_
        thv_pad = detrend(thv_pad);
        % _acasual filtering_
        thv = filtfilt(bfb,bfa,thv_pad);
        %% COMPUTING/PROCESSING DISPLACEMENT
        % _time integration_
        thd_pad = cumtrapz(thv)*dtm;
        % _base-line correction_
        thd_pad = detrend(thd_pad);
        % _applying cosinus taper_
        thd_pad = cos_taper(thd_pad);
        % _acasual filtering_
        thd = filtfilt(bfb,bfa,thd_pad);
    else
        % _time integration_
        thd = cumtrapz(thv)*dtm;
    end
    %% BACK TO ACCELERATION
    % _time differentiation_
    thv = [0;diff(thd)/dtm];
    tha = [0;diff(thv)/dtm];
    % _time integration_
    thv = cumtrapz(tha)*dtm;
    thd = cumtrapz(thv)*dtm;
    %% OUTPUT
    if flag
        varargout{1} = tha(npd+1:ntm+npd,1);
        varargout{2} = thv(npd+1:ntm+npd,1);
        varargout{3} = thd(npd+1:ntm+npd,1);
    else
        varargout{1} = tha(:);
        varargout{2} = thv(:);
        varargout{3} = thd(:);
    end
    return
end