%% *Displacement Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _dis2acc_: function to differentiate displacement record to get velocigram
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
function [varargout] = dis2acc(varargin)
    %% SET-UP
    % _time-step_
    dtm = varargin{1};
    % _accelerogram_
    thd = varargin{2}(:);
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
        %% PROCESSING DISPLACEMENT
        % _pad definition_
        ntm = numel(thd);
        npd = ceil(40/dtm);
        ntm_pad = ntm + 2*npd;
        thd_pad = zeros(ntm_pad,1);
        % _padding_
        thd_pad(:) = padarray(thd,npd,'both');
        % _base-line correction_
        thd_pad = detrend(thd_pad);
        % _applying cosinus taper_
        thd_pad = cos_taper(thd_pad);
        % _acasual filtering_
        thd = filtfilt(bfb,bfa,thd_pad);
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