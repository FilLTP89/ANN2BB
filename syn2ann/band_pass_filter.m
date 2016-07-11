%% *Acceleration Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _band_pass_filter_: function to detrend, filter (band pass) and integrate in
% time input signal (to get velocity and displacement)
%% INPUT:
% * _dtm (sampling time step)_
% * _tha (input accelerogram)_
% * _lfr (corner frequency)_
% * _hfr (cut-off frequency)_
%% OUTPUT:
% * _tha (band-pass filtered acceleration time-history column vector)_
% * _thv (velocity time-history column vector)_
% * _thd (displacement time-history column vector)_
function [varargout] = band_pass_filter(varargin)
    %% SET-UP
    % _time-step_
    dtm = varargin{1};
    % _accelerogram_
    tha = varargin{2}(:);
    % _default corner frequency (high-pass filter)_
    lfr =.01;
    % _default cutoff frequency (low-pass filter)_
    hfr = 25;
    % _default butterworth order_
    bfo = 2;
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
        %% PROCESSING ACCELERATION
        % _pad definition_
        ntm = numel(tha);
        npd = ceil(40/dtm);
        ntm_pad = ntm + 2*npd;
        tha_pad = zeros(ntm_pad,1);
        % _base-line correction_
        tha = detrend(tha);
        % _applying cosinus taper_
%         tha = cos_taper(tha);
        tha = taper_fun(tha,5,1,1);
        % _padding_
        tha_pad(:) = padarray(tha,npd,'both');
        % _acausal Butterworth filtering_
        tha_pad = filtfilt(bfb,bfa,tha_pad);
        %% TIME INTEGRATION
        [tha_pad,thv_pad,thd_pad] = integr_diff_avd(dtm,tha_pad,bfb,bfa);
        %% OUTPUT
        varargout{1} = tha_pad(npd+1:ntm+npd,1);
        varargout{2} = thv_pad(npd+1:ntm+npd,1);
        varargout{3} = thd_pad(npd+1:ntm+npd,1);
    else
        %% TIME INTEGRATION
        [tha,thv,thd] = integr_diff_avd(dtm,tha);
        %% OUTPUT
        varargout{1} = tha;
        varargout{2} = thv;
        varargout{3} = thd;
    end
    return
end

%% *Acceleration Time Integration*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _integr_diff_avd_: function that integrates and differentiate acceleration
% signal
%% INPUT:
% * dtm (sampling time step)
% * tha (BP filtered input signal)
% * bfb (Butterworth's filter b coefficient)
% * bfa (Butterworth's filter a coefficient)
%% OUTPUT:
% * tha (acceleration time-history vector (after differentiation))
% * thv (velocity time-history vector (after differentiation))
% * thd (displacement time-history vector (after differentiation))

function [varargout] = integr_diff_avd(varargin)
    %% SET-UP
    dtm = varargin{1};
    tha = varargin{2}(:);
    flag=false;
    if nargin>2
        bfb = varargin{3};
        bfa = varargin{4};
        flag = true;
    end
    if flag
        %% COMPUTING/PROCESSING VELOCITY
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        % _base-line correction_
        thv = detrend(thv);
        % _acasual filtering_
        thv = filtfilt(bfb,bfa,thv);
        %% COMPUTING/PROCESSING DISPLACEMENT
        % _time integration_
        thd = cumtrapz(thv)*dtm;
        % _base-line correction_
        thd = detrend(thd);
        % _applying cosinus taper_
%         thd = cos_taper(thd);
        thd = taper_fun(thd,5,1,1);
        % _acasual filtering_
        thd = filtfilt(bfb,bfa,thd);
        %% BACK TO ACCELERATION
        % _time differentiation_
        thv = [diff(thd)/dtm;0];
        tha = [diff(thv)/dtm;0];
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        thd = cumtrapz(thv)*dtm;
    else
        %% COMPUTING/PROCESSING VELOCITY
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        %% COMPUTING/PROCESSING DISPLACEMENT
        % _time integration_
        thd = cumtrapz(thv)*dtm;
    end
    
    %% OUTPUT
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end