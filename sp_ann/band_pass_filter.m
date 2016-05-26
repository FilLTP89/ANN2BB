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
    %%
    % _accelerogram_
    tha = varargin{2}(:);
    %%
    % _default corner frequency (high-pass filter)_
    lfr =.01;
    %%
    % _default cutoff frequency (low-pass filter)_
    hfr = 25;
    %%
    % _default butterworth order_
    bfo = 2;
    %%
    % custom corner frequency (high-pass filter)_
    if nargin>=3
        lfr=varargin{3};
    end
    %%
    % custom cutoff frequency (low-pass filter)_
    if nargin>=4
        hfr=varargin{4};
    end
    %%
    % Nyquist frequency
    fNy = 0.5/dtm;
    if hfr>fNy
        hfr = fNy;
    end
    
    %% BUTTERWORTH FILTER
    
    if nargin>3 || nargin<3
        % BP filter definition
        [bfb,bfa] = butter(bfo,[lfr hfr]./fNy,'bandpass');
%         fprintf('\nBP FILTER: f(corner): %.2f Hz - f(cut-off): %.2f Hz\n',...
%             lfr,hfr);
    else
        % LF filter definition
        [bfb,bfa] = butter(bfo,lfr./fNy,'high');
%         fprintf('\nLF FILTER: f(corner): %.2f Hz\n',lfr);
    end
    
    %% PROCESSING ACCELERATION
    % _detrending average value_
    tha = detrend(tha);
    % _applying cosinus taper_
    tha = cos_taper(tha);
    % _acausal Butterworth filtering_
    tha = filtfilt(bfb,bfa,tha);
    
    %% TIME INTEGRATION
    [tha,thv,thd] = integr_diff_avd(dtm,tha,bfb,bfa);
    
    %% OUTPUT
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
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
    bfb = varargin{3};
    bfa = varargin{4};
    
    %% COMPUTING/PROCESSING VELOCITY
    % _time integration_
    thv = cumtrapz(tha)*dtm;
    %%
    % _acasual filtering_
    thv = filtfilt(bfb,bfa,thv);
    %%
    % _detrending_
    thv = detrend(thv);
    
    %% COMPUTING/PROCESSING DISPLACEMENT
    % _integration_
    thd = cumtrapz(thv)*dtm;
    %%
    % _acasual filtering_
    thd = filtfilt(bfb,bfa,thd);
    %%
    % _detrending_
    thd = detrend(thd);
    %%
    % _applying cosinus taper_
    thd = cos_taper(thd);
    
    %% BACK TO ACCELERATION
    thv = [0;diff(thd)/dtm];
    tha = [0;diff(thv)/dtm];
    
    thv = cumtrapz(tha)*dtm;
    thd = cumtrapz(thv)*dtm;
    %% OUTPUT
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end


