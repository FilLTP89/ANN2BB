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
    %%
    % _accelerogram_
    thv = varargin{2}(:);
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
        fprintf('\nBP FILTER: f(corner): %.2f Hz - f(cut-off): %.2f Hz\n',...
            lfr,hfr);
    else
        % LF filter definition
        [bfb,bfa] = butter(bfo,lfr./fNy,'high');
        fprintf('\nLF FILTER: f(corner): %.2f Hz\n',lfr);
    end
    
    %% PROCESSING VELOCITY
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