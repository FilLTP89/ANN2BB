function [varargout] = band_pass_filter(varargin)
    %===============
    % Signal Processing
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % Politecnico di Milano - DICA
    % Copyright 2014-15
    % NOTES
    % band_pass_filter: function to detrend, filter (band pass) and integrate in
    % time input signal (to get velocity and displacement)
    % INPUT:  dtm (sampling time step)
    %         tha (input accelerogram)
    %         lfr (corner frequency)
    %         hfr (cut-off frequency)
    % OUTPUT: tha (band-pass filtered acceleration time-history column vector)
    %         thv (velocity time-history column vector)
    %         thd (displacement time-history column vector)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    dtm=varargin{1};              % time step
    tha=varargin{2}(:);           % accelerogram
    lfr  =.01;                    % default corner frequency
    hfr = 25;                     % default cutoff frequency
    bfo = 2;                      % Butterworth's filter order
    if nargin>=3
        lfr=varargin{3};          % customized corner frequency
    end
    if nargin>=4
        hfr=varargin{4};          % customized corner frequency
    end
    fNy = 0.5/dtm;                % Nyquist frequency
    if hfr>fNy
        hfr = fNy;
    end
    %======================================================================
    %======================================================================
    % BUTTERWORTH FILTER DEFINITION
    %======================================================================
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
    %======================================================================
    %======================================================================
    % INTEGRATION/DIFFERENTIATION/DETRENDING/FILTERING/TAPERING
    %======================================================================
    tha = detrend(tha);                              % detrending average value
    tha = cos_taper(tha);                            % adding cosinus taper
    tha = filtfilt(bfb,bfa,tha);                       % acausal Butterworth filtering
    [tha,thv,thd] = integr_diff_avd(dtm,tha,bfb,bfa);  % integrating
    %======================================================================
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end

function [out_tha,out_thv,out_thd]=integr_diff_avd(varargin)
    %===============
    % NOTES
    % integr_diff_avd: function that integrates and differentiate acceleration
    % signal
    % INPUT:  dtm (sampling time step)
    %         in_acc (BP filtered input signal)
    %         bfb (Butterworth's filter b coefficient)
    %         bfa (Butterworth's filter a coefficient)
    % OUTPUT: out_tha (acceleration time-history vector (after differentiation))
    %         out_thv (velocity time-history vector (after differentiation))
    %         out_thd (displacement time-history vector (after differentiation))
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    dtm = varargin{1};
    in_acc = varargin{2}(:);
    bfb = varargin{3};
    bfa = varargin{4};
    %======================================================================
    %----------------------------------------------------------------------
    % velocity
    %----------------------------------------------------------------------
    out_thv = cumtrapz(in_acc)*dtm;           % integration
    out_thv = filtfilt(bfb,bfa,out_thv);      % filtering
    out_thv = detrend(out_thv);               % detrending
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % displacement
    %----------------------------------------------------------------------
    out_thd = cumtrapz(out_thv)*dtm;          % integration
    out_thd = filtfilt(bfb,bfa,out_thd);      % filtering
    out_thd = detrend(out_thd);               % detrending
    out_thd = cos_taper(out_thd);             % tapering
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % differentiating
    %----------------------------------------------------------------------
    out_thv = [0;diff(out_thd)/dtm];
    out_tha = [0;diff(out_thv)/dtm];
    %======================================================================
    out_thv = cumtrapz(out_tha)*dtm;
    out_thd = cumtrapz(out_thv)*dtm;
    %======================================================================
    return
end


