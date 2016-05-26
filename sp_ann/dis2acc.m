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
function [varargout] = band_pass_filter(varargin)
    
    %% SET-UP
    % _time-step_
    dtm = varargin{1};
    %%
    % _accelerogram_
    thd = varargin{2}(:);
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
    
    %% PROCESSING DISPLACEMENT
    % _acasual filtering_
    thd=filtfilt(bfb,bfa,thd);
    %%
    % _detrending _
    thd=detrend(thd);
    %%
    % _applying cosinus taper_
    thd=cos_taper(thd);
    
    %% BACK TO ACCELERATION
    thv=[0;diff(thd)/dtm];
    tha=[0;diff(thv)/dtm];
    %% OUTPUT
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end