%% *Acceleration Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _bpf_tha_: function to detrend, filter (band pass) and integrate in
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
%% REFERENCES
% @Article{Boore_Bommer_2005,
%   Title                    = {{Processing of strong-motion accelerograms: needs, options and consequences}},
%   Author                   = {Boore, David M. and Bommer, Julian J.},
%   Journal                  = {Soil Dynamics and Earthquake Engineering},
%   Year                     = {2005},
%
%   Month                    = {February},
%   Number                   = {2},
%   Pages                    = {93--115},
%   Volume                   = {25},
%
%   Doi                      = {10.1016/j.soildyn.2004.10.007},
%   File                     = {Boore_Bommer_2005.pdf:Boore_Bommer_2005.pdf:PDF},
%   ISSN                     = {02677261},
%   Keywords                 = {baseline adjustments,filters,instrument corrections,signal-to-noise ratios,strong-motion accelerograms},
%   Url                      = {http://linkinghub.elsevier.com/retrieve/pii/S0267726104001708}
% }
function [varargout] = bpf_tha(varargin)
    %% *SET-UP*
    % time-step
    dtm = varargin{1};
    % accelerogram
    tha = varargin{2}(:);
    % default corner frequency (high-pass filter)
    lfr =.05;
    % default cutoff frequency (low-pass filter)
    hfr = 25;
    % default butterworth order
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
    
    %% *CREATE BUTTERWORTH FILTER*
    [bfb,bfa,flag] = create_butter_filter(bfo,lfr,hfr,fNy);
    
    if flag
        %% *PROCESSING ACCELERATION*
        disp('--->CORRECTING ACCELERATION')
        %
        %  _pad definition_
        %
        % number of padding points (Boore&Bommer,2005)
        npd     = ceil(1.5*bfo./min([lfr;hfr])./dtm);
        %
        % _acceleration base-line correction_
        %
        tha = detrend(tha);
        %
        % _acceleration cosinus tapering_
        %
        tha = cos_taper(tha);
        % EQUIVALENT: tha = taper_fun(tha,2.5,1,1);
        %
        % _padding acceleration_
        %
        tha = padarray(tha,npd,'both');
        %
        % _acceleration acausal Butterworth filtering_
        %
        tha = filtfilt(bfb,bfa,tha);
    end
    %% *TIME INTEGRATION*
    [tha,thv,thd] = integr_diff_avd(dtm,tha,bfb,bfa);
    %% *OUTPUT*
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    varargout{1} = tha;%(npd+1:ntm+npd,1);
    varargout{2} = thv;%(npd+1:ntm+npd,1);
    varargout{3} = thd;%(npd+1:ntm+npd,1);
    varargout{4} = vtm;
    return
end
