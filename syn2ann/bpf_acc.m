%% *Acceleration Processing*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2014-15_
%% NOTES
% _bpf_acc_: function to detrend, filter (band pass) and integrate in
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
function [varargout] = bpf_acc(varargin)
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
        % _pad definition_
        %
        % number of time steps of original record
        ntm     = numel(tha);
        % number of padding points (Boore&Bommer,2005)
        npd     = ceil(1.5*bfo/min([lfr;hfr]));
        % number of time-steps of padded record
        ntm_pad = ntm + 2*npd;
        % padded acceleration
        tha_pad = zeros(ntm_pad,1);
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
        tha_pad(1:ntm_pad,1) = padarray(tha,npd,'both');
        %
        % _acceleration acausal Butterworth filtering_
        %
        tha_pad = filtfilt(bfb,bfa,tha_pad);
        
        %% *TIME INTEGRATION*
        [tha_pad,thv_pad,thd_pad] = ...
            integr_diff_avd(dtm,tha_pad,bfb,bfa);
        %% *OUTPUT*
        varargout{1} = tha_pad(npd+1:ntm+npd,1);
        varargout{2} = thv_pad(npd+1:ntm+npd,1);
        varargout{3} = thd_pad(npd+1:ntm+npd,1);
    else
        %% *TIME INTEGRATION*
        [tha,thv,thd] = integr_diff_avd(dtm,tha);
        %% *OUTPUT*
        varargout{1} = tha;
        varargout{2} = thv;
        varargout{3} = thd;
    end
    return
end
