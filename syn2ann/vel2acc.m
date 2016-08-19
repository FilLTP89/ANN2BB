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
    %% *SET-UP*
    % time-step
    dtm = varargin{1};
    % velocigram
    thv = varargin{2}(:);
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
        %% *PROCESSING VELOCITY*
        disp('--->CORRECTING VELOCITY')
        %
        % _velocity base-line correction_
        %
        thv = detrend(thv);
        %
        % _velocity cosinus tapering_
        %
        thv = cos_taper(thv);
        % EQUIVALENT: thd = taper_fun(thd,2.5,1,1);
        %
        % _pad definition_
        %
        % number of time steps of original record
        ntm     = numel(thv);
        % number of padding points (Boore&Bommer,2005)
        npd     = ceil(1.5*bfo./min([lfr;hfr])./dtm);
        % number of time-steps of padded record
        ntm_pad = ntm + 2*npd;
        %
        % _velocity padding_
        %
        thv_pad = padarray(thv,npd,'both');
        %
        % _velocity acausal Butterworth filtering_
        %
        thv_pad = filtfilt(bfb,bfa,thv_pad);
        
        %% *CORRECTED DISPLACEMENT TIME INTEGRATION*
        % displacement
        thd_pad = cumtrapz(thv_pad)*dtm;
        %
        % _displacement base-line correction_
        %
        thd_pad = detrend(thd_pad);
        %
        % _displacement cosinus tapering_
        %
        thd_pad = cos_taper(thd_pad);
        % EQUIVALENT: thd = taper_fun(thd,2.5,1,1);
        %
        % _displacement acasual filtering_
        %
        thd_pad = filtfilt(bfb,bfa,thd_pad);
        
        %% *BACK TO ACCELERATION*
        %
        % _time differentiation (central differences--->E=o(dtm^2))
        % http://oregonstate.edu/instruct/ch490/lessons/lesson11.htm_
        %
        % velocity 
        thv_pad(2:ntm_pad-1,1) = (thd_pad(3:ntm_pad,1)-thd_pad(1:ntm_pad-2,1))./(2*dtm);
        thv_pad([1,ntm_pad],1) = [0.0;thv_pad(ntm_pad-1,1)];
        % acceleration
        tha = zeros(ntm_pad,1);
        tha(2:ntm_pad-1,1) = (thv_pad(3:ntm_pad,1)-thv_pad(1:ntm_pad-2,1))./(2*dtm);
        tha([1,ntm_pad],1) = [0.0;tha(ntm_pad-1,1)];
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        thd = cumtrapz(thv)*dtm;
    else
        %% *DISPLACEMENT TIME INTEGRATION*
        % displacement
        thd = cumtrapz(thv)*dtm;
        %% *BACK TO ACCELERATION*
        %
        % _time differentiation (central differences--->E=o(dtm^2))
        % http://oregonstate.edu/instruct/ch490/lessons/lesson11.htm_
        %
        % velocity 
        thv(2:ntm-1,1) = (thd(3:ntm,1)-thd(1:ntm-2,1))./(2*dtm);
        thv([1,ntm],1) = [0.0;thv(ntm-1,1)];
        % acceleration
        tha = zeros(ntm,1);
        tha(2:ntm-1,1) = (thv(3:ntm,1)-thv(1:ntm-2,1))./(2*dtm);
        tha([1,ntm],1) = [0.0;tha(ntm-1,1)];
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        thd = cumtrapz(thv)*dtm;
    end
    
    %% *OUTPUT*
    if flag
        varargout{1} = tha;%(npd+1:ntm+npd,1);
        varargout{2} = thv;%(npd+1:ntm+npd,1);
        varargout{3} = thd;%(npd+1:ntm+npd,1);
    else
        varargout{1} = tha(:);
        varargout{2} = thv(:);
        varargout{3} = thd(:);
    end
    return
end
