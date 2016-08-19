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
    %% *SET-UP*
    % time-step
    dtm = varargin{1};
    % accelerogram
    tha = varargin{2}(:);
    flag=false;
    if nargin>2
        bfb = varargin{3};
        bfa = varargin{4};
        flag = true;
    end

    if flag
        %% *CORRECTED ACCELERATION TIME INTEGRATION*
        disp('--->CORRECTED TIME INTEGRATION')
        % number of time steps of original record
        ntm = numel(tha);
        % velocity
        thv = cumtrapz(tha)*dtm;
        %
        % _velocity base-line correction_
        %
        thv = detrend(thv);
        %
        % _velocity acausal Butterworth filtering_
        %
        thv = filtfilt(bfb,bfa,thv);
        %% *DISPLACEMENT TIME INTEGRATION*
        % displacement
        thd = cumtrapz(thv)*dtm;
        %
        % _displacement base-line correction_
        %
        thd = detrend(thd);
        %
        % _displacement cosinus tapering_
        %
        thd = cos_taper(thd);
        % EQUIVALENT: thd = taper_fun(thd,2.5,1,1);
        %
        % _displacement acasual filtering_
        %
        thd = filtfilt(bfb,bfa,thd);
        
        %% *BACK TO ACCELERATION*
        %
        % _time differentiation (central differences--->E=o(dtm^2))
        % http://oregonstate.edu/instruct/ch490/lessons/lesson11.htm_
        %
        % velocity 
        thv(2:ntm-1,1) = (thd(3:ntm,1)-thd(1:ntm-2,1))./(2*dtm);
        thv([1,ntm],1) = [0.0;thv(ntm-1,1)];
        % acceleration
        tha(2:ntm-1,1) = (thv(3:ntm,1)-thv(1:ntm-2,1))./(2*dtm);
        tha([1,ntm],1) = [0.0;tha(ntm-1,1)];
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        thd = cumtrapz(thv)*dtm;
    else
        %% *ACCELERATION TIME INTEGRATION*
        disp('--->TIME INTEGRATION')
        % velocity
        thv = cumtrapz(tha)*dtm;
        % displacement
        thd = cumtrapz(thv)*dtm;
    end
    
    %% *OUTPUT*
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end
