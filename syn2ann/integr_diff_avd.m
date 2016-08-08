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
        disp('TIME INTEGRATION--->FILTER')
        ntm = numel(tha);
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
        thd = cos_taper(thd);
        % _acasual filtering_
        thd = filtfilt(bfb,bfa,thd);
        %% BACK TO ACCELERATION
        % _time differentiation_
        thv(2:ntm-1,1) = (thd(3:ntm,1)-thd(1:ntm-2,1))./(2*dtm);
        thv(1,1) = 0.0;
        thv(ntm,1) = thv(ntm-1,1);
        tha(2:ntm-1,1) = (thv(3:ntm,1)-thv(1:ntm-2,1))./(2*dtm);
        tha(1,1) = 0.0;
        tha(ntm,1) = tha(ntm-1,1);
        % _time integration_
        thv = cumtrapz(tha)*dtm;
        thd = cumtrapz(thv)*dtm;
    else
        disp('TIME INTEGRATION--->NO FILTER')
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
