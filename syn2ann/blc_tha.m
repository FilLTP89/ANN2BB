function [varargout] = blc_tha(varargin)
    dtm = varargin{1};
    tha = varargin{2}(:);
    tp = 0;
    if nargin>2
        tp = varargin{3};
    end
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    %
    % _acceleration 0th order correction (pre-event)_
    %
    idx0 = vtm>=0&vtm<=tp; 
    idx1 = vtm>=tp;
    blc  = mean(tha(idx0));
    tha = tha-blc;
    %
    % _quadratic fitting of velocity trace_
    %
    thv  = cumtrapz(tha)*dtm;
    thv(idx1(1)) = 0;
    pq  = polyfit(vtm(idx1),thv(idx1),2);
    thv = polyval(pq,vtm);
    %
    % _subtracting derivative to 0th order corrected acceleration_
    %
    tha  = tha-avd_diff(dtm,thv);
    %
    % _4th-order Butterworth causal filtering_
    %
    [bfb,bfa,~] = create_butter_filter(4,0.02,[],0.5/dtm);
    tha = filter(bfb,bfa,tha);
    %
    % _integration_
    %
    [~,thv,thd] = integr_diff_avd(dtm,tha);
    
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end