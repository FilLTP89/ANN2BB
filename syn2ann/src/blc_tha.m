function [varargout] = blc_tha(varargin)
    %% *SETUP*
    dtm = varargin{1};
    tha = varargin{2}(:);
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    if nargin>2
        tp = varargin{3};
    else
        tp=0;
        try
            tp = PphasePicker(tha,dtm,'sm','n')-1;
        end
    end
    tp = max([tp,0]);
    %
    % _acceleration 0th order correction (pre-event)_
    %
    idx0 = vtm>=0&vtm<=tp;
    idx1 = vtm>=tp;
    idx2 = find(idx1>=1,1,'first');
    blc  = mean(tha(idx0));
    tha = tha-blc;
    %
    % _quadratic fitting of velocity trace_
    %
    thv  = cumtrapz(tha)*dtm;
    % 'C' is the Vandermonde matrix for 'x'
    n = 2; % Degree of polynomial to fit
    C(:,n+1) = ones(numel(vtm),1,class(vtm));
    for i_ = n:-1:1
        C(:,i_) = vtm.*C(:,i_+1);
    end
    %%
    % We use linear equality constraints to force the curve to hit the required point. In
    % this case, 'Aeq' is the Vandermoonde matrix for 'x0'
    Aeq = vtm(idx2).^(n:-1:0);
    % and 'beq' is the value the curve should take at that point
    beq = 0;
    options = optimoptions('lsqlin','Algorithm','active-set');
    pq = lsqlin(C,thv,[],[],Aeq,beq,[],[],[],options);
    thv = polyval(pq,vtm);
    %
    % _subtracting derivative to 0th order corrected acceleration_
    %
    tha  = tha-avd_diff(dtm,thv);
    %
    % _integration_
    %
    [~,thv,thd] = idc_tha(dtm,tha);
    
    
    %% *OUTPUT*
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end
