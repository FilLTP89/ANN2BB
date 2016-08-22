function [varargout] = blc_tha(varargin)
    %% *SETUP*
    dtm = varargin{1};
    tha = varargin{2}(:);
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    tp = PphasePicker(tha,dtm,'sm','n')-1;
    if nargin>2
        tp = varargin{3};
    end
    %
    % _acceleration 0th order correction (pre-event)_
    %
    idx0 = vtm>=0&vtm<=tp;
    idx1 = vtm>=tp;
    idx2 = find(idx1>=1,1,'first');
    blc  = mean(tha(idx0));
    
    
    figure(1)
    plot(vtm,tha,'b'); hold all;
    hline(blc,'b--');
%     tha = tha-blc;
    tha = detrend(tha);
    tha(idx1) = cos_taper(tha(idx1));
    tha(~idx1) = 0;
    plot(vtm,tha,'r--');
    
    figure(2)
    thv = cumtrapz(tha)*dtm;
    plot(vtm,thv,'b');hold all;
    %
    % _quadratic fitting of velocity trace_
    %
    thv  = cumtrapz(tha)*dtm;
    
    plot(vtm,thv,'r--');hold all;
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
    pq = lsqlin(C,thv,[],[],Aeq,beq);
    thv = polyval(pq,vtm);
    plot(vtm,thv,'g');hold all;
    
    
%     pq  = polyfit(vtm(idx1),thv(idx1),2);
%     thv = polyval(pq,vtm);
    %
    % _subtracting derivative to 0th order corrected acceleration_
    %
    tha  = tha-avd_diff(dtm,thv);
    figure(1)
    plot(vtm,tha,'g');
    
    %
    % _integration_
    %
    [~,thv,thd] = integr_diff_avd(dtm,tha);
    figure(2)
    plot(vtm,thv,'m');
    
    figure(3)
    plot(vtm,thd);
    keyboard
    %% *OUTPUT*
    varargout{1} = tha(:);
    varargout{2} = thv(:);
    varargout{3} = thd(:);
    return
end
