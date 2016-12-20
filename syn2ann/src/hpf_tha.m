function [varargout] = hpf_tha(varargin)
    %% *SETUP*
    dtm = varargin{1};
    tha = varargin{2}(:);
    tp = PphasePicker(tha,dtm,'sm','n')-1;
    if nargin>2
        tp = varargin{3};
    end
    %
    % _base-line correction
    %
    [tha,~,~] = blc_tha(dtm,tha,tp);
    %
    % _4th-order Butterworth causal filtering_
    %
    [bfb,bfa,~] = create_butter_filter(4,0.02,[],0.5/dtm);
    tha = filter(bfb,bfa,tha);
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
