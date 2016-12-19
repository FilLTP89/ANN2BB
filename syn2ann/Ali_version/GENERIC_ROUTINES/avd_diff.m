function [varargout] = avd_diff(varargin)
    dtm = varargin{1};
    ths = varargin{2}(:);
    ntm = numel(ths);
    
    thd(2:ntm-1,1) = (ths(3:ntm,1)-ths(1:ntm-2,1))./(2*dtm);
    thd([1,ntm],1) = [0.0;thd(ntm-1,1)];
    varargout{1} = thd(:);
    return
end