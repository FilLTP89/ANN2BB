function [varargout] = super_ifft(varargin)
    %% *SETUP*
    dtm = varargin{1};
    fss = varargin{2}(:);
    if nargin>2
        error('CHANGE SYNTAX');
    end
    
    nfr = numel(fss);
    ths = ifft(fss,nfr,1,'symmetric')/dtm;
    
    %% *OUTPUT*
    varargout{1} = ths(:);
    return
end