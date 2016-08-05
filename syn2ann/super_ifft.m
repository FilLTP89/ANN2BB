function [varargout] = super_ifft(varargin)
    dtm = varargin{1};
    ntm = varargin{2};
    fss = varargin{3};
    
    nfr = numel(fss);
    fst = zeros(nfr,1);
    
    fst(1:nfr/2) = fss(1:nfr/2);
    fst(nfr/2+2:nfr) = flip(conj(fss(2:nfr/2)));
    ths = real(ifft(fst))./dtm;
%     ths = ths(1:ntm);
    varargout{1} = ths(:);
    return
end