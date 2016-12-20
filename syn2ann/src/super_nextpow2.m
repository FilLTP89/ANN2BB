function [varargout] = super_nextpow2(varargin)
    n = varargin{1};
    
    if mod(n,2)
        p = nextpow2(n);
    else
        p = nextpow2(n+1);
    end
    
    varargout{1} = pow2(p);
    return
end