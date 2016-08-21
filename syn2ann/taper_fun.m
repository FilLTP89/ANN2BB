% Hann cosinus-function tapering
% A = 0.5
% B = 0.5
% w = A + B * cos(ind)
% ind = (ntm - n_) * pi / (ntm - 1) --->end
% ind = n_ * pi / (ntm - 1)         --->start
% ntm tapering points


function [varargout] = taper_fun(varargin)
    
    ths = varargin{1}(:);
    pct = varargin{2};
    idx.before = varargin{3};
    idx.after = varargin{4};
    Ntot=length(ths);
    ntm = round(Ntot * (pct / 100.));
    
    if logical(idx.before)
        for n_ = 1:ntm
            dth = (ntm - n_) * pi / (ntm - 1);
            a = 0.5 * (1.0 + cos(dth));
            ths(n_) = ths(n_) * a ;
        end
    end
    
    if logical(idx.after)
        for n_ = 0:ntm-1
            dth = n_ * pi / (ntm - 1);
            a = 0.5 * (1.0 + cos(dth));
            bb = Ntot - (ntm - 1) + n_;
            ths(bb) = ths(bb) * a;
        end
    end
    
    varargout{1} = ths;
    
    return
end

