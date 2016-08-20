function [varargout] = adjust_pga(varargin)
    %% *SETUP*
    dtm        = varargin{1};
    tha        = varargin{2};
    tar_pga    = abs(varargin{3});
    
    ntm = numel(tha);
    vtm = (0:ntm-1)*dtm;
    [ref_pga,ipga] = max(abs(tha));
    
    nc_pga = 11;
    ic_pga0 = ipga+[-(nc_pga-1)/2;0;(nc_pga-1)/2];
    ic_pga1 = ipga+(-(nc_pga-1)/2:(nc_pga-1)/2)';
    sc_pga  = [1;tar_pga/ref_pga;1];
    
    tha(ic_pga1) = interp1(vtm(ic_pga0),...
        tha(ic_pga0).*sc_pga,vtm(ic_pga1),'pchip');
    
    [ref_pga,ipga] = max(abs(tha));
    
    if abs(ref_pga-tar_pga)>1e-3
        ic_pga0 = ipga+[-(nc_pga-1)/2;0;(nc_pga-1)/2];
        ic_pga1 = ipga+(-(nc_pga-1)/2:(nc_pga-1)/2)';
        sc_pga  = [1;tar_pga/ref_pga;1];
        tha(ic_pga1) = interp1(vtm(ic_pga0),...
            tha(ic_pga0).*sc_pga,vtm(ic_pga1),'pchip');
        ref_pga = max(abs(tha));
    end
       
    % check final PGA
    if abs(ref_pga-tar_pga)/tar_pga>5e-2
        disp('Check PGA adjustment!');
    end
    varargout{1} = tha(:);
    return
end