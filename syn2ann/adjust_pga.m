function [varargout] = adjust_pga(varargin)
    warning off
    dtm = varargin{1};
    tha = varargin{2};
    pga_target = varargin{3};
    scl = varargin{4};
    pga_target = pga_target*scl;
    
    [pga_ref,idx] = max(abs(tha));
    ntm = numel(tha);
    vtm = (0:ntm-1)*dtm;
    dt_cor = 0.05;
    nc = round(dt_cor./dtm);
    if mod(nc,2)==0
        nc = nc+1;
    end
    
    x = [vtm(idx-(nc-1)/2);vtm(idx);vtm(idx+(nc-1)/2)];
    y = [tha(idx-(nc-1)/2);pga_target.*sign(tha(idx));tha(idx+(nc-1)/2)];
    xi = vtm(idx-nc/2:idx+nc/2-1);
    yi_c = interp1(x,y,xi,'cubic');
    %      yi_p = interp1(x,y,xi,'pchip');
    %      yi_l = interp1(x,y,xi,'linear');
    
    tha(idx-(nc-1)/2:idx+(nc-1)/2) = yi_c;
    
    if (max(abs(tha))-pga_target)>1e-3
        % repeat the same operation as before
        [~,idx] = max(abs(tha));
        x = [vtm(idx-(nc-1)/2);vtm(idx);vtm(idx+(nc-1)/2)];
        y = [tha(idx-(nc-1)/2);pga_target.*sign(tha(idx));tha(idx+(nc-1)/2)];
        xi = vtm(idx-nc/2:idx+nc/2-1);
        yi_c = interp1(x,y,xi,'cubic');
        %      yi_p = interp1(x,y,xi,'pchip');
        %      yi_l = interp1(x,y,xi,'linear');
        tha(idx-(nc-1)/2:idx+(nc-1)/2) = yi_c; 
    end
    if abs((max(abs(tha))-pga_target)/(pga_target))>5e-2
        disp('Check PGA adjustment!');
    end
    varargout{1} = tha;
    return
end