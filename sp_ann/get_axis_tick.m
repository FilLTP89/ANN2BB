function [xtk,ytk]=get_axis_tick(xlm,ylm,dx,dy)
    
    % x-tick
    if abs(xlm(1))==abs(xlm(2))
        xtk = xlm(1):dx:0;
        xtk = [xtk,dx:dx:xlm(2)];
    else
        xtk = xlm(1):dx:xlm(2);
    end
    
    if abs(diff(xlm)/dx-round(diff(xlm)/dx))>dx/100;
        xtk(end+1) = xlm(2);
    end
    
    % y-tick
    if abs(ylm(1))==abs(ylm(2))
        ytk = ylm(1):dy:0;
        ytk = [ytk,dy:dy:ylm(2)];
    else
        ytk = ylm(1):dy:ylm(2);
    end
    
    if abs(diff(ylm)/dy-round(diff(ylm)/dy))>dy/100;
        ytk(end+1) = ylm(2);
    end
    
    return
end

