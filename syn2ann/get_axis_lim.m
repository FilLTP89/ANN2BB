function [xlm,ylm]=get_axis_lim(xplt,yplt,sx,sy,varargin)
    
    
    if iscell(xplt)
        xpl=[];
        for i_=1:numel(xplt)
            xpl=[xpl;xplt{i_}(:)];
        end
    else
        xpl = xplt(:);
    end
    
    if iscell(yplt)
        ypl=[];
        for i_=1:numel(yplt)
            ypl=[ypl;yplt{i_}(:)];
        end
    else
        xpl = xplt(:);
    end
    
    xlm_max = (max(xpl));
    ylm_max = (max(ypl));
    xlm_min = (min(xpl));
    ylm_min = (min(ypl));
    if numel(varargin)>0
        ordxx = varargin{1};
        ordyy = varargin{2};
    else
        ordxx = max(abs([ceil(log10(abs(xlm_max)));ceil(log10(abs(xlm_min)))]));
        ordyy = max(abs([ceil(log10(abs(ylm_max)));ceil(log10(abs(ylm_min)))]));
        if isnan(ordxx)
            ordxx=4;
        end
        if isnan(ordyy)
            ordyy=4;
        end
    end
    ordx = 10^ordxx;
    ordy = 10^ordyy;
    
    if xlm_max>=1
        xlm_max=ceil(xlm_max);
    else
        xlm_max=round(ordx*rem(xlm_max,ordx))/ordx;
    end
    
    if ylm_max>=1
        ylm_max=ceil(ylm_max);
    else
        ylm_max=round(ordy*rem(ylm_max,ordy))/ordy;
    end
    
    if abs(xlm_min)>=1
        xlm_min=floor(xlm_min);
    else
        xlm_min=floor(ordx*rem(xlm_min,ordx))/ordx;
    end
    
    if abs(ylm_min)>=1
        ylm_min=floor(ylm_min);
    else
        ylm_min=floor(ordy*rem(ylm_min,ordy))/ordy;
    end
    
    xlm_abs_max = max([abs(xlm_min),abs(xlm_max)]);
    ylm_abs_max = max([abs(ylm_min),abs(ylm_max)]);
    
    if logical(sx)
        xlm = [-xlm_abs_max,xlm_abs_max];
    else
        xlm = [xlm_min,xlm_max];
    end
    if logical(sy)
        ylm = [-ylm_abs_max,ylm_abs_max];
    else
        ylm = [ylm_min,ylm_max];
    end
    return
end

