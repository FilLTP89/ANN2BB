%==== format figures


function format_figures_multi(hax)
    
    if nargin==0
        hax=gca;
    end
    
    set(hax,'box','on','xminortick','on','yminortick','on');
    %
    set(hax,'ticklength',[.02,.02]);
    ti = get(hax,'title');
    set(ti,'fontsize',20,'fontweight','bold');%25,30
    xl = get(hax,'xlabel');
    set(xl,'fontsize',16,'fontweight','bold');%20,24
    yl = get(hax,'ylabel');
    set(yl,'fontsize',16,'fontweight','bold');%,20,24
    zl = get(hax,'zlabel');
    set(zl,'fontsize',16);%20,24
    li = get(hax,'children');
    set(li,'linewidth',2,'markersize',3);%2.5-3
    set(hax,'linewidth',1.6,'fontsize',12);   %2,2-3,20
   
end