ccc;plot_set_up;

batch_folder='/media/filippo/Data/Filippo/PHD_heavyweight/SPEED2D_NL/INPUTS/';
cd(batch_folder);
path_el=sprintf('%selastic/',batch_folder);
path_nl=sprintf('%splastic/',batch_folder);

lab={'dis','vel','acc'};
var={'d','v','a'}; nov=numel(var);
outvar=[]; nm=0;

for i=1:nov
    try
        if nm==0
            nm=numel(dir(sprintf('%smonitor*%s',path_el,var{i})));
        end
        outvar=[outvar;i];
    catch
    end
end
nov=numel(outvar);
xtt=0:5:20; xlt=[xtt(1),xtt(end)+5];
ytt=-1:.5:1; ylt=[ytt(1);ytt(end)];

for i=1:nm
    nl=struct('d',[],'v',[],'a',[]);
    el=struct('d',[],'v',[],'a',[]);
    sol=struct('el',el,'nl',nl);
    
    if i < 10
        file_name=sprintf('monitor0000%u.',i);
    elseif i < 100
        file_name=sprintf('monitor000%u.',i);
    elseif i < 1000
        file_name=sprintf('monitor00%u.',i);
    elseif i < 10000
        file_name=sprintf('monitor0%u.',i);
    elseif i < 100000
        file_name=sprintf('monitor%u.',i);
    end
    file_name_el=sprintf('%s%s',path_el,file_name);
    file_name_nl=sprintf('%s%s',path_nl,file_name);
    
    for j=1:1%nov
        
        eval(sprintf('sol.el.%s=importdata(''%s%s'');',var{outvar(j)},file_name_el,var{outvar(j)}));
        eval(sprintf('sol.nl.%s=importdata(''%s%s'');',var{outvar(j)},file_name_nl,var{outvar(j)}));
        
        fig=figure('name',sprintf('test_sp2_%s_%u',var{outvar(j)},i),...
            'visible','off','position',[1 1 21 14]);
        hax1=subplot(3,1,1,'parent',fig,...
            'xtick',xtt,'yticklabel',[]); hold(hax1,'all');
        hax2=subplot(3,1,2,'parent',fig,...
            'xtick',xtt,'yticklabel',[]); hold(hax2,'all');
        
        title(hax1,sprintf('%s-station: %u',lab{outvar(j)},i));
        xlabel(hax2,'t [s]','fontweight','bold');
        ylabel(hax1,sprintf('%s_X^{[1]}',upper(var{outvar(j)})),'fontweight','bold');
        ylabel(hax2,sprintf('%s_Y^{[1]}',upper(var{outvar(j)})),'fontweight','bold');
        try
            eval(sprintf(['plot(hax1,sol.el.%s(:,1),'...
                'normalize_vector(sol.el.%s(:,2)),''k'',''linewidth'',3);'],...
                var{outvar(j)},var{outvar(j)}));
            eval(sprintf(['plot(hax1,sol.nl.%s(:,1),'...
                'normalize_vector(sol.nl.%s(:,2)),''r'',''linewidth'',1.5);'],...
                var{outvar(j)},var{outvar(j)}));
            eval(sprintf(['plot(hax2,sol.el.%s(:,1),'...
                'normalize_vector(sol.el.%s(:,3)),''k'',''linewidth'',3);'],...
                var{outvar(j)},var{outvar(j)}));
            eval(sprintf(['plot(hax2,sol.nl.%s(:,1),'...
                'normalize_vector(sol.nl.%s(:,3)),''r'',''linewidth'',1.5);'],...
                var{outvar(j)},var{outvar(j)}));
        catch
            keyboard
        end
        saveas(fig,get(fig,'name'),'epsc');
        close(fig);
        
    end
end