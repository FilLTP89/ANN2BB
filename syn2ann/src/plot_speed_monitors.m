%%  *Plot SPEED monitor files*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% plot_speed_monitors: function to plot SPEED output files
%% INPUT:
%          _mon (monitor structure)_
%          mon.pth = path to monitor files (string)
%          mon.idt = monitor identity  (integer vector)
%          mon.ni  = number of monitors
%          mon.var = monitor variables (string array)
%          mon.nv  = number of variables
%          mon.mcp = component (string array)
%          mon.nc  = number of components
%% OUTPUT:
%% N.B.
% Need for _speed_monitor_name.m,fpplot.m_
function plot_speed_monitors(varargin)
    %% SET-UP
    %%
    % _monitor
    mon.pth = varargin{1};
    mon.idt = varargin{2};
    mon.var = varargin{3};
    mon.mcp = varargin{4};
    mon.ni = numel(mon.idt);
    mon.nv = numel(mon.var);
    for j_ = 1:mon.nv
        mon.nc(j_) = numel(mon.mcp{j_});
    end
    
    cl3 = [{'x'};{'y'};{'z'}];
    cl6 = [{'xx'};{'yy'};{'xy'};{'zz'};{'xz'};{'yz'}];
    
    %% PLOT MONITOR RESULTS
    for j_ = 1:mon.nv
        % vector/tensor components
        switch mon.var{j_}
            case {'a','v','d'}
                clb = cl3;
            case {'e','s'}
                clb = cl6;
        end
        for i_ = 1:mon.ni
            dat = importdata(speed_monitor_name(mon.idt(i_),mon.var{j_},...
                mon.pth));
            xpl = repmat({dat(:,1)},mon.nc(j_),1);
            ypl = cell(size(xpl));
            pax = cell(size(xpl));
            leg = cell(size(xpl));
            for kk_ = 1:mon.nc(j_)
                k_ = mon.mcp{j_}(kk_);
                ypl{kk_} = dat(:,k_+1);
                pax{kk_} = kk_;
                leg{kk_} = strcat(num2str(mon.idt(i_)),mon.var{j_},clb{k_});
            end
            [hfg,~,~] = fpplot('nfg',sprintf('%u_%s',mon.idt(i_),mon.var{j_}),...
                'pfg',[0 0 20 16],'xlb',{'t [s]'},'ylb',clb(mon.mcp{j_}),...
                'xpl',xpl,'ypl',ypl,'pax',pax,'leg',leg);
        end
        
    end
    return
end
