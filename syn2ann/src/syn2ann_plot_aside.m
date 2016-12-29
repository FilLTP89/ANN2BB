function syn2ann_plot_aside(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl utd
    
    %% *SET-UP*
    flags = logical(varargin{1});
    xplot = varargin{2};
    identity = varargin{3};
    dlg = varargin{4};
    legplot = varargin{5};
    cpp = varargin{6};
    std = varargin{7};
    fnm = varargin{8};
    time_shift = varargin{9};
    if nargin>9
        mrk = varargin{10};
        lst = varargin{11};
    end
    
    spg = numel(xplot);
    plt.xpl = cell(2*spg,1);
    plt.ypl = cell(2*spg,1);
    %% *FOURIER/PSA-SPECTRA*
    if flags(1)
        leg = cell(2*spg,1);
        plt.mrk = cell(2*spg,1);
        plt.lst = cell(2*spg,1);
        % _fourier's spectra_
        for i_ = 1:spg
            nfr     = numel(xplot{i_}.mon.vfr{identity})/2;
            % x-plot
            plt.xpl{i_,1} = xplot{i_}.mon.vfr{identity}(1:nfr);
            % y-plot
            plt.ypl{i_,1} = abs(xplot{i_}.syn{identity}.fsa.(cpp.nm)(1:nfr))*utd.fsa;
            % legend
            leg(i_,1) = {legplot{i_}};
            plt.mrk{i_,1} = mrk{i_};
            plt.lst{i_,1} = lst{i_};
        end
        % psa spectra
        for i_ = 1:spg
            % x-plot
            plt.xpl{i_+spg,1} = xplot{i_}.mon.vTn;
            % y-plot
            plt.ypl{i_+spg,1} = abs(xplot{i_}.syn{identity}.psa.(cpp.nm))*utd.psa;
            leg(i_+spg,1) = {legplot{i_}};
            plt.mrk{i_+spg,1} = mrk{i_};
            plt.lst{i_+spg,1} = lst{i_};
        end
        % others
        plt.pax = num2cell([1;1;1;2;2;2]);
        plt.leg = cell(2,1);
        plt.xlm = cell(2,1);
        plt.xtk = cell(2,1);
        plt.xlb = cell(2,1);
        plt.ylm = cell(2,1);
        plt.ytk = cell(2,1);
        plt.ylb = cell(2,1);
        plt.scl = cell(2,1);
        plt.grd = cell(2,1);
        plt.tit = cell(2,1);
        plt.scl(1,1) = scl.fsa(1);
        plt.grd(1,1) = grd.fsa(1);
        plt.tit(1,1) = {strcat(std,'-',dlg)};
        plt.xlm(1,1) = xlm.fsa(1);
        plt.xtk(1,1) = xtk.fsa(1);
        plt.xlb(1,1) = xlb.fsa(1);
        plt.ylm(1,1) = ylm.fsa(1);
        plt.ytk(1,1) = ytk.fsa(1);
        plt.ylb(1,1) = ylb.fsa(1);
        plt.scl(2,1) = scl.psa(1);
        plt.grd(2,1) = grd.psa(1);
        plt.tit(2,1) = {strcat(std,'-',dlg)};
        plt.xlm(2,1) = xlm.psa(1);
        plt.xtk(2,1) = xtk.psa(1);
        plt.xlb(2,1) = xlb.psa(1);
        plt.ylm(2,1) = ylm.psa(1);
        plt.ytk(2,1) = ytk.psa(1);
        plt.ylb(2,1) = ylb.psa(1);
        plt.leg(1,1) = {leg(1:spg,1)};
        plt.leg(2,1) = {leg(spg+1:end,1)};
        plt.spg = [1,2];
        
        fpplot('xpl',plt.xpl,'ypl',plt.ypl,...
            'pfg',pfg.fsp,'spg',plt.spg,'pax',plt.pax,...
            'scl',plt.scl,'grd',plt.grd,'tit',plt.tit,...
            'xlm',plt.xlm,'xtk',plt.xtk,'xlb',plt.xlb,...
            'ylm',plt.ylm,'ytk',plt.ytk,'ylb',plt.ylb,...
            'leg',plt.leg,'mrk',plt.mrk,'lst',plt.lst);
        saveas(gcf,strcat(fnm,sprintf('_fsp_%s',cpp.nm)),'jpg');
    end
    
    %% *TIME HISTORIES*
    if flags(2)
        % others
        plt.spg = [spg,2];
        plt.pax = num2cell([1;3;5;2;4;6]);
        plt.leg = cell(2*spg,1);
        plt.xlm = cell(2*spg,1);
        plt.xtk = cell(2*spg,1);
        plt.xlb = cell(2*spg,1);
        plt.ylm = cell(2*spg,1);
        plt.ytk = cell(2*spg,1);
        plt.ylb = cell(2*spg,1);
        plt.scl = cell(2*spg,1);
        plt.grd = cell(2*spg,1);
        plt.tit = cell(2*spg,1);
        plt.leg = cell(2*spg,1);
        % _velocities_
        for i_ = 1:spg
            % x-plot
            plt.xpl{i_,1} = xplot{i_}.mon.vtm{identity}(:)-time_shift(i_);
            plt.xlm(i_,1) = xlm.thv(cpp.nb);
            plt.xtk(i_,1) = xtk.thv(cpp.nb);
            plt.xlb(spg,1) = xlb.thv(1);
            % y-plot
            plt.ypl{i_,1} = xplot{i_}.syn{identity}.thv.(cpp.nm)(:)*utd.thv;
            plt.ylm(i_,1) = ylm.thv(cpp.nb);
            plt.ytk(i_,1) = ytk.thv(cpp.nb);
            plt.ylb(i_,1) = ylb.thv(cpp.nb);
            % others
            plt.scl(i_,1) = scl.thv(1);
            plt.grd(i_,1) = grd.thv(1);
            plt.tit(i_,1) = {legplot{i_}};
            plt.leg(i_,1) = {strcat(std,'-',dlg)};
        end
        % accelerations
        for i_ = 1:spg
            % x-plot
            plt.xpl{i_+spg,1} = xplot{i_}.mon.vtm{identity}(:)-time_shift(i_);
            plt.xlm(i_+spg,1) = xlm.tha(cpp.nb);
            plt.xtk(i_+spg,1) = xtk.tha(cpp.nb);
            plt.xlb(spg+spg,1) = xlb.tha(1);
            % y-plot
            plt.ypl{i_+spg,1} = xplot{i_}.syn{identity}.tha.(cpp.nm)(:)*utd.tha;
            plt.ylm(i_+spg,1) = ylm.tha(cpp.nb);
            plt.ytk(i_+spg,1) = ytk.tha(cpp.nb);
            plt.ylb(i_+spg,1) = ylb.tha(1);
            % others
            plt.scl(i_+spg,1) = scl.tha(1);
            plt.grd(i_+spg,1) = grd.tha(1);
            plt.tit(i_+spg,1) = {legplot{i_}};
            plt.leg(i_+spg,1) = {strcat(std,'-',dlg)};
        end
        
        fpplot('xpl',plt.xpl,'ypl',plt.ypl,...
            'pfg',pfg.fth,'spg',plt.spg,'pax',plt.pax,...
            'scl',plt.scl,'grd',plt.grd,...
            'xlm',plt.xlm,'xtk',plt.xtk,'xlb',plt.xlb,...
            'ylm',plt.ylm,'ytk',plt.ytk,'ylb',plt.ylb,...
            'tit',plt.tit,'leg',plt.leg);
        saveas(gcf,strcat(fnm,sprintf('_fth_%s',cpp.nm)),'jpg');
    end
    close all;
    %
    return
end