function syn2ann_plot_aside(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl mrk mrka tit utd
    
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
    
    spg = numel(xplot);
    
    plt.spg = [spg,2];
    plt.pax = num2cell([1;3;5;2;4;6]);
    
    plt.xpl = cell(2*spg(1),1);
    plt.xlm = cell(2*spg(1),1);
    plt.xtk = cell(2*spg(1),1);
    plt.xlb = cell(2*spg(1),1);
    plt.ypl = cell(2*spg(1),1);
    plt.ylm = cell(2*spg(1),1);
    plt.ytk = cell(2*spg(1),1);
    plt.ylb = cell(2*spg(1),1);
    plt.scl = cell(2*spg(1),1);
    plt.grd = cell(2*spg(1),1);
%     plt.leg = cell(2*spg(1),1);
    plt.tit = cell(2*spg(1),1);
    
    %% *FOURIER/PSA-SPECTRA*
    if flags(1)
        % _fourier's spectra_
        for i_ = 1:plt.spg(1)
            nfr     = numel(xplot{i_}.mon.vfr{identity})/2;
            % x-plot
            plt.xpl{i_,1} = xplot{i_}.mon.vfr{identity}(1:nfr);
            plt.xlm(i_,1) = xlm.fsa(1);
            plt.xtk(i_,1) = xtk.fsa(1);
            plt.xlb(i_,1) = xlb.fsa(1);
            % y-plot
            plt.ypl{i_,1} = abs(xplot{i_}.syn{identity}.fsa.(cpp.nm)(1:nfr))*utd.fsa;
            plt.ylm(i_,1) = ylm.fsa(1);
            plt.ytk(i_,1) = ytk.fsa(1);
            plt.ylb(i_,1) = ylb.fsa(1);
            % others
            plt.scl(i_,1) = scl.fsa(1);
            plt.grd(i_,1) = grd.fsa(1);
            plt.tit(i_,1) = {legplot{i_}};
            plt.leg(i_,1) = {strcat(std,'-',dlg)};
        end
        % psa spectra
        for i_ = 1:plt.spg(1)
            % x-plot
            plt.xpl{i_+plt.spg(1),1} = xplot{i_}.mon.vTn;
            plt.xlm(i_+plt.spg(1),1) = xlm.psa(1);
            plt.xtk(i_+plt.spg(1),1) = xtk.psa(1);
            plt.xlb(i_+plt.spg(1),1) = xlb.psa(1);
            % y-plot
            plt.ypl{i_+plt.spg(1),1} = abs(xplot{i_}.syn{identity}.psa.(cpp.nm))*utd.psa;
            plt.ylm(i_+plt.spg(1),1) = ylm.psa(1);
            plt.ytk(i_+plt.spg(1),1) = ytk.psa(1);
            plt.ylb(i_+plt.spg(1),1) = ylb.psa(1);
            % others
            plt.scl(i_+plt.spg(1),1) = scl.psa(1);
            plt.grd(i_+plt.spg(1),1) = grd.psa(1);
            plt.tit(i_+plt.spg(1),1) = {legplot{i_}};
            plt.leg(i_+plt.spg(1),1) = {strcat(std,'-',dlg)};
        end

        
        fpplot('xpl',plt.xpl,'ypl',plt.ypl,...
            'pfg',pfg.fsp,'spg',plt.spg,'pax',plt.pax,...
            'scl',plt.scl,'grd',plt.grd,...
            'xlm',plt.xlm,'xtk',plt.xtk,'xlb',plt.xlb,...
            'ylm',plt.ylm,'ytk',plt.ytk,'ylb',plt.ylb,...
            'tit',plt.tit,'leg',plt.leg);
        saveas(gcf,strcat(fnm,sprintf('_fsp_%s',cpp.nm)),'jpg');
    end
    
    %% *TIME HISTORIES*
    if flags(2)
        % _velocities_
        for i_ = 1:plt.spg(1)
            % x-plot
            plt.xpl{i_,1} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            plt.xlm(i_,1) = xlm.thv(cpp.nb);
            plt.xtk(i_,1) = xtk.thv(cpp.nb);
            plt.xlb(i_,1) = xlb.thv(1);
            % y-plot
            plt.ypl{i_,1} = xplot{i_}.syn{identity}.thv.(cpp.nm)*utd.thv;
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
        for i_ = 1:plt.spg(1)
            % x-plot
            plt.xpl{i_+plt.spg(1),1} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            plt.xlm(i_+plt.spg(1),1) = xlm.tha(cpp.nb);
            plt.xtk(i_+plt.spg(1),1) = xtk.tha(cpp.nb);
            plt.xlb(i_+plt.spg(1),1) = xlb.tha(1);
            % y-plot
            plt.ypl{i_+plt.spg(1),1} = xplot{i_}.syn{identity}.tha.(cpp.nm)*utd.tha;
            plt.ylm(i_+plt.spg(1),1) = ylm.tha(cpp.nb);
            plt.ytk(i_+plt.spg(1),1) = ytk.tha(cpp.nb);
            plt.ylb(i_+plt.spg(1),1) = ylb.tha(1);
            % others
            plt.scl(i_+plt.spg(1),1) = scl.tha(1);
            plt.grd(i_+plt.spg(1),1) = grd.tha(1);
            plt.tit(i_+plt.spg(1),1) = {legplot{i_}};
            plt.leg(i_+plt.spg(1),1) = {strcat(std,'-',dlg)};
        end

        
        fpplot('xpl',plt.xpl,'ypl',plt.ypl,...
            'pfg',pfg.fth,'spg',plt.spg,'pax',plt.pax,...
            'scl',plt.scl,'grd',plt.grd,...
            'xlm',plt.xlm,'xtk',plt.xtk,'xlb',plt.xlb,...
            'ylm',plt.ylm,'ytk',plt.ytk,'ylb',plt.ylb,...
            'tit',plt.tit,'leg',plt.leg);
        saveas(gcf,strcat(fnm,sprintf('_fth_%s',cpp.nm)),'jpg');
    end
    
    %     if any(flags(3:end))
    %         mrka = [repmat(mrk.tha,[spg,1]);repmat(mrk.pga,[spg,1])];
    %     end
    %     if flags(3)
    %         %
    %         % * _ACCELERATION TIME HISTORY_
    %         %
    %         for i_ = 1:spg(1)
    %             
    %         end
    %         xlbt = cell(spg(1),1);
    %         xlbt(end)= xlb.tha;
    %         for i_ = 1:spg(1)
    %             xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp.nm)(1)-time_shift(i_);
    %             ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp.nm)(2)*utd.tha;
    %             [~,ytk.tha{i_}]=get_axis_tick(xlm.tha{i_},ylm.tha{i_},10,diff(ylm.tha{1})/4);
    %         end
    %         %
    %         % # _COMPARE_
    %         %
    %         fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
    %             'scl',scl.tha,'grd',grd.tha,...
    %             'xlm',xlm.tha(end),'xlb',xlb.tha,'xtk',xtk.tha(end),...
    %             'ylm',ylm.tha(end),'ylb',ylb.tha,'ytk',ytk.tha(end),...
    %             'leg',leg,'tit',strcat(std,'-',dlg));
    %         saveas(gcf,strcat(fnm,sprintf('_tha_c_%s',cpp)),'jpg');
    %     end
    %
    %     if flags(4)
    %         %
    %         % * _VELOCITY TIME HISTORY_
    %         %
    %         for i_ = 1:spg(1)
    %             xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    %             ypl{i_} = xplot{i_}.syn{identity}.thv.(cpp.nm)*utd.thv;
    %         end
    %         xlbt = cell(spg(1),1);
    %         xlbt(end)= xlb.thv;
    %         for i_ = 1:spg(1)
    %             xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp.nm)(1)-time_shift(i_);
    %             ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp.nm)(2)*utd.thv;
    %
    %         end
    %         %
    %         % # _COMPARE_
    %         %
    %         fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
    %             'scl',scl.thv,'grd',grd.thv,...
    %             'xlm',xlm.thv(end),'xlb',xlb.thv,'xtk',xtk.thv(end),...
    %             'ylm',ylm.thv(end),'ylb',ylb.thv,'ytk',ytk.thv(end),...
    %             'leg',leg,'tit',strcat(std,'-',dlg));
    %         saveas(gcf,strcat(fnm,sprintf('_thv_c_%s',cpp)),'jpg');
    %     end
    %
    %     if flags(5)
    %         %
    %         % * _DISPLACEMENT TIME HISTORY_
    %         %
    %         for i_ = 1:spg(1)
    %             xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    %             ypl{i_} = xplot{i_}.syn{identity}.thd.(cpp.nm)*utd.thd;
    %         end
    %         xlbt = cell(spg(1),1);
    %         xlbt(end) = xlb.thd;
    %         for i_ = 1:spg(1)
    %             xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp.nm)(1)-time_shift(i_);
    %             ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp.nm)(2)*utd.thd;
    %             if isempty(ylm.thd{i_})
    %                 [~,ylm.thd{i_}]=get_axis_lim(xpl,ypl,0,1);
    %             end
    %             [~,ytk.thd{i_}]=get_axis_tick(xlm.thd{i_},ylm.thd{i_},10,diff(ylm.thd{1})/4);
    %         end
    %         %
    %         % * _COMPARE_
    %         %
    %         fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
    %             'scl',scl.thd,'grd',grd.thd,...
    %             'xlm',xlm.thd(end),'xlb',xlb.thd,'xtk',xtk.thd(end),...
    %             'ylm',ylm.thd(end),'ylb',ylb.thd,'ytk',ytk.thd(end),...
    %             'leg',leg,'tit',strcat(std,'-',dlg));
    %         saveas(gcf,strcat(fnm,sprintf('_thd_c_%s',cpp)),'jpg');
    %     end
    close all;
    %
    return
end