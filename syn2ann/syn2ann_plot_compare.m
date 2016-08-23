function syn2ann_plot_compare(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl mrk mrka tit utd
    
    %% *SET-UP*
    flags = logical(varargin{1});
    xplot = varargin{2};
    identity = varargin{3};
    legplot = varargin{4};
    cpp = varargin{5};
    fn = varargin{6};
    time_shift = varargin{7};
    if nargin>7
        mrkd = varargin{8};
        lstd = varargin{9};
    end
    
    leg = {legplot};
    spg = numel(xplot);
    
    spg = [spg,1];
    pax = num2cell([1:spg(1),1:spg(1)]');
    
    if flags(1)
        %
        % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vTn;
            ypl{i_} = abs(xplot{i_}.syn{identity}.psa.(cpp))*utd.psa;
        end
        try
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
                'scl',scl.psa,'grd',grd.psa,'mrk',mrkd,'lst',lstd,...
                'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
                'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
                'leg',leg,'tit',tit.psa);
        catch
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
                'scl',scl.psa,'grd',grd.psa,...
                'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
                'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
                'leg',leg,'tit',tit.psa);
        end
        saveas(gcf,strcat(fn,sprintf('_psa_c_%s',cpp)),'epsc');
    end
    
    if flags(2)
        %
        % * _FOURIER SPECTRA_
        %
        for i_ = 1:spg(1)
            
            nfr     = numel(xplot{i_}.mon.vfr{identity})/2;
            xpl{i_} = xplot{i_}.mon.vfr{identity}(1:nfr);
            ypl{i_} = abs(xplot{i_}.syn{identity}.fsa.(cpp)(1:nfr))*utd.fsa;
        end
        try
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.fsa,...
                'scl',scl.fsa,'grd',grd.fsa,'mrk',mrkd,'lst',lstd,...
                'xlm',xlm.fsa,'xlb',xlb.fsa,'xtk',xtk.fsa,...
                'ylm',ylm.fsa,'ylb',ylb.fsa,'ytk',ytk.fsa,...
                'leg',leg,'tit',tit.fsa);
        catch
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.fsa,...
                'scl',scl.fsa,'grd',grd.fsa,...
                'xlm',xlm.fsa,'xlb',xlb.fsa,'xtk',xtk.fsa,...
                'ylm',ylm.fsa,'ylb',ylb.fsa,'ytk',ytk.fsa,...
                'leg',leg,'tit',tit.fsa);
        end
        saveas(gcf,strcat(fn,sprintf('_fsa_c_%s',cpp)),'epsc');
    end
    
    if any(flags(3:end))
        mrka = [repmat(mrk.tha,spg,1);repmat(mrk.pga,spg,1)];
    end
    if flags(3)
        %
        % * _ACCELERATION TIME HISTORY_
        %
        
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.tha.(cpp)*utd.tha;
        end
        xlbt = cell(spg(1),1);
        xlbt(end)= xlb.tha;
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp)(2)*utd.tha;
            if isempty(ylm.tha{i_})
                [~,ylm.tha{i_}]=get_axis_lim(xpl,ypl,0,1);
            end
            [~,ytk.tha{i_}]=get_axis_tick(xlm.tha{i_},ylm.tha{i_},10,diff(ylm.tha{1})/4);
        end
        
        %
        % # _ASIDE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.tha,'grd',grd.tha,...
            'xlm',xlm.tha,'xlb',xlbt,'xtk',xtk.tha,'mrk',mrka,'spg',spg,...
            'ylm',ylm.tha,'ylb',ylb.tha,'ytk',ytk.tha,'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_tha_s_%s',cpp)),'epsc');
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.tha,'grd',grd.tha,...
            'xlm',xlm.tha(end),'xlb',xlb.tha,'xtk',xtk.tha(end),...
            'ylm',ylm.tha(end),'ylb',ylb.tha,'ytk',ytk.tha(end),...
            'leg',leg,'tit',tit.tha);
        saveas(gcf,strcat(fn,sprintf('_tha_c_%s',cpp)),'epsc');
    end
    
    if flags(4)
        %
        % * _VELOCITY TIME HISTORY_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.thv.(cpp)*utd.thv;
        end
        xlbt = cell(spg(1),1);
        xlbt(end)= xlb.thv;
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp)(2)*utd.thv;
            if isempty(ylm.thv{i_})
                [~,ylm.thv{i_}]=get_axis_lim(xpl,ypl,0,1);
            end
            [~,ytk.thv{i_}]=get_axis_tick(xlm.thv{i_},ylm.thv{i_},10,diff(ylm.thv{1})/4);
        end
        %
        % # _ASIDE_
        %
        try
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.thv,'grd',grd.thv,...
            'xlm',xlm.thv(i_),'xlb',xlbt,'xtk',xtk.thv(i_),'mrk',mrka,'spg',spg,...
            'ylm',ylm.thv(i_),'ylb',ylb.thv,'ytk',ytk.thv(i_),'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_thv_s_%s',cpp)),'epsc');
        catch
            keyboard
        end
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thv,'grd',grd.thv,...
            'xlm',xlm.thv(end),'xlb',xlb.thv,'xtk',xtk.thv(end),...
            'ylm',ylm.thv(end),'ylb',ylb.thv,'ytk',ytk.thv(end),...
            'leg',leg,'tit',tit.thv);
        saveas(gcf,strcat(fn,sprintf('_thv_c_%s',cpp)),'epsc');
    end
    
    if flags(5)
        %
        % * _DISPLACEMENT TIME HISTORY_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.thd.(cpp)*utd.thd;
        end
        xlbt = cell(spg(1),1);
        xlbt(end) = xlb.thd;
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp)(2)*utd.thd;
            if isempty(ylm.thd{i_})
                [~,ylm.thd{i_}]=get_axis_lim(xpl,ypl,0,1);
            end
            [~,ytk.thd{i_}]=get_axis_tick(xlm.thd{i_},ylm.thd{i_},10,diff(ylm.thd{1})/4);
        end
        
        
        %
        % # _ASIDE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.thd,'grd',grd.thd,...
            'xlm',xlm.thd(i_),'xlb',xlbt,'xtk',xtk.thd(i_),'mrk',mrka,'spg',spg,...
            'ylm',ylm.thd(i_),'ylb',ylb.thd,'ytk',ytk.thd(i_),'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_thd_s_%s',cpp)),'epsc');
        %
        % * _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thd,'grd',grd.thd,...
            'xlm',xlm.thd(end),'xlb',xlb.thd,'xtk',xtk.thd(end),...
            'ylm',ylm.thd(end),'ylb',ylb.thd,'ytk',ytk.thd(end),...
            'leg',leg,'tit',tit.thd);
        saveas(gcf,strcat(fn,sprintf('_thd_c_%s',cpp)),'epsc');
    end
    close all;
    %
    return
end