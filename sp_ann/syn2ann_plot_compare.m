function syn2ann_plot_compare(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl mrk mrka tit
    
    %% *SET-UP*
    flags = logical(varargin{1});
    xplot = varargin{2};
    identity = varargin{3};
    legplot = varargin{4};
    cpp = varargin{5};
    fn = varargin{6};
    time_shift = varargin{7};
    
    leg = {legplot};
    spg = numel(xplot);
    mrka = [repmat(mrk.tha,spg,1);repmat(mrk.pga,spg,1)];
    spg = [spg,1];
    pax = num2cell([1:spg(1),1:spg(1)]');
    
    if flags(1)
        %
        % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vTn;
            ypl{i_} = abs(xplot{i_}.syn{identity}.psa.(cpp));
        end
        try
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
                'scl',scl.psa,'grd',grd.psa,...
                'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
                'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
                'leg',leg,'tit',tit.psa);
        catch
            keyboard
        end
        saveas(gcf,strcat(fn,sprintf('_psa_c_%s',cpp)),'epsc');
    end
    
    if flags(2)
        %
        % * _FOURIER SPECTRA_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vfr{identity};
            ypl{i_} = abs(xplot{i_}.syn{identity}.fsa.(cpp));
        end
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.fsa,...
            'scl',scl.fsa,'grd',grd.fsa,...
            'xlm',xlm.fsa,'xlb',xlb.fsa,'xtk',xtk.fsa,...
            'ylm',ylm.fsa,'ylb',ylb.fsa,'ytk',ytk.fsa,...
            'leg',leg,'tit',tit.fsa);
        saveas(gcf,strcat(fn,sprintf('_fsa_c_%s',cpp)),'epsc');
    end
    
    if flags(3)
        %
        % * _ACCELERATION TIME HISTORY_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.tha.(cpp);
        end
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pga.(cpp)(2);
        end
        
        [~,ylm.tha{identity}]=get_axis_lim(xpl(1:2),ypl(1:2),0,1);
        [~,ytk.tha{identity}]=get_axis_tick(xlm.tha{identity},ylm.tha{identity},10,diff(ylm.tha{1})/4);
        %
        % # _ASIDE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.tha,'grd',grd.tha,...
            'xlm',xlm.tha(identity),'xlb',xlb.tha,'xtk',xtk.tha(identity),'mrk',mrka,'spg',spg,...
            'ylm',ylm.tha(identity),'ylb',ylb.tha,'ytk',ytk.tha(identity),'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_tha_s_%s',cpp)),'epsc');
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.tha,'grd',grd.tha,...
            'xlm',xlm.tha(identity),'xlb',xlb.tha,'xtk',xtk.tha(identity),...
            'ylm',ylm.tha(identity),'ylb',ylb.tha,'ytk',ytk.tha(identity),...
            'leg',leg,'tit',tit.tha);
        saveas(gcf,strcat(fn,sprintf('_tha_c_%s',cpp)),'epsc');
    end
    
    if flags(4)
        %
        % * _VELOCITY TIME HISTORY_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.thv.(cpp);
        end
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgv.(cpp)(2);
        end
        
        [~,ylm.thv{identity}]=get_axis_lim(xpl(1:2),ypl(1:2),0,1);
        [~,ytk.thv{identity}]=get_axis_tick(xlm.thv{identity},ylm.thv{identity},10,diff(ylm.thv{1})/4);
        %
        % # _ASIDE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.thv,'grd',grd.thv,...
            'xlm',xlm.thv(identity),'xlb',xlb.thv,'xtk',xtk.thv(identity),'mrk',mrka,'spg',spg,...
            'ylm',ylm.thv(identity),'ylb',ylb.thv,'ytk',ytk.thv(identity),'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_thv_s_%s',cpp)),'epsc');
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thv,'grd',grd.thv,...
            'xlm',xlm.thv(identity),'xlb',xlb.thv,'xtk',xtk.thv(identity),...
            'ylm',ylm.thv(identity),'ylb',ylb.thv,'ytk',ytk.thv(identity),...
            'leg',leg,'tit',tit.thv);
        saveas(gcf,strcat(fn,sprintf('_thv_c_%s',cpp)),'epsc');
    end
    
    if flags(5)
        %
        % * _DISPLACEMENT TIME HISTORY_
        %
        for i_ = 1:spg(1)
            xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
            ypl{i_} = xplot{i_}.syn{identity}.thd.(cpp);
        end
        for i_ = 1:spg(1)
            xpl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp)(1)-time_shift(i_);
            ypl{spg(1)+i_} = xplot{i_}.syn{identity}.pgd.(cpp)(2);
        end
        
        [~,ylm.thd{identity}]=get_axis_lim(xpl(1:2),ypl(1:2),0,1);
        [~,ytk.thd{identity}]=get_axis_tick(xlm.thv{identity},ylm.thv{identity},10,diff(ylm.thd{1})/4);
        %
        % # _ASIDE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.s,...
            'pax',pax,'scl',scl.thd,'grd',grd.thd,...
            'xlm',xlm.thd(identity),'xlb',xlb.thd,'xtk',xtk.thd(identity),'mrk',mrka,'spg',spg,...
            'ylm',ylm.thd(identity),'ylb',ylb.thd,'ytk',ytk.thd(identity),'tit',legplot);
        saveas(gcf,strcat(fn,sprintf('_thd_s_%s',cpp)),'epsc');
        %
        % * _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thd,'grd',grd.thd,...
            'xlm',xlm.thd(identity),'xlb',xlb.thd,'xtk',xtk.thd(identity),...
            'ylm',ylm.thd(identity),'ylb',ylb.thd,'ytk',ytk.thd(identity),...
            'leg',leg,'tit',tit.thd);
        saveas(gcf,strcat(fn,sprintf('_thd_c_%s',cpp)),'epsc');
    end
    close all;
    %
    return
end