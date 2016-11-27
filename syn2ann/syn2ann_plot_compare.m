function syn2ann_plot_compare(varargin)
    global pfg xlm xlb xtk ylm ylb ytk grd scl mrk mrka utd tit
    
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
        mrkd = varargin{10};
        lstd = varargin{11};
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
                'leg',leg,'tit',strcat(std,'-',dlg));
        catch
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
                'scl',scl.psa,'grd',grd.psa,...
                'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
                'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
                'leg',leg,'tit',strcat(std,'-',dlg));
        end
        saveas(gcf,strcat(fnm,sprintf('_psa_%s',cpp)),'epsc');
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
                'leg',leg,'tit',strcat(std,'-',dlg));
        catch
            fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.fsa,...
                'scl',scl.fsa,'grd',grd.fsa,...
                'xlm',xlm.fsa,'xlb',xlb.fsa,'xtk',xtk.fsa,...
                'ylm',ylm.fsa,'ylb',ylb.fsa,'ytk',ytk.fsa,...
                'leg',leg,'tit',strcat(std,'-',dlg));
        end
        saveas(gcf,strcat(fnm,sprintf('_fsa_%s',cpp)),'epsc');
    end
    
    if any(flags(3:end))
        mrka = [repmat(mrk.tha,[spg,1]);repmat(mrk.pga,[spg,1])];
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
            [~,ytk.tha{i_}]=get_axis_tick(xlm.tha{i_},ylm.tha{i_},10,diff(ylm.tha{1})/4);
        end
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.tha,'grd',grd.tha,...
            'xlm',xlm.tha(end),'xlb',xlb.tha,'xtk',xtk.tha(end),...
            'ylm',ylm.tha(end),'ylb',ylb.tha,'ytk',ytk.tha(end),...
            'leg',leg,'tit',strcat(std,'-',dlg));
        saveas(gcf,strcat(fnm,sprintf('_tha_%s',cpp)),'epsc');
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
            
        end
        %
        % # _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thv,'grd',grd.thv,...
            'xlm',xlm.thv(end),'xlb',xlb.thv,'xtk',xtk.thv(end),...
            'ylm',ylm.thv(end),'ylb',ylb.thv,'ytk',ytk.thv(end),...
            'leg',leg,'tit',strcat(std,'-',dlg));
        saveas(gcf,strcat(fnm,sprintf('_thv_%s',cpp)),'epsc');
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
        % * _COMPARE_
        %
        fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,...
            'scl',scl.thd,'grd',grd.thd,...
            'xlm',xlm.thd(end),'xlb',xlb.thd,'xtk',xtk.thd(end),...
            'ylm',ylm.thd(end),'ylb',ylb.thd,'ytk',ytk.thd(end),...
            'leg',leg,'tit',strcat(std,'-',dlg));
        saveas(gcf,strcat(fnm,sprintf('_thd_%s',cpp)),'epsc');
    end
    close all;
    %
    return
end
