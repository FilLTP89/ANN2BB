function trann_test_psa_performance(varargin)
    %% *SET-UP*
    rec = varargin{1};
    ann = varargin{2};
    identity = varargin{3};
    stdd = varargin{4};
    fnm = varargin{5};
    col = [0,0,0;0.4,0.4,0.4;0,0,0];
    set(0,'defaultaxescolororder',col);
    mrk = {'none','s','v'};
    lwd = 1;
    plot_set_up;
    xlm = [0;1];
    ylm = [-1;1];
    [xtk,ytk] = get_axis_tick(xlm,ylm,abs(diff(xlm))/4,abs(diff(ylm))/4);
    
    %% *COMPUTE TRENDS*
    xpl = cell(numel(ann),1);
    ypl = cell(numel(ann),1);
    leg = cell(numel(ann),1);
    
    
    for i_=1:numel(ann)
        cpn = seismo_dir_conversion(ann{i_}.cpp);
        [inp.idx,tar.idx] = trann_check_vTn(ann{i_}.inp,ann{i_}.tar,rec.mon,1e-8);
        xpl{i_,1} = ann{i_}.tar.vTn(:)./ann{i_}.TnC;
        ypl{i_,1} = ones(numel(tar.idx),1);
        
        for j_=1:numel(cpn)
            fprintf('GEOMETRIC MEAN ON RECORDED %s-PSA\n',upper(cpn{j_}));
            ypl{i_,1} = ypl{i_,1}.*rec.syn{identity}.psa.(cpn{j_})(tar.idx,1);
        end
        if numel(cpn)>1
            disp('---> averaging over horizontal components');
            ypl{i_,1} = ypl{i_,1}.^(1./numel(cpn));
        end
        
        idx = ann{i_}.mon.vTn < ann{i_}.TnC;
        ypl{i_,1} = log10(ypl{i_,1}./...
            ann{i_}.syn{identity}.psa.(ann{i_}.cpp)(idx,1));
        leg{i_,1} = strcat('$',ann{i_}.scl,'\left(',upper(ann{i_}.cpp),...
            ';T^\star=',num2str(ann{i_}.TnC,'%.2f'),'s\right)$');
    end
    
    fpplot('xpl',xpl,'ypl',ypl,'pfg',[0,0,10,10],...
        'mrk',mrk,'lwd',lwd,'tit',{stdd},...
        'xlb',{'T/T* [1]'},'xlm',{xlm},'xtk',{xtk},...
        'ylb',{'log_{10}(Sa_{REC}/Sa_{ANN}) [1]'},'ylm',{ylm},'ytk',{ytk},...
        'leg',{leg});
    
    saveas(gcf,strcat(fnm,'_prf'),'epsc');
    
    return
end
