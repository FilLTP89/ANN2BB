%% *SA SPECTRA*
if flags(1)==1
    for i_ = 1:spg(1)
        xpl{i_} = xplot{i_}.mon.vTn;
        ypl{i_} = abs(xplot{i_}.syn{identity}.psa.(cpp))*utd.psa;
    end
    title = std{1};
elseif flags(1)==2
    for i_ = 1:spg(1)
        disp('PLOTTING GEOEMTRIC MEAN PSA!')
        xpl{i_} = xplot{i_}.mon.vTn;
        ypl{i_} = abs(xplot{i_}.syn{identity}.psa.gh)*utd.psa;
    end
    title = std{1};
    title(end-1:end) = 'GM';
    title=strcat(title,'H');
end
leg = {legplot};
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
        'scl',scl.psa,'grd',grd.psa,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',{xlm.psa},'xlb',xlb.psa,'xtk',{xtk.psa},...
        'ylm',{ylm.psa},'ylb',ylb.psa,'ytk',{ytk.psa},...
        'leg',leg,'tit',{title});
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
        'scl',scl.psa,'grd',grd.psa,...
        'xlm',{xlm.psa},'xlb',xlb.psa,'xtk',{xtk.psa},...
        'ylm',{ylm.psa},'ylb',ylb.psa,'ytk',{ytk.psa},...
        'leg',leg,'tit',{title});
end
if flags(1)==1
    saveas(gcf,strcat(fnm,sprintf('_psa_%s',cpp)),'epsc');
else flags(1)==2
    disp('HOTFIX: PLOT TWICE THE PSA FOR GMH');
    disp('HOTFIX: PLOT TWICE THE PSA FOR GMH');
    disp('HOTFIX: PLOT TWICE THE PSA FOR GMH');
    saveas(gcf,strcat(fnm,sprintf('_psa_gh')),'epsc'); 
end
