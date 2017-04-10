function lfhf_filter_plot(varargin)
    dtm = varargin{1};
    vtm = (0:dtm:varargin{2})';
    ntm = numel(vtm);
    nfr = 2^nextpow2(ntm);
    dfr = 1/dtm/nfr;
    vfr = dfr*(1:nfr)';
    fNy = 0.5/dtm;
    fhp = varargin{3};
    flp = varargin{4};
%     nfa = round(fhp/dfr);
%     nfb = round(flp/dfr);
%     fac = pi./(dfr*(nfb-nfa-1));
%     WLF = zeros(nfr,1);
%     nfa = round(fhp/dfr);
%     nfb = round(flp/dfr);
%     WLF(1:nfa) = 1.0;
%     WLF(nfa+1:nfb) = 0.5*(1+cos(fac*dfr*(0:nfb-nfa-1)));
%     WHF = 1-WLF;
    color_order = [rgb('blue');rgb('red');0,0,0;0,0,0;0,0,0];
    set(0,'defaultaxescolororder',color_order);
%     fpplot('xpl',{vfr;vfr;[fhp;fhp];[flp;flp]},'ypl',{WLF;WHF;[0,1];[0,1]},...
%         'pfg',[0 0 12 12],'scl',{'slx'},'xlb',{'f [Hz]'},'xlm',{[1e-1,5e1]},...
%         'xtk',{[fhp,flp]},'ylb',{'[1]'},'ylm',{[-.05,1.3]},...
%         'ytk',{0:.25:1},'leg',{{'$w^{LF}\left(f\right)$','$w^{HF}\left(f\right)$'}},...
%         'tit',{'LF-HF HYBRIDIZATION'},'lst',{'-';'-';'--';'--'},'vfg','on');
%     set(gca,'xticklabel',{'$f_{LF}$';'$f_{HF}$'});
%     set(gca,'TickLabelInterpreter', 'latex');
%     
%     
%     color_order = [rgb('Blue');rgb('Red');rgb('IntenseGreen');rgb('IntenseGreen')];
%     set(0,'defaultaxescolororder',color_order);
%     fpplot('xpl',{vfr;vfr;[fhp;fhp];[flp;flp]},'ypl',{WLF;WHF;[0,1];[0,1]},...
%         'pfg',[0 0 12 12],'scl',{'slx'},'xlb',{'f [Hz]'},'xlm',{[1e-1,5e1]},...
%         'xtk',{[fhp,flp]},'ylb',{'[1]'},'ylm',{[-.05,1.3]},...
%         'ytk',{0:.25:1},'leg',{{'$w^{LF}\left(f\right)$','$w^{HF}\left(f\right)$'}},...
%         'tit',{'LF-HF HYBRIDIZATION'},'lst',{'-';'-';'--';'--'},'vfg','on');
%     set(gca,'xticklabel',{'$f_{LF}$';'$f_{HF}$'});
%     set(gca,'TickLabelInterpreter', 'latex');
    bfo = 3;
    [bwb.hp,bwa.hp,~] = create_butter_filter(bfo,fhp,[],fNy);
    [bwb.lp,bwa.lp,~] = create_butter_filter(bfo,[],flp,fNy);
    omg = linspace(0, pi, nfr);
    ker = exp(-1i.*omg(:));
    cfc = 0:bfo;
    nom.hp   = -999*ones(nfr,bfo+1);
    denom.hp = -999*ones(nfr,bfo+1);
    nom.lp   = -999*ones(nfr,bfo+1);
    denom.lp = -999*ones(nfr,bfo+1);
    for i_=1:numel(cfc)
        %
        nom.hp(:,i_)   = bwb.hp(i_).*(ker(:).^cfc(i_));
        denom.hp(:,i_) = bwa.hp(i_).*(ker(:).^cfc(i_));
        %
        nom.lp(:,i_)   = bwb.lp(i_).*(ker(:).^cfc(i_));
        denom.lp(:,i_) = bwa.lp(i_).*(ker(:).^cfc(i_));
    end
    
    bwf.hp = sum(nom.hp,2)./sum(denom.hp,2);
    bwf.lp = sum(nom.lp,2)./sum(denom.lp,2);
    
    close all;
    fpplot('xpl',{vfr;vfr;[flp;flp];2.*[flp;flp];3.*[flp;flp]},...
        'ypl',{abs(bwf.lp);abs(bwf.hp);[0,1];[0,1];[0,1]},...
        'pfg',[0 0 18 13],'scl',{'slx'},'xlb',{'f [Hz]'},'xlm',{[1e0,fNy]},...
        'xtk',{[1e0,fhp,2*fhp,3*fhp,4*fhp,5*fhp,fNy]},'ylb',{'[1]'},'ylm',{[-.05,1.02]},...
        'ytk',{0:.25:1},'leg',{{'\boldmath$w^{LF}\left(f\right)$',...
        '\boldmath$w^{HF}\left(f\right)$'}},...
        'lst',{'-';'-';'--';'--';'--';'--';'--'},'vfg','on',...
        'lwd',[5,5,2,2,2]);
    set(gca,'xticklabel',{'\boldmath$1$';'\boldmath$f_{C}$';...
        '\boldmath$2f_{C}$';'\boldmath$3f_{C}$';'\boldmath$f_{Nyq}$'});
    set(gca,'TickLabelInterpreter', 'latex');
    obj = get(gcf,'children');
    obj(1).Location='northoutside';
    obj(1).Orientation='horizontal';
    obj(1).FontSize = 17;
    
%     color_order = [rgb('Blue');rgb('Red');rgb('IntenseGreen');rgb('IntenseGreen')];
%     set(0,'defaultaxescolororder',color_order);
%     fpplot('xpl',{vfr;vfr;[fhp;fhp];[flp;flp]},'ypl',{WLF;WHF;[0,1];[0,1]},...
%         'pfg',[0 0 12 12],'scl',{'slx'},'xlb',{'f [Hz]'},'xlm',{[1e-1,5e1]},...
%         'xtk',{[fhp,flp]},'ylb',{'[1]'},'ylm',{[-.05,1.3]},...
%         'ytk',{0:.25:1},'leg',{{'$w^{LF}\left(f\right)$','$w^{HF}\left(f\right)$'}},...
%         'tit',{'LF-HF HYBRIDIZATION'},'lst',{'-';'-';'--';'--'},'vfg','on');
%     set(gca,'xticklabel',{'$f_{LF}$';'$f_{HF}$'});
%     set(gca,'TickLabelInterpreter', 'latex');

    if nargin>4
        sp = varargin{5};
        saveas(gcf,sp,'epsc');
    end
    
    return
end
