function lfhf_filter_plot(varargin)
    dtm = varargin{1};
    vtm = (0:dtm:varargin{2})';
    ntm = numel(vtm);
    nfr = 2^nextpow2(ntm);
    dfr = 1/dtm/(nfr-1);
    vfr = dfr*(0:nfr-1)';
    fr_max = 0.5/dtm;
    fa = varargin{3};
    fb = varargin{4};
    nfa = round(fa/dfr);
    nfb = round(fb/dfr);
    fac = pi./(dfr*(nfb-nfa-1));
    WLF = zeros(nfr,1);
    nfa = round(fa/dfr);
    nfb = round(fb/dfr);
    WLF(1:nfa) = 1.0;
    WLF(nfa+1:nfb) = 0.5*(1+cos(fac*dfr*(0:nfb-nfa-1)));
    WHF = 1-WLF;
    fpplot('xpl',{vfr;vfr;[fa;fa];[fb;fb]},'ypl',{WLF;WHF;[0,1];[0,1]},...
        'pfg',[0 0 16 16],'scl',{'slx'},'xlb',{'f [Hz]'},'xlm',{[1e-1,1e1]},...
        'xtk',{[1e-1,1,1e1]},'ylb',{'[1]'},'ylm',{[-.05,1.05]},...
        'ytk',{0:.25:1},'vfg','on','leg',{{'LF-filter','HF-filter'}},...
        'tit',{'LF-HF HYBRIDIZATION'},'lst',{'-';'-';'--';'--'});
    return
end
