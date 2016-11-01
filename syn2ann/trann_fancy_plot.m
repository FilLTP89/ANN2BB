%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_plot_res_station_: function to plot ANN test results per station
%% *N.B.*
% Need for:
% _syn2ann_plot_compare.m_

global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});
set(0,'defaultaxescolororder',clr0f)
%% *POST-PROCESS - RECORDS vs ANN*
plt.spg = [2,3];
plt.pax = num2cell([1;2;3;4;4;5;5;6;6]);
plt.mrk = {'none';'none';'none';'none';'o';'none';'o';'none';'o'};
plt.lst = {'-';'-';'-';'-';':';'-';':';'-';':'};

plt.xpl = cell(9,1);
plt.xlm = cell(6,1);
plt.xtk = cell(6,1);
plt.xlb = cell(6,1);
plt.ypl = cell(9,1);
plt.ylm = cell(6,1);
plt.ytk = cell(6,1);
plt.ylb = cell(6,1);
plt.scl = cell(6,1);
plt.grd = cell(6,1);
plt.leg = cell(6,1);
plt.tit = cell(6,1);
for j_ = 1:numel(cpp.rec)
    flag.rec = seismo_dir_conversion(cpp.rec{j_});
    for kk_ = 1:tst.mtd.nr
        cpp.ann = tst.mtd.cpp{kk_};
        flag.ann = seismo_dir_conversion(cpp.ann);
        if any(strcmpi(flag.rec,flag.ann))
            [~,tha_lab] = get_axis_tick(tha_lim,tha_lim,1,200);
            % x-plot
            plt.xpl{j_,1}   = rec.mon.vtm{mm_};
            plt.xpl{3+2*j_-1,1} = rec.mon.vTn;
            plt.xpl{3+2*j_,1}   = ann.tst{kk_}.mon.vTn;
            plt.xlm(j_,1)   = {vtm_lim};
            plt.xlm(j_+3,1) = xlm.psa(1);
            plt.xtk(j_,1)   = {vtm_lab};
            plt.xtk(j_+3,1) = xtk.psa(1);
            plt.xlb(j_,1)   = {'t [s]'};
            plt.xlb(j_+3,1) = xlb.psa(1);
            % y-plot
            plt.ypl{j_,1}   = rec.syn{mm_}.tha.(cpp.rec{j_})*100;
            plt.ypl{3+2*j_-1,1} = rec.syn{mm_}.psa.(cpp.rec{j_})*100;
            plt.ypl{3+2*j_,1}   = ann.tst{kk_}.syn{mm_}.psa.(cpp.rec{j_})*100;
            plt.ylm(j_,1)   = {tha_lim};
            plt.ylm(j_+3,1) = ylm.psa(1);
            plt.ytk(j_,1)   = {tha_lab};
            plt.ytk(j_+3,1) = ytk.psa(1);
            plt.ylb(j_,1)   = {'a(t) [cm/s/s]'};
            plt.ylb(j_+3,1) = ylb.psa(1);
            % other
            plt.scl(j_,1)   = {'lin'};
            plt.scl(j_+3,1) = scl.psa(j_);
            plt.grd(j_,1)   = {'minor'};
            plt.grd(j_+3,1) = grd.psa;
            plt.grd(j_+6,1) = grd.psa;
            plt.leg(j_,1)   = {strcat('REC-',dlg(j_))};
            plt.leg(j_+3,1) = {{sprintf('REC-%s',dlg{j_});sprintf('ANN-%s',dlg{j_})}}; 
            
            fnm = fullfile(spp,sprintf('%s_rec_%s_tha_%s',st,tst.mtd.scl{kk_},cpp.rec{j_}));
            fpplot('xpl',plt.xpl(j_,1),'ypl',plt.ypl(j_,1),...
                'pfg',[0 0 14 5.25],...
                'scl',plt.scl(j_,1),'grd',plt.grd(j_,1),...
                'xlm',plt.xlm(j_,1),'xtk',plt.xtk(j_,1),'xlb',plt.xlb(j_,1),...
                'ylm',plt.ylm(j_,1),'ytk',plt.ytk(j_,1),'ylb',plt.ylb(j_,1),...
                'leg',{{'REC'}},'tit',dlg(j_));
%             saveas(gcf,fnm,'jpg');
            saveas(gcf,fnm,'epsc');
            
            fnm = fullfile(spp,sprintf('%s_rec_%s_psa_%s',st,tst.mtd.scl{kk_},cpp.rec{j_}));
            fpplot('xpl',plt.xpl(3+2*j_-1:3+2*j_,1),'ypl',plt.ypl(3+2*j_-1:3+2*j_,1),...
                'pfg',[0 0 14 14],...
                'scl',plt.scl(j_+3,1),'grd',plt.grd(j_+3,1),...
                'xlm',plt.xlm(j_+3,1),'xtk',plt.xtk(j_+3,1),'xlb',plt.xlb(j_+3,1),...
                'ylm',plt.ylm(j_+3,1),'ytk',plt.ytk(j_+3,1),'ylb',plt.ylb(j_+3,1),...
                'leg',{{'REC';'ANN'}},'tit',dlg(j_),'mrk',{'none';'o'},'lst',{'-';':'});
%             saveas(gcf,fnm,'jpg');
            saveas(gcf,fnm,'epsc');
        end
    end
end

fnm = fullfile(spp,sprintf('%s_rat_%s',st,tst.mtd.scl{kk_}));
fpplot('xpl',plt.xpl,'ypl',plt.ypl,...
    'pfg',pfg.fat,'spg',plt.spg,'pax',plt.pax,...
    'scl',plt.scl,'grd',plt.grd,...
    'xlm',plt.xlm,'xtk',plt.xtk,'xlb',plt.xlb,...
    'ylm',plt.ylm,'ytk',plt.ytk,'ylb',plt.ylb,...
    'leg',plt.leg,'mrk',plt.mrk,'lst',plt.lst);
% saveas(gcf,fnm,'jpg');
saveas(gcf,fnm,'epsc');
