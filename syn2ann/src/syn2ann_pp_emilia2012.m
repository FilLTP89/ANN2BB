ccc;
delta=importdata('/home/filippo/Data/Filippo/aeolus/ann_new/Summary_Simulated_PG_Values.csv');
datum = load('/home/filippo/Data/Filippo/aeolus/ann_new/syn2ann_res_emilia2012_ALL.mat');
stk = 95*pi/180;

fid_pga=fopen('/home/filippo/Data/Filippo/aeolus/ann_new/pga_emilia2012_network.csv','w+');
fid_pgv=fopen('/home/filippo/Data/Filippo/aeolus/ann_new/pgv_emilia2012_network.csv','w+');
fprintf(fid_pga,'%10s,%20s,%20s,%20s,%20s,%20s,%20s,%20s,\n',...
    'ID','PGA-EW (m/s2)','PGA-NS (cm/s2)','PGA-FP (cm/s2)','PGA-FN (cm/s2)',...
    'PGA-GH (cm/s2)','PGA-LEA-2016 (cm/s2)','PGA-UD (cm/s2),');
fprintf(fid_pgv,'%10s,%20s,%20s,%20s,%20s,%20s,%20s,%20s,\n',...
    'ID','PGV-EW (cm/s)','PGV-NS (cm/s)','PGV-FP (cm/s)','PGV-FN (cm/s)',...
    'PGV-GH (cm/s)','PGV-LEA-2016 (cm/s)','PGV-UD (cm/s),');

for i_=1:datum.mon.na
    %
    sts{i_,1} = datum.bhr.nm{i_};
    [~,~,idx(i_,1)] = intersect(sts{i_,1},delta.textdata(2:end,:));
end

%% COMPARE WITH GMPES
rhypo = delta.data(idx,:);
lons  = str2double(delta.textdata(1+idx,3));
lats  = str2double(delta.textdata(1+idx,2));
idxA = find(cell2mat(regexpi(delta.textdata(1+idx,4),'A'))==1);
idxB = find(cell2mat(regexpi(delta.textdata(1+idx,4),'B'))==1);
idxC = find(cell2mat(regexpi(delta.textdata(1+idx,4),'C'))==1);
ec8c = repmat('D',[datum.mon.na,1]);
ec8c(idxA) = 'A';
ec8c(idxB) = 'B';
ec8c(idxC) = 'C';
gmpe = generate_gmpe('mod','lanzano_et_al_2016_rfocal',...
    'Mw',repmat(6,datum.mon.na,1),'Delta',rhypo,...
    'lons',lons,'lats',lats,'fsof','TF','EC8',ec8c,'bas',ones(datum.mon.na,1));
gmpe.log10pga = 10.^(gmpe.log10pga);
gmpe.log10pgv = 10.^(gmpe.log10pgv);

rjb = 1:50;
njb = numel(rjb);
gmpe_rjb = generate_gmpe('mod','lanzano_et_al_2016_rjb',...
    'Mw',repmat(6,njb,1),'Delta',rjb,...
    'lons',repmat(lons(1),[njb,1]),'lats',repmat(lats(1),[njb,1]),...
    'fsof','TF','EC8',repmat('C',[njb,1]),...
    'bas',ones(njb,1));
save('lanzano_2016.mat','rjb','gmpe_rjb');
return
for i_=1:datum.mon.na
    %
    dtm    = datum.spm.sps.ew.mon.dtm(i_);
    %
    thv_ew = datum.spm.sps.ew.syn{i_}.thv.ew*100;
    thv_ns = datum.spm.sps.ns.syn{i_}.thv.ns*100;
    tha_ew = datum.spm.sps.ew.syn{i_}.tha.ew*100;
    tha_ns = datum.spm.sps.ns.syn{i_}.tha.ns*100;
    %
    [thv_fp,thv_fn]=rotationComponent2D(thv_ew,thv_ns,stk);
    [tha_fp,tha_fn]=rotationComponent2D(tha_ew,tha_ns,stk);
    %
    pgv_ew = datum.spm.sps.ew.syn{i_}.pgv.ew(2)*100;
    pgv_ns = datum.spm.sps.ns.syn{i_}.pgv.ns(2)*100;
    pgv_ud = datum.spm.sps.ud.syn{i_}.pgv.ud(2)*100;
    pga_ew = datum.spm.sps.ew.syn{i_}.pga.ew(2)*100;
    pga_ns = datum.spm.sps.ns.syn{i_}.pga.ns(2)*100;
    pga_ud = datum.spm.sps.ud.syn{i_}.pga.ud(2)*100;
    %
    [~,pga_fp,~,pgv_fp,~,~] = PGAVD_eval(dtm,tha_fp,thv_fp);
    [~,pga_fn,~,pgv_fn,~,~] = PGAVD_eval(dtm,tha_fn,thv_fn);
    %
    fprintf(fid_pga,'%10s,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,\n',...
        sts{i_,1},abs(pga_ew),abs(pga_ns),abs(pga_fp),abs(pga_fn),...
        sqrt(abs(pga_ew*pga_ns)),gmpe.log10pga(i_),abs(pga_ud));
    fprintf(fid_pgv,'%10s,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,%20.5f,\n',...
        sts{i_,1},abs(pgv_ew),abs(pgv_ns),abs(pgv_fp),abs(pgv_fn),...
        sqrt(abs(pgv_ew*pgv_ns)),gmpe.log10pgv(i_),abs(pgv_ud));
end
fclose all;