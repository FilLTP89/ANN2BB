%% *LOAD TRAINED ANN*
fprintf('---------------------\n1. ANN (JUST PSA)\n---------------------\n');
% parse ann networks on each motion component
fprintf('--> Parsing \n');

for j_ = 1:tst.mtd.nr
    fn = fullfile(wd,'training',tst.mtd.nl{j_});
    TnC  = tst.mtd.TnC{j_};
    ann.tst{j_} = syn2ann_ann_parser(TnC,fn,tst.mtd.cpp{j_},tst.mtd.scl{j_});
end
%% *LOAD TEST SET*
fprintf('---------------------\n2. LOAD TEST SET\n---------------------\n');

Tn_min   = 0;
% maximum natural period
Tn_max   = 5;
% natural period step
dTn      = 0.05;
rec.org.mon.vTn  = (Tn_min:dTn:Tn_max)';
rec.org.mon.nT   = numel(rec.org.mon.vTn);
rec.org.mon.zeta = 0.05;
rec.org.mon.nc = 2;
rec.org.mon.cp = {'ew','ns'};
rec.org.mon.na = numel(bhr.sns);

for i_=1:rec.org.mon.na
    [~,rec.org.idx_sns,~]=intersect(rec.org.mon.vTn,bhr.sns{i_}.vTn);
    rec.org.syn{i_}.psa.ew=zeros(rec.org.mon.nT,1);
    rec.org.syn{i_}.psa.ns=zeros(rec.org.mon.nT,1);
    rec.org.syn{i_}.psa.ew(rec.org.idx_sns,1)=bhr.sns{i_}.psa;
    rec.org.syn{i_}.psa.ns(rec.org.idx_sns,1)=bhr.sns{i_}.psa;
end