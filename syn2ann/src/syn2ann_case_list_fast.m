%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_case_list_: function to load all the metadata required to the 
% analyses to be run.
%% *N.B.*
% Need for:

%% *RECORDING STATION: bhrr*
% _station identity_
% MIR01
bhrr.st{1}.id = {'MIR'};
bhrr.st{1}.ni = {'TV';'HN'};
% MRN
bhrr.st{2}.id = {'MRN'};
bhrr.st{2}.ni = {'CIT';'HN'};
% MIR02
bhrr.st{3}.id = {'MIR'};
bhrr.st{3}.ni = {'TV';'HN'};
% SAN0
bhrr.st{4}.id = {'SAN'};
bhrr.st{4}.ni = {'TV';'HN'};
% MIR08
bhrr.st{5}.id = {'MIR'};
bhrr.st{5}.ni = {'TV';'HN'};
% T0802
bhrr.st{6}.id = {'T08'};
bhrr.st{6}.ni = {'TV';'HN'};
% T0813
bhrr.st{7}.id = {'T08'};
bhrr.st{7}.ni = {'TV';'HN'};
% MIR03
bhrr.st{8}.id = {'MIR'};
bhrr.st{8}.ni = {'TV';'HN'};
% T0818
bhrr.st{9}.id = {'T08'};
bhrr.st{9}.ni = {'TV';'HN'};
% MIR04
bhrr.st{10}.id = {'MIR'};
bhrr.st{10}.ni = {'TV';'HN'};
% T0814
bhrr.st{11}.id = {'T08'};
bhrr.st{11}.ni = {'TV';'HN'};
% T0800
bhrr.st{12}.id = {'T08'};
bhrr.st{12}.ni = {'TV';'HN'};
% T0811
bhrr.st{13}.id = {'T08'};
bhrr.st{13}.ni = {'TV';'HN'};
% T0812
bhrr.st{14}.id = {'T08'};
bhrr.st{14}.ni = {'TV';'HN'};
% MIR05
bhrr.st{15}.id = {'MIR'};
bhrr.st{15}.ni = {'TV';'HN'};
% SMS0
bhrr.st{16}.id = {'SMS'};
bhrr.st{16}.ni = {'TV';'HN'};
% RAV0
bhrr.st{17}.id = {'RAV'};
bhrr.st{17}.ni = {'TV';'HN'};
% FIN0
bhrr.st{18}.id = {'FIN'};
bhrr.st{18}.ni = {'TV';'HN'};
% T0824
bhrr.st{19}.id = {'T0824'};
bhrr.st{19}.ni = {'TV';'HN'};
% MOG0
bhrr.st{20}.id = {'MOG'};
bhrr.st{20}.ni = {'TV';'HN'};
% CRP
bhrr.st{21}.id = {'CRP'};
bhrr.st{21}.ni = {'TV';'HN'};
% T0805
bhrr.st{22}.id = {'T08'};
bhrr.st{22}.ni = {'TV';'HN'};
% MIR06
bhrr.st{23}.id = {'MIR'};
bhrr.st{23}.ni = {'TV';'HN'};
% CNT
bhrr.st{24}.id = {'CNT'};
bhrr.st{24}.ni = {'TV';'HN'};
% T0803
bhrr.st{25}.id = {'T08'};
bhrr.st{25}.ni = {'TV';'HN'};
% SERM
bhrr.st{26}.id = {'SERM'};
bhrr.st{26}.ni = {'TV';'HN'};
% SAG0
bhrr.st{27}.id = {'SAG'};
bhrr.st{27}.ni = {'TV';'HN'};
% CAS05
bhrr.st{28}.id = {'CAS'};
bhrr.st{28}.ni = {'TV';'HN'};
% CAS0
bhrr.st{29}.id = {'CAS'};
bhrr.st{29}.ni = {'TV';'HN'};
% BON0
bhrr.st{30}.id = {'BON'};
bhrr.st{30}.ni = {'TV';'HN'};
% MODE
bhrr.st{31}.id = {'MODE'};
bhrr.st{31}.ni = {'TV';'HN'};
% MDN
bhrr.st{32}.id = {'MDN'};
bhrr.st{32}.ni = {'TV';'HN'};
% NVL
bhrr.st{33}.id = {'NVL'};
bhrr.st{33}.ni = {'TV';'HN'};
% AQK
bhrr.st{34}.id = {'AQK'};
bhrr.st{34}.ni = {'IT';'HN'};
% AQU
bhrr.st{35}.id = {'AQU'};
bhrr.st{35}.ni = {'MN';'HL'};

% _recorded events_
bhrr.st{1}.ev  = {'20120529.070002'};
bhrr.st{2}.ev  = {'20120529.070002'};
bhrr.st{3}.ev  = {'20120529.070002'};
bhrr.st{4}.ev  = {'20120529.070002'};
bhrr.st{5}.ev  = {'20120529.070002'};
bhrr.st{6}.ev  = {'20120529.070002'};
bhrr.st{7}.ev  = {'20120529.070002'};
bhrr.st{8}.ev  = {'20120529.070002'};
bhrr.st{9}.ev  = {'20120529.070002'};
bhrr.st{10}.ev = {'20120529.070002'};
bhrr.st{11}.ev = {'20120529.070002'};
bhrr.st{12}.ev = {'20120529.070002'};
bhrr.st{13}.ev = {'20120529.070002'};
bhrr.st{14}.ev = {'20120529.070002'};
bhrr.st{15}.ev = {'20120529.070002'};
bhrr.st{16}.ev = {'20120529.070002'};
bhrr.st{17}.ev = {'20120529.070002'};
bhrr.st{18}.ev = {'20120529.070002'};
bhrr.st{19}.ev = {'20120529.070002'};
bhrr.st{20}.ev = {'20120529.070002'};
bhrr.st{21}.ev = {'20120529.070002'};
bhrr.st{22}.ev = {'20120529.070002'};
bhrr.st{23}.ev = {'20120529.070002'};
bhrr.st{24}.ev = {'20120529.070002'};
bhrr.st{25}.ev = {'20120529.070002'};
bhrr.st{26}.ev = {'20120529.070002'};
bhrr.st{27}.ev = {'20120529.070002'};
bhrr.st{28}.ev = {'20120529.070002'};
bhrr.st{29}.ev = {'20120529.070002'};
bhrr.st{30}.ev = {'20120529.070002'};
bhrr.st{31}.ev = {'20120529.070002'};
bhrr.st{32}.ev = {'20120529.070002'};
bhrr.st{33}.ev = {'20120529.070002'};
bhrr.st{34}.ev = {'20090406.013240'};
bhrr.st{35}.ev = {'20090406.013240'};
% _device list_
bhrr.st{1}.dv = {'01'};
bhrr.st{2}.dv = {''};
bhrr.st{3}.dv = {'02'};
bhrr.st{4}.dv = {'0'};
bhrr.st{5}.dv = {'08'};
bhrr.st{6}.dv = {'02'};
bhrr.st{7}.dv = {'13'};
bhrr.st{8}.dv = {'03'};
bhrr.st{9}.dv = {'18'};
bhrr.st{10}.dv = {'04'};
bhrr.st{11}.dv = {'14'};
bhrr.st{12}.dv = {'00'};
bhrr.st{13}.dv = {'11'};
bhrr.st{14}.dv = {'12'};
bhrr.st{15}.dv = {'05'};
bhrr.st{16}.dv = {'0'};
bhrr.st{17}.dv = {'0'};
bhrr.st{18}.dv = {'0'};
bhrr.st{19}.dv = {'24'};
bhrr.st{20}.dv = {'0'};
bhrr.st{21}.dv = {''};
bhrr.st{22}.dv = {'05'};
bhrr.st{23}.dv = {'06'};
bhrr.st{24}.dv = {''};
bhrr.st{25}.dv = {'03'};
bhrr.st{26}.dv = {''};
bhrr.st{27}.dv = {'0'};
bhrr.st{28}.dv = {'05'};
bhrr.st{29}.dv = {'0'};
bhrr.st{30}.dv = {'0'};
bhrr.st{31}.dv = {''};
bhrr.st{32}.dv = {''};
bhrr.st{33}.dv = {''};
bhrr.st{34}.dv = {''};
bhrr.st{35}.dv = {''};

% _database_
bhrr.st{1}.tp =  {'itaca'};
bhrr.st{2}.tp =  {'itaca'};
bhrr.st{3}.tp =  {'itaca'};
bhrr.st{4}.tp =  {'itaca'};
bhrr.st{5}.tp =  {'itaca'};
bhrr.st{6}.tp =  {'itaca'};
bhrr.st{7}.tp =  {'itaca'};
bhrr.st{8}.tp =  {'itaca'};
bhrr.st{9}.tp =  {'itaca'};
bhrr.st{10}.tp = {'itaca'};
bhrr.st{11}.tp = {'itaca'};
bhrr.st{12}.tp = {'itaca'};
bhrr.st{13}.tp = {'itaca'};
bhrr.st{14}.tp = {'itaca'};
bhrr.st{15}.tp = {'itaca'};
bhrr.st{16}.tp = {'itaca'};
bhrr.st{17}.tp = {'itaca'};
bhrr.st{18}.tp = {'itaca'};
bhrr.st{19}.tp = {'itaca'};
bhrr.st{20}.tp = {'itaca'};
bhrr.st{21}.tp = {'itaca'};
bhrr.st{22}.tp = {'itaca'};
bhrr.st{23}.tp = {'itaca'};
bhrr.st{24}.tp = {'itaca'};
bhrr.st{25}.tp = {'itaca'};
bhrr.st{26}.tp = {'itaca'};
bhrr.st{27}.tp = {'itaca'};
bhrr.st{28}.tp = {'itaca'};
bhrr.st{29}.tp = {'itaca'};
bhrr.st{30}.tp = {'itaca'};
bhrr.st{31}.tp = {'itaca'};
bhrr.st{32}.tp = {'itaca'};
bhrr.st{33}.tp = {'itaca'};
bhrr.st{34}.tp = {'itaca'};
bhrr.st{35}.tp = {'itaca'};

% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *MONITOR STATION: monn*
% _station identity_
% monn.id = [16928,15045,1,2,18446,18437,18446
monn.id = [...
    18446;
    16928; 
    18437;
    17670;
    15045;
    16909;
    16532;
    13249;
    13781;
    13916;
    16236;
    17724;
    16294;
    12052;
    10265;
    13497;
    12306;
    17869;
    14652;
    14120;
    13346;
    14521;
     3298;
    12720;
    15196;
     6833;
    15406;
    14506;
     5576;
    14279;
     3695;
     4506;
     5782;
        1;
        2];

% _monn field names_
fni.monn = fieldnames(monn);
fnn.monn = numel(fni.monn);

%% *EMPIRICAL ANALYSIS: mtdd*
% _Sabetta&Pugliese 1996 - metadata_
% STORED IN COLUMNS
mtdd.sp96.mw = [6.0*ones(33,1);6.3;6.3];
mtdd.sp96.dtm_sp96 = 0.01*ones(35,1);
% site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m));
mtdd.sp96.scc = [2*ones(33,1);1;2];
mtdd.sp96.sst = zeros(35,1);
mtdd.sp96.scl = 0.01*ones(35,1);
% _hybridization frequencies_
% horizontal components (stored in (na x 2) matrix)
%mtdd.sp96.ew = [1.5,1.5,1.5,1.5,1.5,1.5;...
%    1.5,1.5,1.5,1.5,1.5,1.5]'; % in Hz
%mtdd.sp96.ns = [1.5,1.5,1.5,1.5,1.5,1.5;...
%    1.5,1.5,1.5,1.5,1.5,1.5]'; % in Hz
mtdd.sp96.ew = 1.5*ones(35,2);
mtdd.sp96.ns = 1.5*ones(35,2);
% vertical components (stored in (na x 2) matrix)  
%mtdd.sp96.ud = [1.5,1.5,1.5,1.5,1.5,1.5;...
%    1.5,1.5,1.5,1.5,1.5,1.5]'; % in Hz
mtdd.sp96.ud = 1.5*ones(35,2);
%% _Exsim - reference files_
%mtdd.exsim.fnm{1} = fullfile(wd,'exsim_old','exsim_emilia','MRN_new');
%mtdd.exsim.fnm{2} = fullfile(wd,'exsim_old','exsim_emilia','MIR08_new');
%mtdd.exsim.fnm{3} = fullfile(wd,'exsim_old','exsim_aquila','walters','aqk');
%mtdd.exsim.fnm{4} = fullfile(wd,'exsim_old','exsim_aquila','walters','aqu');
%mtdd.exsim.pf{1} = strcat(mtdd.exsim.fnm{1},filesep,'MRN_exsim_');
%mtdd.exsim.pf{2} = strcat(mtdd.exsim.fnm{2},filesep,'MIR08_exsim_');
%mtdd.exsim.pf{3} = strcat(mtdd.exsim.fnm{3},filesep,'aqk_0604_');
%mtdd.exsim.pf{4} = strcat(mtdd.exsim.fnm{4},filesep,'aqu_0604_');

% _mtdd field names_
fni.mtdd.sp96 = fieldnames(mtdd.sp96);
fnn.mtdd.sp96 = numel(fni.mtdd.sp96);
%fni.mtdd.exsim = fieldnames(mtdd.exsim);
%fnn.mtdd.exsim = numel(fni.mtdd.exsim);
