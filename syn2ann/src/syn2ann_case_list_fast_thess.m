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

% _recorded events_
bhrr.st{1}.ev  = {'20120529.070002'};
bhrr.st{2}.ev  = {'20120529.070002'};

% _database_
bhrr.st{1}.tp =  {'itaca'};
bhrr.st{2}.tp =  {'itaca'};

% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *MONITOR STATION: monn*
% _station identity_
% monn.id = [16928,15045,1,2,18446,18437,18446
monn.id = [...
    1406;
    1776];

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
mtdd.sp96.ew = 3.0*ones(35,2);
mtdd.sp96.ns = 3.0*ones(35,2);
% vertical components (stored in (na x 2) matrix)  
%mtdd.sp96.ud = [1.5,1.5,1.5,1.5,1.5,1.5;...
%    1.5,1.5,1.5,1.5,1.5,1.5]'; % in Hz
mtdd.sp96.ud = 3.0*ones(35,2);
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
