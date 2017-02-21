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

%% *MONITOR STATION: monn*
% _station identity_
monn.id = (1:17974)';
job_nb = 2;
sel_id = 1:5;
%monn.id = (1:17974)';
% monn.id = [
%     18446;   
%     16928;
%     18437;
%     17670;
%     15045;
%     16909;
%     16532;
%     13249;
%     13781;
%     13916;
%     16236;
%     17724;
%     16294;
%     12052;
%     10265;
%     13497;
%     12306;
%     17869;
%     14652;
%     14120;
%     13346;
%     14521;
%      3298;
%     12720;  
%     15196;
%      6833;
%     15406;
%     14506;
%      5576;
%     14279;
%      3695;
%      4506;
%      5782];

for i_=1:numel(monn.id)
    bhrr.st{i_}.id = {'';''};
    bhrr.st{i_}.ni = {'';''};
    bhrr.st{i_}.ev = {'20120529.070002'};
    bhrr.st{i_}.dv = {''};
    bhrr.st{i_}.dv = {'itaca'};
end

% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

% _monn field names_
fni.monn = fieldnames(monn);
fnn.monn = numel(fni.monn);

%% *EMPIRICAL ANALYSIS: mtdd*
% _Sabetta&Pugliese 1996 - metadata_
% STORED IN COLUMNS
mtdd.sp96.mw = 6*ones(numel(monn.id),1);
mtdd.sp96.dtm_sp96 = 0.01*ones(numel(monn.id),1);
mtdd.sp96.scc = 2*ones(numel(monn.id),1);
mtdd.sp96.sst = zeros(numel(monn.id),1);
mtdd.sp96.scl = 0.01*ones(numel(monn.id),1);
% _hybridization frequencies_
% horizontal components (stored in (na x 2) matrix)
mtdd.sp96.ew = 1.5*ones(numel(monn.id),2); % in Hz
mtdd.sp96.ns = 1.5*ones(numel(monn.id),2); % in Hz
mtdd.sp96.ud = 1.5*ones(numel(monn.id),2); % in Hz
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
