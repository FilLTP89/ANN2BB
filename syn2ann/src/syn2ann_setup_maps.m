%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_setup_maps_: function to select the set of analyses to the
% analyses to be run.
%% *N.B.*
% Need for:
% _ccc.m,syn2ann_case_list.m_

ccc;
fprintf('============================\n');
fprintf('----------0. SETUP----------\n');
fprintf('============================\n');

%% *DEFINE WORKDIR*
% _main workdir_
wd = '/media/filippo/Data/Filippo/PHD_passing_through_polimi/syn2ann/database';
% _save path_
sp = '/home/filippo/Scrivania/ann';

if exist(wd,'dir')~=7
    % _main workdir_
    wd = '/workdir/gattif/syn2ann_maps';
    % _save path_
    sp = '/workdir/gattif/syn2ann_maps';
end

%% *LOAD ALL METADATA AVAILABLE*
fprintf('Workdir: %s\n',wd);
syn2ann_case_list_maps;

hybrid_type='sp96';
MAXIT = 5;
%% *DEFINE ANN METADATA*
% ANN vector follow the reference system mon.rs
% _site class considered : ALL,AB,CD_
ann.mtd.scl = {'ALL';'ALL';'ALL'};
% _corner period for each ANN_
ann.mtd.TnC = {0.75;0.75;0.75};
% _ANN motion component : gh,ud (geometric mean horizontal, updip)_
ann.mtd.cpn = {'gh';'gh';'ud'};
for i_ = 1:numel(ann.mtd.TnC)
    ann.mtd.nl(i_) = {sprintf('net_%u_%s_%s_new.mat',...
        round(ann.mtd.TnC{i_}*100),ann.mtd.scl{i_},ann.mtd.cpn{i_})};
    ann.mtd.tol(i_).low.psa = 0.15;
    ann.mtd.tol(i_).low.pga = 0.15;
    ann.mtd.tol(i_).hgh.psa = 0.15;
    ann.mtd.tol(i_).hgh.pga = 0.15;
end
% tolerances
ann.mtd.nit = [10;10;40];
