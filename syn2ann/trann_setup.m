%% *SET-UP*
ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
% ann.wd = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
%     'EMILIA_2905','training');
ann.wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
    'EMILIA_2905','training');
fprintf('Workdir: %s\n',ann.wd);
ann.nr = 2;
ann.mtd(1).tc = 0.75;
ann.mtd(1).cp = 'gh';
ann.mtd(2).tc = 0.75;
ann.mtd(2).cp = 'ud';
ann.mtd(3).tc = 0.75;
ann.mtd(3).cp = 'gh';

% dbn = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
%     'EMILIA_2905','simbad_v04','SIMBAD_v04.mat');
dbn = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
    'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
ftr = fullfile(ann.wd,'index_train.txt');
fva = fullfile(ann.wd,'index_val.txt');
fprintf('Training Database: %s\n',dbn);
fprintf('Training Values: %s\n',ftr);
fprintf('Validation Database: %s\n',fva);

for i_ = 1:ann.nr
    ann.mtd(i_).dbn = dbn;
    ann.mtd(i_).ftr = ftr;
    ann.mtd(i_).fva = fva;
end