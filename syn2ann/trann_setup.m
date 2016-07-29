%% *SET-UP*
ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
ann.wd = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
    'EMILIA_2905','training');
%ann.wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%    'EMILIA_2905','training');
fprintf('Workdir: %s\n',ann.wd);

ann.nr = 2;
ann.mtd(1).TnC = 0.75;
ann.mtd(1).cp = 'gh';
ann.mtd(2).TnC = 0.75;
ann.mtd(2).cp = 'ud';

dbn = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
    'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
%dbn = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%    'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
fprintf('Training Database: %s\n',dbn);

for i_ = 1:ann.nr
    ann.mtd(i_).dbn = dbn;
end
