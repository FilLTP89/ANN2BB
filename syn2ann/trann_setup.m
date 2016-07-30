%% *SET-UP*
ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
ann.wd = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
    'EMILIA_2905','training');
%ann.wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%    'EMILIA_2905','training');
fprintf('Workdir: %s\n',ann.wd);

ann.nr = 1;
% gh <---> AB
ann.mtd(1).TnC = 0.75;
ann.mtd(1).cp = 'gh';
ann.mtd(1).cl = 'ALL';
% % ud <---> AB
% ann.mtd(2).TnC = 0.75;
% ann.mtd(2).cp = 'ud';
% ann.mtd(2).cl = 'AB';
% % gh <---> CD
% ann.mtd(3).TnC = 0.75;
% ann.mtd(3).cp = 'gh';
% ann.mtd(3).cl = 'CD';
% % ud <---> CD
% ann.mtd(4).TnC = 0.75;
% ann.mtd(4).cp = 'ud';
% ann.mtd(4).cl = 'CD';
% % gh <---> ALL
% ann.mtd(5).TnC = 0.75;
% ann.mtd(5).cp = 'gh';
% ann.mtd(5).cl = 'ALL';
% % ud <---> ALL
% ann.mtd(6).TnC = 0.75;
% ann.mtd(6).cp = 'ud';
% ann.mtd(6).cl = 'ALL';

dbn = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
    'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
%dbn = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%    'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
fprintf('Training Database: %s\n',dbn);

for i_ = 1:ann.nr
    ann.mtd(i_).dbn = dbn;
end
