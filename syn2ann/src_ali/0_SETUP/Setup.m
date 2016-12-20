%% *SET-UP*
fprintf('---------------------\n0. SETUP\n---------------------\n');

%% *WORKDIR*
%wd = fullfile(filesep,'Users','maria','Documents','PHD','Work','Generation_BB_Groundmotion_with_ANN');
wd='/media/filippo/Data/Filippo/PHD_passing_through_polimi/syn2ann';
wd_data =  fullfile(wd,'database_earthqks','EMILIA_120529');
wd_results = fullfile(wd,'database_earthqks','EMILIA_120529','PPT_SPEED_results');
mkdir(wd_results);

%% *DEFINE INPUT PARAMETERS*
%MRN
Mw = 6;
lon = 11.062;
lat = 44.878;
R_epi = 4.1; % epicentral distance [Km]
scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
% %AQK
% Mw = 6.3;
% lon = 13.4009;
% lat = 42.3450;
% R_epi = 5.7; % epicentral distance [Km]
% scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
% %AQU
% Mw = 6.3;
% lon = 13.4019;
% lat = 42.3539;
% R_epi = 6; % epicentral distance [Km]
% scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
% %Airport Istanbul
% Mw = 7;
% lon = 28.81058;
% lat = 40.98558;
% R_epi = 12; % epicentral distance [Km] !!!CORRECT IS THE APPROXIMATELY Rjb
% scc = 2; % site conditions (0=rock, (Vs,30>800m/s); 1=shallow all. (H<=20); 2=deep alluvium (H>20m)); 
%% *NUMERICAL SIMULATION*
% simulationID
% %MRN
num_sim.simID = 'E00001'; 
num_sim.monID = '16928';
% %AQK
% num_sim.simID = 'E00001'; 
% num_sim.monID = '1';
%AQU
% num_sim.simID = 'E00001'; 
% num_sim.monID = '2';
% %Airport Istanbul
% num_sim.simID = 'E00508'; 
% num_sim.monID = '4560';
%% *REAL RECORDS*
cfr_record = 1; % decide if compare the simulations (0) with real records or not (1) 
if cfr_record == 1
    record.station = 'MRN';
%    record.station = 'AQK';
%   record.station = 'AQU';
   record.lon = lon;
   record.lat = lat;
end 

%% *ANN NETWORK*
%% **net_075** 
load('net_75_ALL_gh_new.mat');net_h=net;
load('net_75_ALL_ud_new.mat');net_v=net;
inp_Tn = [0.75,0.8:0.1:1.0,1.25:0.25:5.0];
tar_Tn = [0,0.05,0.1:0.1:0.7];
