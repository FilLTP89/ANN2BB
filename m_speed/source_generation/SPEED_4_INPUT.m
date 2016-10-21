%**************************************************************************
% Copyright (C) 2012 The SPEED FOUNDATION
% Author: SPEED-Team (POLIMI)
%         Politecnico di Milano 
%         P.zza Leonardo da Vinci, 32 
%         20133 Milano 
%         Italy                                        
%
% This file is part of SPEED.
%
% SPEED is free software; you can redistribute it and/or modify it
% under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% SPEED is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with SPEED.  If not, see <http://www.gnu.org/licenses/>.
%************************************************************************** 


function varargout = SPEED_4_INPUT(id)

% MATLAB code for the generation of SPEED inputs 
% Set up the following variables according to your data 

addpath src

% MODIFY THE FOLLOWING PARAMETERS ACCORDING WITH 
% THE SCENARIOS DATABASE

handles.OS = 'Unix';

% SCENARIO ID
handles.SCENARIO_ID = id;

% GRIDFILE & MATEFILE

%% FROM 1 TO 100 --> SANTIAGO
handles.GRID = 'Santiago_SanRamon_official';
handles.MATE = ['Santiago_SanRamon_official_', num2str(id)];
% NOT-HONORING SANTIAGO
handles.DEPTH = -20;
handles.WR = 'R';
handles.NH_MATE = 1;
handles.NH_TOL = 10;
handles.TAG_CASE = 8;

%% FROM 101 TO 200 ---> EMILIA F11
% handles.GRID = 'EmiliaHisada_2905_Dip60';
% handles.MATE = ['Emilia_2905_', num2str(id)];
% % NOT-HONORING EMILIA
% handles.DEPTH = -11000;
% handles.WR = 'W';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10';
% handles.TAG_CASE = 12;

%% FROM 101 TO 200 ---> EMILIA F12
% handles.GRID = 'Emilia_F12';
% handles.MATE = ['Emilia_F12_', num2str(id)];
%% NOT-HONORING EMILIA
% handles.DEPTH = -11000;
% handles.WR = 'W';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 12;

%% FROM 201 TO 300 ---> CHRISTCHURCH
% handles.GRID = 'Chch_Feb22';
% handles.GRID = 'Chch_June13';
% handles.GRID = 'Chch_Dec23';
% handles.GRID = 'Chch_Feb22';
% handles.MATE = ['Chch_Feb22', num2str(id)];
% handles.MATE = ['Chch_June13', num2str(id)];
% handles.MATE = ['Chch_Dec23', num2str(id)];
% NOT-HONORING CHRISTCHURCH
% handles.DEPTH = -5000;
% handles.WR = 'W';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 11;


%% FROM 301 TO 400 ---> WELLINGTON
% handles.GRID = 'Wellington';
% handles.MATE = ['Wellington_', num2str(id)];
% %% NOT-HONORING WELLINGTON
% handles.DEPTH = -20;
% handles.WR = 'W';
% handles.NH_MATE = 10;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 14;

%% FROM 401 TO 500 ---> MARSICA
% handles.GRID = 'Marsica';
% handles.MATE = ['Marsica_', num2str(id)];
% % NOT-HONORING MARSICA
% handles.DEPTH = -1000;
% handles.WR = 'R';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 15;

%% FROM 501 TO 600 ---> ISTANBUL
% handles.GRID = 'Istanbul_with_topo_4';
% % handles.GRID = 'Istanbul_notopo_testcase';
% handles.MATE = ['Istanbul_', num2str(id)]; 
% % NOT-HONORING ISTANBUL
% handles.DEPTH = -1300;
% handles.WR = 'R';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 16;

%% 601 ---> GRENOBLE
% handles.GRID = 'Grenoble_empty';
% handles.MATE = ['grenoble_', num2str(id)];
% % NOT-HONORING SANTIAGO
% handles.DEPTH = -3000;
% handles.WR = 'R';
% handles.NH_MATE = 0;
% handles.NH_TOL = 0;
% handles.TAG_CASE = 0;


%% 701 ---> BEIJING
% handles.GRID = 'BEIJING_SHOUNY';
% handles.MATE = ['BEIJING_SHOUNY_', num2str(id)]; 
% % NOT-HONORING BEIJING
% handles.DEPTH = -100;
% handles.WR = 'W';
% handles.NH_MATE = 1;
% handles.NH_TOL = 10;
% handles.TAG_CASE = 21;



% NOT-HONORING WITH VS30 PROFILES 0-NO, 1-YES
handles.VS30 = 1;

% SPACE DISCRETIZATION: Set 'DG ' OR 'SEM'
% handles.SPACE_DIS = 'DG ';
handles.SPACE_DIS = 'SEM';

% TIME INTEGRATION
handles.TIME_SCHEME = 'Leap-Frog';
handles.START = 0;
handles.STOP =  30;
handles.DELTAT = 0.001;
handles.TMONITOR = 10;

% OUTPUTS
handles.DISP = 1; 
handles.VEL = 0;
handles.ACC = 0;
handles.STRESS = 0;
handles.STRAIN = 0;
handles.ROTATION = 0;


% NON-LINEAR
handles.NL_MATE = '';
handles.NLE_DEPTH = '';
handles.PEAK = '';
handles.NL_DAMP = '';
handles.NL_GG0 = '';

% TABLE OF MATERIAL PROPERTIES/ BOUNDARY CONDITIONS / DG INTERFACE CONDITIONS
handles.MATE_TABLE = 'materials.ini';

% FOLDER IN WHICH THE INPUT RESULTS WILL BE STORED
handles.FOLDER = ['E00',num2str(id)];



% OPTIONS IN SPEED.input
handles.MONITORS = 'MONITORS';
handles.MPI = 'FILES_MPI';
handles.BACKUP = '';
handles.SNAP = 0;

% SLIP DISTRIBUTION
% 1 STANDARD, 2 ARCHULETA, 3 GALLOVIC
handles.SLIPDSTR = 1;

% DAMPING
% 1 FREQ PROPORTIONAL, 2 FREQ CONSTANT, 3 RAYLEIGH
handles.DAMPING = 1;



%******************************************************************************
% DEFAULT PARAMETERS : 
% DO NOT CHANGE THESE PARAMETERS FOR SELECTED SCENARIOS FROM DATABASE
%******************************************************************************
handles.MAGNITUDE = 0;
handles.SLIP = '';
handles.SCALE = 0;
handles.RANDOM = 0;
handles.RISETIME = 0;
handles.RAKE = 0;
handles.WWMAX = 0;
handles.W = 0;
handles.LLMAX = 0;
handles.L = 0;
handles.Z = 0;
handles.Y = 0;
handles.X = 0;

handles.FID = 'F00000'; 
handles.LFAULT = 0;
handles.WFAULT = 0;
handles.STRIKEFAULT = 0;
handles.DIPFAULT = 0;
handles.FOX = 0;
handles.FOY = 0;
handles.FOZ = 0;

%******************************************************************************
% DEFAULT PARAMETERS END 
% DO NOT CHANGE THESE PARAMETERS FOR SELECTED SCENARIOS FROM DATABASE
%******************************************************************************



% MAE INPUT FILES
Make_input_file_SPEED(handles)

close all;
disp('END');



