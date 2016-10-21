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


function [] = Make_input_file_SPEED(data)
%%  WRITE input files for SPEED:
%
% 1) SPEED.input
%    Header file containing the main parameters for the simulation
%
% 2) FILENAME.mate
%    File containing: Material properties
%                     Boundary conditions
%                     DG interface conditions (optional)
%                     Forcing terms (slip)
%                     Not-Honoring (optional)
%                     Non linear Dynamics (optional)
%
% 3) LS.input
%    File containing the list of monitores points
%
% 4) ALL.out/XYZ.out/VS_RS.out
%    Files containing the description of not-honoring surfaces (optional)
%

delete speed.log
warning off
addpath src
clc
delete speed.log;
diary speed.log

disp('Starting input files generation...')




%% SELECT THE OPERATING SISTEM
%UnixOrDos = 'Unix';             %ADMISSIBILE KEYWORDS: 'Unix', 'Dos '
UnixOrDos = data.OS;
if(isempty(data.OS) )
    disp('Plese Select the Operating System');
end

if UnixOrDos == 'Unix'
    dirSeparator = '/';
    disp('on Unix Filesystem')
    disp(' ')
elseif UnixOrDos == 'Dos '
    dirSeparator = '\';
    disp('on Windows Filesystem')
    disp(' ')
end


%% DEFINE SEISMIC FAULT, OUTPUT FOLDERS AND PATHNAMES

%SELECT THE SCENARIOS: Emilia, Santiago, Christchurch, Wellington
SCENARIO = data.SCENARIO_ID;
if((SCENARIO == 0) && ...
        (isempty(data.FID) ||  data.LFAULT == 0 ...
        || data.WFAULT == 0 || isempty(data.STRIKEFAULT) ...
        || isempty(data.DIPFAULT) ...
        || isempty(data.FOX) || isempty(data.FOY) ...
        || isempty(data.FOZ) ))
    disp('Plese insert a Scenario ID or fill all data for the creation of a new scenario');
    return
end

if (SCENARIO ~= 0)
    index = 0;
    while ((100*index - SCENARIO) < 0 )
        index = index + 1;
    end
    %index = index-1;
    %if (index == 2) index = 1; end   %%ONLY FOR SANTIAGO
    %if (index >= 7) index = 6; end   %%ONLY FOR ISTANBUL
    
    disp([num2str(SCENARIO), ' th scenario selected from the database'])
    disp(' ')
    
    
    fid = fopen('SCENARIOS/utmzone.dat','r');
    for iline = 1:index
        tline = fgetl(fid);
    end
    utmzone = tline(8:10);
    
    if (index/10) < 1
        loc_str = ['L00',num2str(index)];
    elseif (index/10) < 10
        loc_str = ['L0',num2str(index)];
    else
        loc_str = ['L',num2str(index)];
    end
    
    disp(['Location number: ', loc_str]);
    
    scen_par = load(['SCENARIOS/scen_par.dat']);
    row = find(scen_par(:,1) == SCENARIO);
    
    SCEN.ID = scen_par(row,2);
    if length(num2str(SCEN.ID))==1
        fault_str = ['F0000',num2str(SCEN.ID)];
    elseif length(num2str(SCEN.ID))==2
        fault_str = ['F000',num2str(SCEN.ID)];
    elseif length(num2str(SCEN.ID))==3
        fault_str = ['F00',num2str(SCEN.ID)];
    elseif length(num2str(SCEN.ID))==4
        fault_str = ['F0',num2str(SCEN.ID)];
    else
        fault_str = ['F',num2str(SCEN.ID)];
    end
    
    disp(['Fault number: ', fault_str]);
    disp(' ')
    
    FAULTNAME = 'fault';
    
    
    SCEN.SLIP = scen_par(row,3);
    SCEN.HYPO = scen_par(row,4);
    hypocenters = load(['SCENARIOS',dirSeparator,loc_str,dirSeparator,fault_str,...
        dirSeparator,'sism',dirSeparator,'hypocenter.dat']);
    
    
    SCEN.HYPO_X = hypocenters(SCEN.HYPO,1);
    SCEN.HYPO_Y = hypocenters(SCEN.HYPO,2);
    SCEN.HYPO_Z = hypocenters(SCEN.HYPO,3);
    SCEN.LOLMAX = scen_par(row,5);
    SCEN.WOWMAX = scen_par(row,6);
    SCEN.L = scen_par(row,7);
    SCEN.W = scen_par(row,8);
    SCEN.MW = scen_par(row,9);
    SCEN.VR = scen_par(row,10);
    SCEN.TAU = scen_par(row,11);
    SCEN.RAKE = scen_par(row,12);
    SCEN.ISCALE = scen_par(row,13);
    SCEN.RANDOM = scen_par(row,14);
    
    
    fault_par = ['fault.dat'];
    
    % Loading source parameters
    % source_par = input('Source parameters filename:','s');
    src_par = load(['SCENARIOS',dirSeparator,loc_str,dirSeparator,fault_str,...
        dirSeparator,'sism',dirSeparator,fault_par]);
    
    % FAULT GEOMETRY
    FAULT.ID = src_par(1);                   % fault id number
    FAULT.L = src_par(2);                    % fault length (along strike) [m]
    FAULT.W = src_par(3);                    % fault width (down dip) [m]
    FAULT.STRIKE = src_par(4);               % measured clockwise from north [deg]
    FAULT.DIP = src_par(5);                  % measured counter clockwise from horiz. strike direc.,
    % as in Aki&Richard(1980) [deg]
    FAULT.RAKE0 = src_par(6);               % rake angle (for scenario 0) [deg]

    FAULT.FO_x = src_par(7);                 % X coor of Fault Origin F0 (zero strike and zero dip) [m]
    FAULT.FO_y = src_par(8);                 % Y coor of Fault Origin F0 (zero strike and zero dip) [m]
    FAULT.FO_z = src_par(9);                 % Z coor of fault Origin F0 (zero strike and zero dip)  [m]
%     FAULT.Mw_scen0 = src_par(9);             % moment magnitude (for scenario 0)
%     FAULT.Mw_max = FAULT.Mw_scen0;
%     FAULT.slip_scen0 = src_par(10);          % slip pattern model (for scenario 0)
%     FAULT.hypo_scen0 = src_par(11);          % hypocenter (for scenario 0)
%     FAULT.velocity_scen0 = src_par(12);      % rupture velocity (for scenario 0) [m/s]
%     FAULT.tau_scen0 = src_par(13);           % rise time (for scenario 0) [s]
    
    
else
    
    SCEN.ID = 0;
    nomeslip = 'SCENARIOS/slip_models.dat';
    fids = fopen(nomeslip,'r');
    tline = fgetl(fids);
    id_slip = 0;
    found_slip = 0;
    while (found_slip == 0)
        if (length(tline) == length(data.SLIP))
            if(tline == data.SLIP)
                found_slip = 1;
            end
        end
        tline = fgetl(fids);
        id_slip = id_slip + 1;
    end
    SCEN.SLIP = id_slip;
    
    
    SCEN.HYPO = 0;
    SCEN.HYPO_X = data.X;
    SCEN.HYPO_Y = data.Y;
    SCEN.HYPO_Z = data.Z;
    SCEN.LOLMAX = data.LLMAX;
    SCEN.WOWMAX = data.WWMAX;
    SCEN.L = data.L;
    SCEN.W = data.W;
    SCEN.MW = data.MAGNITUDE;
    SCEN.VR = data.VR;
    SCEN.TAU = data.RISETIME;
    SCEN.RAKE = data.RAKE;
    SCEN.ISCALE = data.SCALE;
    SCEN.RANDOM = data.RANDOM;
    
    % FAULT GEOMETRY
    FAULT.ID = data.FID;                              % fault id number
    FAULT.L = data.LFAULT;                            % fault length (along strike) [m]
    FAULT.W = data.WFAULT;                            % fault width (down dip) [m]
    FAULT.STRIKE = data.STRIKEFAULT;                  % measured clockwise from north [deg]
    FAULT.DIP = data.DIPFAULT;                        % measured counter clockwise from horiz. strike direc.,
    % as in Aki&Richard(1980) [deg]
    
    FAULT.FO_x = data.FOX;                            % X coor of Fault Origin F0 (zero strike and zero dip) [m]
    FAULT.FO_y = data.FOX;                            % Y coor of Fault Origin F0 (zero strike and zero dip) [m]
    FAULT.FO_z = data.FOX;                            % Z coor of fault Origin F0 (zero strike and zero dip)  [m]
    FAULT.Mw_scen0 = data.MAGNITUDE;                  % moment magnitude (for scenario 0)
    FAULT.Mw_max = FAULT.Mw_scen0;
    FAULT.slip_scen0 = data.SLIP;                     % slip pattern model (for scenario 0)
    FAULT.hypo_x_scen0 = data.X;                      % hypocenter (for scenario 0)
    FAULT.hypo_y_scen0 = data.Y;
    FAULT.hypo_z_scen0 = data.Z;
    FAULT.velocity_scen0 = data.VR;                   % rupture velocity (for scenario 0) [m/s]
    FAULT.tau_scen0 = data.RISETIME;                  % rise time (for scenario 0) [s]
    FAULT.RAKE0 = data.RAKE;                          % rake angle (for scenario 0) [deg]
    
    FAULTNAME = ['fault'];
    
    fault_str = FAULT.ID;
    loc_str = 'L000';
    utmzone = '000';
end


if (isempty(data.FOLDER))
    if (SCENARIO/10 < 1)
        data.FOLDER = ['E0000',num2str(SCENARIO)];
    elseif (SCENARIO/10 < 10)
        data.FOLDER = ['E000',num2str(SCENARIO)];
    elseif (SCENARIO/10 < 100)
        data.FOLDER = ['E00',num2str(SCENARIO)];
    elseif (SCENARIO/10 < 1000)
        data.FOLDER = ['E0',num2str(SCENARIO)];
    else
        data.FOLDER = ['E',num2str(SCENARIO)];
    end
end

path = [cd,dirSeparator,'SCENARIOS',dirSeparator,loc_str,dirSeparator,fault_str,...
    dirSeparator];

mkdir([path,data.FOLDER,dirSeparator]);
path_out = [path,data.FOLDER,dirSeparator];

mkdir([path_out,'input4speed'])
path_out_speed = [path_out,'input4speed',dirSeparator];

mkdir([path,data.FOLDER,dirSeparator,'figures',dirSeparator]);
path_f = [path,data.FOLDER,dirSeparator,'figures',dirSeparator];

path_in = [path,'cub',dirSeparator];

disp(' ')
disp('---------------------------------------------------------------')
disp('SPEED INPUT FILES WILL BE STORED IN FOLDER:')
disp(' ')
disp(path_out)
disp(' ')
disp('---------------------------------------------------------------')
disp('Scenario properties:')
disp(['SLIP ID: ', num2str(SCEN.SLIP)])
disp(['Mw: ', num2str(SCEN.MW)])
disp(['HYPO ID: ', num2str(SCEN.HYPO)])
disp(['HYPO LOC: ' num2str(SCEN.HYPO_X), '  ',num2str(SCEN.HYPO_Y), '  ',num2str(SCEN.HYPO_Z)])
disp(['L [m]: ', num2str(SCEN.L), ' L0/LMAX: ', num2str(SCEN.LOLMAX)])
disp(['W [m]: ', num2str(SCEN.W), ' W0/WMAX: ', num2str(SCEN.WOWMAX)])
disp(['RUPTURE VEL [m/s]: ', num2str(SCEN.VR)])
disp(['RISE TIME [s]: ', num2str(SCEN.TAU)])
disp(['RAKE [deg]: ', num2str(SCEN.RAKE)])
disp(['SCALING: ', num2str(SCEN.ISCALE)])
disp('---------------------------------------------------------------')

disp('Fault properties:' )
disp(['FAULT ID: ', num2str(FAULT.ID)])
disp(['L MAX [m]: ', num2str(FAULT.L)])
disp(['W MAX [m]: ', num2str(FAULT.W)])
disp(['STRIKE [deg]: ', num2str(FAULT.STRIKE)])
disp(['DIP [deg]: ', num2str(FAULT.DIP)])
disp(['F0X: ',num2str(FAULT.FO_x),'   F0Y: ',num2str(FAULT.FO_y),'   F0Z: ',num2str(FAULT.FO_z)])
disp('---------------------------------------------------------------')
disp(' ')



%% DEFINITION OF INPUT PARAMETERS TO GENERATE MATFILE.MATE FILES

%DEFINE THE NAME OF THE OUTPUT FILE CONTAINING THE MATERIAL DESCRIPTION
MATFILE = data.MATE;

%DEFINE THE NAME OF THE OUTPUT FILE CONTAINING THE GRID
GRIDFILE = data.GRID;
disp(' ')
disp(['Reading material file: ',path,'mate', dirSeparator,MATFILE,'.mate'])


%FILE CONTRINING THE DESCRIPTION OF THE MATERIAL PROPERTIES TO BE READ FOR THE SELECTED MODEL
MATERIALS = data.MATE_TABLE;        %The file has to be store in the folder 'mate'
if (isempty(data.MATE_TABLE))
    disp('Material file not found!!!')
    disp('Insert material filename before running')
end


%SELECT THE NOT-HONORING CASE DESCRIBING THE SCENARIO
%TAGCASE: 11-Christchurch, 12-Emila Romagna, ??-Santiago

TAGCASE = data.TAG_CASE;
NH_TOL =data.NH_TOL;                %Tolerance value for finding nodes inside the alluvial basin
%according to the not-honoring technique

if(TAGCASE~=0)
    disp(' ')
    disp(['Not honoring case: ', num2str(TAGCASE)])
    disp(['Tolerance: ', num2str(NH_TOL)])
end

%SELECT THE NON LINEAR BEHAVIOR
%FILES CONTAINING THE G/G0 AND DAMPING CURVES FOR NONLINEAR MATERIALS (if present)

CURVE_NLE1 = data.NL_GG0;       %The file has to be store in the folder 'mate'
TYPE_FUNC_NLE1 = 60;            %Time function type for the curve G/G0

CURVE_NLE2 = data.NL_DAMP;      %The file has to be store in the folder 'mate'
TYPE_FUNC_NLE2 = 61;            %Time function type for the damping curve

DEPTH_NLE = data.NLE_DEPTH;     %Depth (in meters) from the surface in which a nonlinear dynamics is applied
TAG_FUNC_NLE = 11;              %Tag for the functions described in CURVE_NLE1 and CURVE_NLE2


if(DEPTH_NLE ~=0)
    disp(' ')
    disp('Non linear elasticity: ')
    disp(['Reading '])
    disp([path,'mate', dirSeparator, data.NL_GG0])
    disp([path,'mate', dirSeparator,data.NL_DAMP])
    
end





%% DEFINITION OF INPUT PARAMETERS TO GENERATE SPEED.input FILE
MPIFILE = data.MPI; %'FILES_MPI';          %Folder in which files *.mpi are stored
MONFILE = data.MONITORS;  %'MONITORS';           %Folder in which MONITOR* files are stored
BKPFILE = data.BACKUP;

DGMODEL = 0;                        %0-SEM discretization, 1-DGSEM discretization
if (data.SPACE_DIS == 'DG ')
    DGMODEL = 1;
    DGMETHOD = -1;                  %Select the DG method (-1 SIPG, 0 IIPG, 1 NIPG)
    PENALIZC = 250;                 %Constant value for the Penalty term in the DG method
end

OPTIOUT = [data.DISP data.VEL data.ACC data.STRESS data.STRAIN data.ROTATION 1 6];    %Option for the output (See Manual)

TIMESTEP = data.DELTAT;                     %Time step used for the time scheme
TMONITOR = data.TMONITOR;                   %The solution is saved every TMONITOR step
                                            %(or every discrete times = TIMESTEP*TMONITOR)
STOPTIME = data.STOP; %50;                  %Final time
STARTIME = data.START;

MLST_depth = data.DEPTH; %-10;               %Depth (in meters) from which start finding monitored points
MLST_tag = 0;                                %Tag for MLST.input file. 0-Write, 1-Read MLST.input file
if (data.WR == 'R') MLST_tag = 1; end

NB_SNAP = data.SNAP;
if (data.SNAP ~= 0) SNAPSHOTS = linspace(data.START,data.STOP,NB_SNAP+1); end



%% FIRST CONVERT EXODUS FILE IN cub FOLDER TO ASCII FILES
% START
% Reading file from cubit
% File .mesh, ALL.out, XYZ.out generation

% SELECT INPUT FILES:
% List of files exported from CUBIT (.e). All these files have to be
% stored in the folder 'cub'

ALLUVIAL = 'ALL';               %File containing the mesh of the alluvial basin
TOPO = 'XYZ';                   %File containing the mesh of the topography
MONITOR = 'monitors';           %File containing the list of monitored point, here definied as a mesh

disp('---------------------------------------------------------------')

disp('Start reading files exported from cubit...')
disp(['Opening folder :' path,'cub'])
disp(' ')
dummy = 0;
if (UnixOrDos == 'Dos ')
    disp('Do you want to convert EXODUS files in ASCII files? (1=YES)');
    dummy = 0;
    dummy = input('[default = NO]: ');
    
    if (dummy == 1)
        file1 = [path,'cub',dirSeparator,ALLUVIAL,'.e'];
        file2 = [path,'cub',dirSeparator,ALLUVIAL,'.txt'];
        system(['ncdumps.exe ',file1,' > ',file2]);
        
        file1 = [path,'cub',dirSeparator,TOPO,'.e'];
        file2 = [path,'cub',dirSeparator,TOPO,'.txt'];
        system(['ncdumps.exe ',file1,' > ',file2]);
        
        file1 = [path,'cub',dirSeparator,MONITOR,'.e'];
        file2 = [path,'cub',dirSeparator,MONITOR,'.txt'];
        system(['ncdumps.exe ',file1,' > ',file2]);
        
        file1 = [path,'cub',dirSeparator,GRIDFILE,'.e'];
        file2 = [path,'cub',dirSeparator,GRIDFILE,'.txt'];
        system(['ncdumps.exe ',file1,' > ',file2]);
    end
else
    dummy = 0;
    disp('Do you want to convert EXODUS files in ASCII files? (1=YES)');
    dummy = 0;
    %dummy = input('[default = NO]: ');
    if (dummy == 1)
        disp('Please convert EXODUS files before running this program!')
        return
    end
end



%CONVERT GRID FROM .txt TO .mesh
dummy = 0;
disp('Do you want to convert ASCII files to .mesh files ? (1=YES)');
dummy = 0;
%dummy = input('[default = NO]: ');
if dummy == 1 %(size(dummy) ~= 0)
    exo2gid_3D_tria(path_in,path_in,ALLUVIAL);
    file1 = [path,'cub',dirSeparator,ALLUVIAL,'.out'];
    copyfile(file1,path_out_speed);
    
    exo2gid_3D_tria(path_in,path_in,TOPO);
    file1 = [path,'cub',dirSeparator,TOPO,'.out'];
    copyfile(file1,path_out_speed);
    
    exo2ucdx_3D_nconf(path_in,path_in,GRIDFILE);
    file1 = [path,'cub',dirSeparator,GRIDFILE,'.mesh'];
    copyfile(file1,path_out_speed);
    
else
    file1 = [path,'cub',dirSeparator,ALLUVIAL,'.out'];
    %copyfile(file1,path_out_speed);
    file1 = [path,'cub',dirSeparator,TOPO,'.out'];
    %copyfile(file1,path_out_speed);
    file1 = [path,'cub',dirSeparator,GRIDFILE,'.mesh'];
    %copyfile(file1,path_out_speed);
    
end

if (data.VS30 == 1)
    file1 = [path,'cub',dirSeparator,'VS_RS.out'];
    %copyfile(file1,path_out_speed);
end


% END

%% FILE FILENAME.mate GENERATION

% START
nomefiler = [path,'mate',dirSeparator,MATERIALS];
file_mate = [path_out_speed,MATFILE,'.mate'];

fid = fopen(nomefiler,'r');

% READ NUMBER OF MATERIALS
nmate = 0; boundary_cond = 0;
tline = fgetl(fid);
tline = fgetl(fid);
while  strcmp(tline(1:4),'STOP') ~= 1
    tline = fgetl(fid);
    length(tline)
    if(length(tline) >= 50)
        nmate = nmate+1;
    else
        boundary_cond = boundary_cond +1;
    end
    
end
nmate = nmate-1;
nmate_all = nmate + boundary_cond;
fclose(fid);

fid = fopen(nomefiler,'r');
tline = fgetl(fid);
tline = fgetl(fid);
mate = zeros(nmate-1,9);
fidw = fopen(file_mate,'w');

for k=1:nmate
    tline = fgetl(fid);
    scratch = str2num(tline);
    mate_no = scratch(1);
    SD  = scratch(2);
    rho = scratch(3);
    Vp = scratch(4);
    Vs = scratch(5);
    Qp = scratch(6);
    Qs = scratch(7);
    fmax = scratch(8);  
       
    
    fprintf(fidw,['MATE    %d    %d   %6.2f  %6.2f    %6.2f    %d    %d \n'],mate_no,...
        SD,  rho, Vs, Vp, Qs, Qp );
    
end
fprintf(fidw,'\n');
fprintf(fidw,['ABSO    %d'],[mate_no+1]);
fprintf(fidw,'\n');

if(DGMODEL == 1)
    for i = mate_no+2 : nmate_all
        fprintf(fidw,['DGIC    %d %d'],[i mod(i - mate_no,2)]);
        fprintf(fidw,'\n');
    end
end

fprintf(fidw,'\n');
fprintf(fidw,['FMAX    %d'], fmax);
fprintf(fidw,'\n');


fclose(fid);

MATE_NLE = data.NL_MATE;
MATE_NH = data.NH_MATE;

if (MATE_NH ~=0 )
    fprintf(fidw,['CASE %d  %d %3.2f \n'],TAGCASE,MATE_NH,NH_TOL);
    fprintf(fidw,'\n');
end

if (MATE_NLE ~= 0)
    
    nomefile_1 = [path,'mate',dirSeparator,CURVE_NLE1];
    nomefile_2 = [path,'mate',dirSeparator,CURVE_NLE2];
    
    lett1 = load(nomefile_1);
    strain1 = lett1(:,1);
    G = lett1(:,2);
    npun1 = length(strain1);
    
    lett2 = load(nomefile_2);
    strain2 = lett2(:,1);
    damp = lett2(:,2);
    npun2 = length(strain2);
    
    figure(51)
    set(gca,'fontsize',18);
    semilogx(strain1*100,G,'r','linewidth',3);
    xlabel('\gamma (%)');
    ylabel('G/G0 (-)');
    %grid on
    hold on
    set(gca,'fontsize',18);
    semilogx(strain2*100,damp,':b','linewidth',3);
    xlabel('\gamma (%)');
    ylabel('\zeta (-)');
    grid on
    
    fprintf(fidw,['MATN %d %d %d %3.2f\n'],MATE_NLE,mate(MATE_NLE,2),TAG_FUNC_NLE,DEPTH_NLE);
    fprintf(fidw,'\n');
    fprintf(fidw,['FPEK %6.4f\n'],data.PEAK);
    fprintf(fidw,'\n');
    
    fprintf(fidw,'FUNC %d %d %d ',[TAG_FUNC_NLE,TYPE_FUNC_NLE1,npun1]);
    for k = 1:npun1
        fprintf(fidw,'%13.9f %13.9f',[strain1(k) G(k)]);
    end
    fprintf(fidw,'\n');
    
    
    fprintf(fidw,'FUNC %d %d %d ',[TAG_FUNC_NLE,TYPE_FUNC_NLE2,npun2]);
    for k = 1:npun1
        fprintf(fidw,'%13.9f %13.9f',[strain2(k) damp(k)]);
    end
    fprintf(fidw,'\n');
    
end

% % END


%% GENERATION OF SEISMIC FAULT

if (UnixOrDos == 'Dos ')
    file1 = [path,'cub',dirSeparator,FAULTNAME,'.e'];
    file2 = [path,'cub',dirSeparator,FAULTNAME,'.txt'];
    system(['ncdumps.exe ',file1,' > ',file2]);
end

%file1 = [path,'cub',dirSeparator,FAULTNAME,'.e'];
file2 = [path,'cub',dirSeparator,FAULTNAME,'.txt'];
copyfile(file2,[path,'sism',dirSeparator])

if(data.SLIPDSTR==1)
    fprintf(fidw,'SLIP %s \n', 'STD');
elseif(data.SLIPDSTR==2)
    fprintf(fidw,'SLIP %s \n', 'ARC');
elseif(data.SLIPDSTR==3)
    fprintf(fidw,'SLIP %s \n', 'GAL');
end
fprintf(fidw,'\n');
    

if (data.SLIPDSTR == 1)
    
    %SOURCE TIME FUNCTION FOR THE EXTENDED FAULT
    TAG_FUNC_STF = 15;              %Tag for the selected time function (See Manual for the Keyword FUNC)
    TYPE_FUNC_STF = 50;             %Type of the selected time function (See Manual for the Keyword FUNC)
    
    %DEFINE PARAMETES FOR SELECTED TIME FUNCTOIN FOR THE EXTENDED FAULT (See Manual)
    switch TYPE_FUNC_STF
        case 50
            PAR_1 = 0.7;
            PAR_2 = 2;
        otherwise
            PAR_1 = 0;
            PAR_2 = 0;
    end
    
    
    % .SISM for fault model
    fprintf(fidw,'FUNC %d %d %4.2f %4.2f\n',[TAG_FUNC_STF,TYPE_FUNC_STF,PAR_1,PAR_2]);
    fprintf(fidw,'\n');
    
    scen_generator(path,path_f,path_out,dirSeparator,utmzone,loc_str,...
        MATFILE,fault_str,SCENARIO,...
        TAG_FUNC_STF,SCEN,FAULT,data.FOLDER,TAGCASE)
    
else
    
    %SOURCE TIME FUNCTION FOR THE EXTENDED FAULT
    TAG_FUNC_STF = 21;              %Tag for the selected time function (See Manual for the Keyword FUNC)
    TYPE_FUNC_STF = 55;             %Type of the selected time function (See Manual for the Keyword FUNC)
    
    %DEFINE PARAMETES FOR SELECTED TIME FUNCTOIN FOR THE EXTENDED FAULT (See Manual)
    switch TYPE_FUNC_STF
        case 55
            PAR_1 = 0;
            PAR_2 = 0;
        otherwise
            PAR_1 = 0;
            PAR_2 = 0;
    end
    
    
    % .SISM for fault model
    fprintf(fidw,'FUNC %d %d %4.2f %4.2f\n',[TAG_FUNC_STF,TYPE_FUNC_STF,PAR_1,PAR_2]);
    fprintf(fidw,'\n');
    
    scen_generator_archuleta(path,path_f,path_out,dirSeparator,utmzone,loc_str,...
        MATFILE,fault_str,SCENARIO,...
        TAG_FUNC_STF,SCEN,FAULT,data.FOLDER,TAGCASE)
    
    
    
end

%
% fclose(fidw);
%


%% FILE LS.input GENERATION

exo2gid_3D_tria([path_in],path_out_speed,MONITOR); % cambia qui non in output
file_monitor = [path_out_speed,'LS.input'];


disp(' ' )
disp(['Write LS.input from ', MONITOR, '.txt'])

% START  %%%%%%%%%%%%%%
fid = fopen([path_in,MONITOR,'.monitor'],'r');
tline = fgetl(fid);
nmon = 0;
while  strcmp(tline(1:3),'MON') == 1
    tline = fgetl(fid);
    nmon = nmon+1;
end
nmon = nmon-1;
fclose(fid);

fid = fopen([path_in,MONITOR,'.monitor'],'r');
tline = fgetl(fid);
Coord_utm = zeros(nmon,3);
fidw2 = fopen(file_monitor,'w');
fprintf(fidw2,'%d\n',[nmon]);
for k=1:nmon
    tline = fgetl(fid);
    trash = tline(1:9);
    %DOS X = str2num(tline(10:26));
    %DOS Y = str2num(tline(27:43));
    %DOS Z = str2num(tline(44:58));
    X = str2num(tline(10:25));
    Y = str2num(tline(26:41));
    Z = str2num(tline(42:55));
    
    Coord_utm(k,1) = X;
    Coord_utm(k,2) = Y;
    Coord_utm(k,3) = Z;
    fprintf(fidw2,'%d   %e   %e   %e \n',[k,Coord_utm(k,1),Coord_utm(k,2),...
        Coord_utm(k,3)]);
    
end
fclose(fidw2);
fclose(fid);

% system(['del ',path_out_speed,dirSeparator,'monitors.out']);
delete([path_out_speed,dirSeparator,'monitors.out']);

% END

%% FILE SPEED.input GENERATION
file_speed = [path_out_speed,'SPEED.input'];
disp('  ')
disp('---------------------------------------------------------------')
disp('Start writing SPEED.input file')


% START
fidw3 = fopen(file_speed,'w');
fprintf(fidw3,'GRIDFILE  %s\n',GRIDFILE);
fprintf(fidw3,'MATFILE   %s\n',MATFILE);
if (MPIFILE ~= ' ')
    fprintf(fidw3,'MPIFILE   %s\n',MPIFILE);
    command = ['mkdir ',path_out_speed,MPIFILE];
    system(command);
end

if (MONFILE ~= ' ')
    fprintf(fidw3,'MONFILE   %s\n',MONFILE);
    command = ['mkdir ',path_out_speed,MONFILE];
    system(command);
end
if (BKPFILE ~= ' ')
    fprintf(fidw3,'BKPFILE   %s\n',BKPFILE);
    command = ['mkdir ',path_out_speed,BKPFILE];
    system(command);
end


if(DGMODEL == 1) fprintf(fidw3,'DGMETHOD   %d\n',DGMETHOD); end;
if(DGMODEL == 1) fprintf(fidw3,'PENALIZC   %d\n',PENALIZC); end;

fprintf(fidw3,'\n');
fprintf(fidw3,'OPTIOUT  %d %d %d  %d %d %d    %d %d\n',OPTIOUT(1),...
    OPTIOUT(2),OPTIOUT(3),OPTIOUT(4),OPTIOUT(5),OPTIOUT(6),...
    OPTIOUT(7),OPTIOUT(8));
fprintf(fidw3,'\n');

fprintf(fidw3,'DAMPING  %d \n', data.DAMPING);

if(data.TIME_SCHEME(1:5) == 'Runge')  fprintf(fidw3,'TIMESCHM RUNGEKUTTA 4 4 \n'); end


fprintf(fidw3,'\n');

if (STARTIME ~= 0) fprintf(fidw3,'STARTIME   %6.2f\n',STARTIME); end
fprintf(fidw3,'TIMESTEP   %6.5f\n',TIMESTEP);
fprintf(fidw3,'TMONITOR   %d\n',TMONITOR);
fprintf(fidw3,'STOPTIME   %6.2f\n',STOPTIME);
fprintf(fidw3,'\n');

fprintf(fidw3,'MLST   %8.2f  %d\n',MLST_depth,MLST_tag);

if (data.SNAP ~= 0)
    for i = 2 : NB_SNAP+1
        fprintf(fidw3,'SNAPSHOT   %6.2f\n',SNAPSHOTS(i));
    end
end
fclose(fidw3);

% END
% pause
disp('END')
diary 'speed.log'

