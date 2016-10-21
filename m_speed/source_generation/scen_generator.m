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


function scen_generator(path,path_f,path_out,dirSeparator,utmzone,loc_str,...
    file_mate,fault_str,SCENARIO,...
    nfunc,SCEN,FAULT,dirname,tagcase)
%
% Generation of the slip, rake angle, rupture time and rise time
% distributions on the fault plane according to your data

disp('---------------------------------------------------------------')
disp(['Generating Extended Fault...'])

file_name = ['fault'];
prop_par  = ['mec_prop.dat'];


% ======================================================================= %
% ========================= INPUT PARAMETERS ============================ %
% ======================================================================= %
% Loading source parameters

% FAULT GEOMETRY
L = FAULT.L;
W = FAULT.W;
STRIKE = FAULT.STRIKE;
DIP = FAULT.DIP;
FO_x = FAULT.FO_x;
FO_y = FAULT.FO_y;
FO_z = FAULT.FO_z;
% Mw_scen0 = FAULT.Mw_scen0;
% Mw_max = Mw_scen0;

% Fault Origin lat-lon
[Lat_FO,Lon_FO] = utm2wgs(FO_x,FO_y,utmzone);

% SEISMIC MOMENT TENSOR
% SCENARIO

if SCENARIO ~= 0
    
    scen_num =  SCENARIO;
    %fault_id =  SCEN.ID;
    slip_scen = SCEN.SLIP;
    Hypocenter = SCEN.HYPO;
    L0_scen_norm = SCEN.LOLMAX;
    W0_scen_norm = SCEN.WOWMAX;
    L_scen = SCEN.L;
    W_scen = SCEN.W;
    Mw = SCEN.MW;
    velocity = SCEN.VR;
    tau = SCEN.TAU;
    RAKE = SCEN.RAKE;
    iscale = SCEN.ISCALE;
    icrsp = SCEN.RANDOM;
else
    scen_num = 0;
    %     slip_scen = 0;
    Mw = SCEN.MW;
    velocity = SCEN.VR;
    tau = SCEN.TAU;
    RAKE = SCEN.RAKE;
    slip_scen = SCEN.SLIP;
    Hypocenter = SCEN.HYPO;
    iscale = SCEN.ISCALE;
    icrsp = SCEN.RANDOM;
end

Magnitude0 = 10^( 3/2*(Mw+6) );  % M0 in Nm
%----------------------------------------------------------------------
% Mechanical Properties (hyp.: we assume flat layered stratigraphy in
%                        where the fault is located):
%----------------------------------------------------------------------
% list_prop = input('Mechanical Properties filename:','s');
prop = load([char(path),'sism',dirSeparator,prop_par]);
dim = size(prop);
nlay = dim(1);

disp('Mechanical properties:')
disp(' ')
disp('Number of layers:');
disp(num2str(nlay));

Vp = prop(:,1);
Vs = prop(:,2);
rho = prop(:,3);

for i=1:nlay
    disp(['Layer #',num2str(i)]);
    disp(['Vp = ',num2str(Vp(i)),' m/s']);
    disp(['Vs = ',num2str(Vs(i)),' m/s']);
    disp(['rho = ',num2str(rho(i)),' kg/m3']);
end

% HYPOCENTER LOCATION
XC =  SCEN.HYPO_X;
YC =  SCEN.HYPO_Y;
ZC =  SCEN.HYPO_Z;

% hypocenter lat-lon
[Lat_HY,Lon_HY] = utm2wgs(XC,YC,utmzone);

% Slip vector (uses slip_vec.m)
slip_mu = slip_vec(STRIKE,RAKE,DIP);
slip_x_mu = slip_mu(1);
slip_y_mu = slip_mu(2);
slip_z_mu = slip_mu(3);

% Fault normal vector (uses norm_vec.m)
norm_mu = norm_vec(STRIKE,RAKE,DIP);
norm_x_mu  = norm_mu(1);
norm_y_mu  = norm_mu(2);
norm_z_mu  = norm_mu(3);

disp(' ')
disp('Loading slip distribution')
% SLIP DISTRIBUTION
nomeslip = [cd,dirSeparator,'SCENARIOS',dirSeparator,'slip_models.dat'];
fids = fopen(nomeslip,'r');

for is = 1:slip_scen
    tline = fgetl(fids);
end
nomer2 = char(tline);

disp(['Slip pattern model: ',nomer2]);
disp(' ')
fclose(fids);

nomefiler2 = [cd,dirSeparator,'SCENARIOS', dirSeparator,'sism.srcmod',dirSeparator,nomer2,'.srcmod'];


% % SLIP DISTRIBUTION (Orientation)
% disp('                                         ');
% disp('                                         ');
% disp('---------------------------------------------');
% disp('Orientation of the Slip Pattern...');
%
% default = 1;
% disp('distribution (ORIG=1, LR=2, UD=3, UD(LR)=4) ');
% dummy = input('[default = 1]: ');
% if size(dummy) == 0
%     type_orient = default;
% else
%     type_orient = dummy;
% end
%
% disp('---------------------------------------------');
% disp('                                         ');
%
% if type_orient==1
%     disp('    FROM        ->          TO');
%     disp('+---------+             +---------+');
%     disp('| 1  |  2 |             | 1  |  2 |');
%     disp('+---------+    ORIG     +---------+');
%     disp('| 4  |  3 |             | 4  |  3 |');
%     disp('+---------+             +---------+');
%
% elseif type_orient==2
%     disp('    FROM        ->          TO');
%     disp('+---------+             +---------+');
%     disp('| 1  |  2 |             | 2  |  1 |');
%     disp('+---------+     LR      +---------+');
%     disp('| 4  |  3 |             | 3  |  4 |');
%     disp('+---------+             +---------+');
%
% elseif type_orient==3
%
%     disp('    FROM        ->          TO');
%     disp('+---------+             +---------+');
%     disp('| 1  |  2 |             | 4  |  3 |');
%     disp('+---------+     UD      +---------+');
%     disp('| 4  |  3 |             | 1  |  2 |');
%     disp('+---------+             +---------+');
%
% elseif type_orient==4
%     disp('    FROM        ->          TO');
%     disp('+---------+             +---------+');
%     disp('| 1  |  2 |             | 3  |  4 |');
%     disp('+---------+   UD(LR))   +---------+');
%     disp('| 4  |  3 |             | 2  |  1 |');
%     disp('+---------+             +---------+');
% end
%
% disp('                                         ');
% disp('---------------------------------------------');
% disp('---------------------------------------------');

% CORRELATED RANDOM SOURCE PARAMETERS
if icrsp == 1
    nomecrsp = [cd,dirSeparator,'SCENARIOS',dirSeparator,'crsp.dat'];
    crsp_param = load(nomecrsp);
    
    % rake
    tag_crsp_rake = crsp_param(1);
    cov_rake = crsp_param(2);
    type_psd_rake = crsp_param(3);
    eta_rake = crsp_param(4);
    mean_rake = RAKE;
    std_rake = cov_rake.*mean_rake;
    % rise time
    tag_crsp_rt = crsp_param(5);
    cov_rt = crsp_param(6);
    type_psd_rt = crsp_param(7);
    eta_rt = crsp_param(8);
    mean_rt = tau;
    std_rt = cov_rt.*mean_rt;
    
    %rupture velocity
    tag_crsp_vr = crsp_param(9);
    cov_vr = crsp_param(10);
    type_psd_vr = crsp_param(11);
    eta_vr = crsp_param(12);
    vel_threshold = crsp_param(13);
    mean_vr = velocity;
    std_vr = cov_vr.*mean_vr;
    
    ax = crsp_param(14);
    ay = crsp_param(15);
    % ax default = 10^(-2.5+1/2*Mw)*1000; % Mai & Beroza (2002,JGR)
    % ay default = 10^(-1.5+1/3*Mw)*1000; % Mai & Beroza (2002,JGR)
else
    tag_crsp_rake = 0;
    tag_crsp_rt = 0;
    tag_crsp_vr = 0;
end

% $$$$$$$$$$$$$$$ %

tag_crsp_vr = 0;




if length(num2str(scen_num))==1
    eqk_str = ['E0000',num2str(scen_num)];
elseif length(num2str(scen_num))==2
    eqk_str = ['E000',num2str(scen_num)];
elseif length(num2str(scen_num))==3
    eqk_str = ['E00',num2str(scen_num)];
elseif length(num2str(scen_num))==4
    eqk_str = ['E0',num2str(scen_num)];
else
    eqk_str = ['E',num2str(scen_num)];
end

if length(num2str(slip_scen))==1
    slip_str = ['S000',num2str(slip_scen)];
elseif length(num2str(slip_scen))==2
    slip_str = ['S00',num2str(slip_scen)];
elseif length(num2str(slip_scen))==3
    slip_str = ['S0',num2str(slip_scen)];
elseif length(num2str(slip_scen))==4
    slip_str = ['S',num2str(slip_scen)];
end



nomefiler=[char(path),'sism',dirSeparator,char(file_name),'.txt'];
mkdir([char(path_out),dirSeparator,'sism']);
nomefilew4=[char(path_out),dirSeparator,'sism',dirSeparator,char(eqk_str),'.SISM'];
nomefilew5=[char(path_out),dirSeparator,'sism',dirSeparator,char(eqk_str),'.param'];
nomefilew6=[char(path),'scenario_',fault_str,'.param'];
nomefile_mate=[char(path),dirname,dirSeparator,'input4speed',dirSeparator,char(file_mate),'.mate'];

fid = fopen(nomefiler,'r');
tini=cputime;

disp(['Reading and Converting ', nomefiler])
disp(' ' )
tab=7;
ctria=0;
grid=0;
num_ctria=0;

for i= 1:8
    tline = fgetl(fid);
end

num_nodes = str2num(tline(20-tab:length(tline)-1));
tline = fgetl(fid);
num_elem = str2num(tline(19-tab:length(tline)-1));
tline = fgetl(fid);
num_el_blk=str2num(tline(21-tab:length(tline)-1));
tline = fgetl(fid);
tline = fgetl(fid);

for i=1:num_el_blk
    num_el_in_blk(i)=str2num(tline(25-tab+fix(i/10):length(tline)-1));
    tline = fgetl(fid);
    num_nod_per_el(i)=str2num(tline(26-tab+fix(i/10):length(tline)-1));
    tline = fgetl(fid);
    if num_nod_per_el(i)==3
        num_ctria = num_ctria +  num_el_in_blk(i);
    end
    tline = fgetl(fid);
end

num_ctria;

l_tline=8;
tline = fgetl(fid);
i=0;

con_ctria=zeros(num_ctria,3);

grid_cord=zeros(num_nodes*3,1);

ctria_tag=zeros(num_ctria,1);
ctria_id=zeros(num_ctria,1);

grid_cord=zeros(num_nodes,3);
grid_id=zeros(num_nodes,1);
grid_x=zeros(num_nodes,1);
grid_y=zeros(num_nodes,1);
grid_z=zeros(num_nodes,1);

ctria_1=zeros(num_ctria,1);
ctria_2=zeros(num_ctria,1);
ctria_3=zeros(num_ctria,1);


for i=1:num_el_blk
    
    if i==1
        tline = fgetl(fid);
        %         disp([' '])
        disp(['Reading ',char(num2str(i)),'st block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
    elseif i==2
        disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
        %         disp([' '])
        disp(['Reading ',char(num2str(i)),'nd block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
    elseif i==3
        disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
        %         disp([' '])
        disp(['Reading ',char(num2str(i)),'rd block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
    else
        disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
        %         disp([' '])
        disp(['Reading ',char(num2str(i)),'th block  (',char(num2str(num_el_in_blk(i))),' elem.)'])
    end
    
    tini=cputime;
    
    while strcmp(tline(1:l_tline),' connect')~=1
        tline = fgetl(fid);
        if length(tline)<8
            l_tline=length(tline);
        elseif length(tline)==0
            l_tline=1;
        else
            l_tline=8;
        end
    end
    
    if num_nod_per_el(i)==3
        for j=1:num_el_in_blk(i)
            tline = fgetl(fid);
            %pos=findstr(tline,',');
            ctria=ctria+1;
            ctria_id(ctria)=ctria;
            ctria_tag(ctria)=i;
            con_ctria(ctria,:)=str2num(tline);
        end
    end
end

% cont tag materials
num_tag_ctria = 1;
for i = 1:ctria-1
    if ctria_tag(i) ~= ctria_tag(i+1)
        num_tag_ctria = num_tag_ctria +1;
    end
end


disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
% disp([' '])

tini=cputime;
% disp([' '])
% disp([' '])
disp(['*** Storing data informations ***'])

if ctria>0
    ctria_1(1:ctria)=con_ctria(1:ctria,1);
    ctria_2(1:ctria)=con_ctria(1:ctria,2);
    ctria_3(1:ctria)=con_ctria(1:ctria,3);
end

disp(['completed in ',char(num2str(cputime-tini)),' sec.'])
% disp([' '])

tini=cputime;
% disp([' '])
disp(['BEGIN - Reading nodes coordinates'])


tline = fgetl(fid);
tline = fgetl(fid);
vero=[];
grid=1;


while strcmp(tline(1:l_tline),' coord =')~=1
    tline = fgetl(fid);
    if length(tline)<8
        l_tline=length(tline);
    elseif length(tline)==0
        l_tline=1;
    else
        l_tline=8;
    end
end

while length(vero)==0
    tline = fgetl(fid);
    vero=findstr(tline,';');
    
    if length(vero)==0
        howmanynumb=length(str2num(tline));
        grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
        grid=grid+howmanynumb;
        
    else
        howmanynumb=length(str2num(tline));
        grid_cord(grid:grid+howmanynumb-1)=str2num(tline);
        grid=grid+howmanynumb;
    end
end


grid=num_nodes;
for i=1:num_nodes
    grid_id(i)=i;
end
grid_x(1:num_nodes)=grid_cord(1:num_nodes);
grid_y(1:num_nodes)=grid_cord(num_nodes+1:2*num_nodes);
grid_z(1:num_nodes)=grid_cord(2*num_nodes+1:3*num_nodes);

fclose(fid);


disp(['END - Reading nodes coordinates in ',char(num2str(cputime-tini)),' sec.'])
% disp([' '])


%---------------------------------------------------------------------
% FAULT
%----------------------------------------------------------------------
%     fid = fopen(nomefilew,'w');
%
%     blank_space=length(num2str(grid_id(grid)));
%     fprintf(fid,'#mesh of a table\n');
%     fprintf(fid,'MESH "model" Dimension 3 ElemType Triangle Nnode 3\n');
%     fprintf(fid,'Coordinates\n');
%     fprintf(fid,'# node_number coordinate_x coordinate_y coordinate_z\n');
%

tini=cputime;
% disp([' '])
disp(' ')
disp(['BEGIN - Writing input format'])

figure(1);
set(gca,'fontsize',16);
plot3(grid_x,grid_y,grid_z,'b.','markersize',8);hold on;
plot3(XC,YC,ZC,'pr','linewidth',2.0,'markersize',20);hold on;
xlabel('X - East (m)','fontsize',16);
ylabel('Y - North (m)','fontsize',16);
zlabel('Z - Elevation (m)','fontsize',16);
hold on
plot3(grid_x(1),grid_y(1),grid_z(1),'^r');
text(grid_x(1),grid_y(1),grid_z(1),'1','color','r');
plot3(grid_x(2),grid_y(2),grid_z(2),'^r');
text(grid_x(2),grid_y(2),grid_z(2),'2','color','r');
plot3(grid_x(3),grid_y(3),grid_z(3),'^r');
text(grid_x(3),grid_y(3),grid_z(3),'3','color','r');

dist12 = sqrt((grid_x(1)-grid_x(2))^2+(grid_y(1)-grid_y(2))^2 + ...
    (grid_z(1)-grid_z(2))^2);
dist13 = sqrt((grid_x(1)-grid_x(3))^2+(grid_y(1)-grid_y(3))^2 + ...
    (grid_z(1)-grid_z(3))^2);
dist23 = sqrt((grid_x(2)-grid_x(3))^2+(grid_y(2)-grid_y(3))^2 + ...
    (grid_z(2)-grid_z(3))^2);

dist = min([dist12 dist13 dist23]);
dist = round(dist./2);

x_vec1 = grid_x(2)-grid_x(1);
y_vec1 = grid_y(2)-grid_y(1);
z_vec1 = grid_z(2)-grid_z(1);
x_vec2 = grid_x(2)-grid_x(3);
y_vec2 = grid_y(2)-grid_y(3);
z_vec2 = grid_z(2)-grid_z(3);
vec1 = [x_vec1,y_vec1,z_vec1];
vec2 = [x_vec2,y_vec2,z_vec2];
normale1 = cross(vec2,vec1);
ass  = norm(normale1);
normale = normale1/ass;
check = normale' - norm_mu;
checki = -normale' - norm_mu;
tol = 1e-2;
%%% CHIARA CODE
if ((abs(check(1))>tol)||(abs(check(2))>tol)||(abs(check(3))>tol)) & ...
        ((abs(checki(1))>tol)||(abs(checki(2))>tol)||(abs(checki(3))>tol))
    disp('ERROR: check fault !!!!!!');
    stop
else
    disp('CORRECT FAULT !');
end

%%%MODIFIED BY ILARIO
% err = abs(dot(normale,norm_mu));
% tol = 1;
% if (err > tol)
%     disp('ERROR: check fault !!!!!!');
%     stop
% else
%     disp('CORRECT FAULT !');
% end

% ====================================================================== %
%    3D rotation & translation to map the fault onto the X-Y plane     %
% ====================================================================== %
% disp([' '])
disp(['*** Reading and interpolating data about slip distribution***'])
[v1,v2] = rotrasl3D(STRIKE,DIP,L,W);
grid_X = zeros(num_nodes,1);
grid_Y = zeros(num_nodes,1);
for i=1:num_nodes
    P = [grid_x(i) grid_y(i) grid_z(i)];
    [Pp] = fault_points(v1,v2,P);
    grid_X(i) = Pp(1);
    grid_Y(i) = Pp(2);
end
% hypocenter coordinates
HY = [XC YC ZC];
[HYp] =  fault_points(v1,v2,HY);
HY_X = HYp(1);
HY_Y = HYp(2);

% centroid (Xg,Yg) of each triangle
Xg = zeros(num_ctria,1);
Yg = zeros(num_ctria,1);
X1 = grid_X(ctria_1(1:ctria)')';
Y1 = grid_Y(ctria_1(1:ctria)')';
X2 = grid_X(ctria_2(1:ctria)')';
Y2 = grid_Y(ctria_2(1:ctria)')';
X3 = grid_X(ctria_3(1:ctria)')';
Y3 = grid_Y(ctria_3(1:ctria)')';
Xg = ( X1+X2+X3 )./3;
Yg = ( Y1+Y2+Y3 )./3;

% translation of quantity (t1,t2)
% ORIGIN 0 = zero strike and uppermost (w.r.t. ground surface) point !
t1 = min([X1,X2,X3]); % t1 = min(Xg);
t2 = max([Y1,Y2,Y3]); % t2 = max(Yg);
%     t1 = min(Xg);
%     t2 = max(Yg);

Xg = Xg-t1;
Yg = Yg-t2;
X1 = X1-t1; X2 = X2-t1; X3 = X3-t1;
Y1 = Y1-t2; Y2 = Y2-t2; Y3 = Y3-t2;
HY_X = HY_X -t1;
HY_Y = HY_Y -t2;


if (SCEN.ID == 0)   %iscenario
    Leff = L;
    Weff = W;
    L0 = 0;
    W0 = 0;
else
    L0 = L0_scen_norm*L;
    W0 = W0_scen_norm*W;
    Leff = L_scen;
    Weff = W_scen;
end

xf_scen = [L0,L0+Leff,L0+Leff,L0,L0];
yf_scen = [-W0,-W0,-W0-Weff,-W0-Weff,-W0];
[in on] = inpolygon(Xg,Yg,xf_scen,yf_scen);

ctria_scen = find(in == 1);

[inhy onhy] = inpolygon(HY_X,HY_Y,xf_scen,yf_scen);
% if inhy == 0
%     disp('CHECK HYPO LOCATION !');
%     stop
% end

figure(1000)
plot(X1,Y1,'.k');
hold on
plot(X2,Y2,'.k');
hold on
plot(X3,Y3,'.k');
hold on
plot(Xg,Yg,'.k');
hold on
plot(xf_scen,yf_scen,'-r','linewidth',2);
hold on
plot(Xg(ctria_scen),Yg(ctria_scen),'^g','markersize',6);
hold on
plot(HY_X,HY_Y,'pr','markersize',15);


% ====================================================================== %
%   INTERPOLATION: mapping slip pattern from SRCMOD database (M. Mai)    %
% ====================================================================== %

slip_scrmod1 = load(nomefiler2);
%     if type_orient==2
%         slip_scrmod1 = fliplr(slip_scrmod);
%     elseif type_orient==3
%         slip_scrmod1 = flipud(slip_scrmod);
%     elseif type_orient==4
%         slip_scrmod1 = flipud(fliplr(slip_scrmod));
%     end

dim_scrmod1 = size(slip_scrmod1);
Nstr1 = dim_scrmod1(2);
Ndip1 = dim_scrmod1(1);
dxm_eff = Leff/Nstr1;
dym_eff = Weff/Ndip1;

% plot of original slip distribution from source inversion study
xme = [L0:dxm_eff:L0+Leff-dxm_eff]+dxm_eff/2;
yme = [W0:dym_eff:W0+Weff-dym_eff]+dym_eff/2;

u_mean_orig = mean(mean(slip_scrmod1));
u_max_orig = max(max(slip_scrmod1));

figure(10)
set(gca,'fontsize',20);
imagesc([1:Nstr1]/Nstr1-1/2/Nstr1,[1:Ndip1]/Ndip1-1/2/Ndip1,slip_scrmod1/max(max(slip_scrmod1))); %shading flat; hold on;
hold on
view(2)
colormap(jet)
hold on
xlabel('Along Strike (km)');
ylabel('Down Dip (km)');
colorbar('SouthOutside','fontsize',20)
xlim([0  1]);
ylim([0 1]);


fig_name = [char(path_f),[eqk_str,'_SlipOrig']];
print('-f10', '-djpeg','-r100',fig_name);



if (SCEN.ID == 0)   %iscenario
    Nstr = (Nstr1);
    Ndip = (Ndip1);
    slip_scrmod = slip_scrmod1;
else
    dx_plus = dxm_eff;
    dy_plus = dym_eff;
    Ndip_up = round(W0/dy_plus)+1;
    Ndip_bottom = round((W-Weff-W0)/dy_plus)+1;
    Ndip_plus = Ndip_up+Ndip_bottom;
    Nstr_left = round(L0/dx_plus)+1;
    Nstr_right = round((L-Leff-L0)/dx_plus)+1;
    Nstr_plus = Nstr_left+Nstr_right;
    
    Nstr = (Nstr1+Nstr_plus);
    Ndip = (Ndip1+Ndip_plus);
    slip_scrmod = zeros(Ndip,Nstr);
    slip_scrmod(Ndip_up+1:Ndip1+Ndip_up,Nstr_left+1:Nstr1+Nstr_left) = ....
        slip_scrmod1(1:Ndip1,1:Nstr1);
end

%-------------------------------------------------------------------------
% Interpolation of slip distribution
%----------------------------------------------------------------------

xm = zeros(Nstr*Ndip,1);
ym = zeros(Nstr*Ndip,1);
zm = zeros(Nstr*Ndip,1);
rakem = zeros(Nstr*Ndip,1);
cont = 1;

if (SCEN.ID == 0)
    for j=1:Ndip
        for i=1:Nstr
            xm(cont) = dxm_eff*(i-1) + dxm_eff/2;
            ym(cont) = -dym_eff*(j-1) - dym_eff/2;
            zm(cont) = slip_scrmod(j,i);
            cont = cont +1;
        end
    end
else
    for j=1:Ndip
        for i=1:Nstr
            xm(cont) = dxm_eff*(i-1)+(L0+dxm_eff/2-Nstr_left*dxm_eff);
            ym(cont) = -dym_eff*(j-1) -(W0+dym_eff/2-Ndip_up*dym_eff);
            zm(cont) = slip_scrmod(j,i);
            cont = cont +1;
        end
    end
    
end

if (SCEN.ID == 0)
    lon0 = 0;
    lonm = Leff;
    lat0 = 0;
    latm = -Weff;
else
    lon0 = min(xm);
    lonm = max(xm);
    lat0 = max(ym);
    latm = min(ym);
end
dlon = 100;
dlat = 100;
lon = [lon0:dlon:lonm];
lat = [lat0:-dlat:latm];

% $$$$$$$$$$$$$$$ %
if mod(length(lon),2) ~= 0
    lon = [lon0:dlon:lonm+dlon];
end
if mod(length(lat),2) ~= 0
    lat = [lat0:-dlat:latm-dlat];
end
% $$$$$$$$$$$$$$$ %

[loni,lati] = meshgrid(lon,lat);
elevi = griddata(xm,ym,zm,loni,lati,'nearest');

Sd = zeros(num_ctria,1);
Sd = interp2(loni,lati,elevi,Xg,Yg,'nearest');
for i = 1 : length(Sd)
    if(isnan(Sd(i))); Sd(i) = 0; end;
end



figure(2);
set(gca,'fontsize',20);
surf(loni./1000,lati./1000,elevi);
hold on
view(2)
shading flat;%shading interp ;
colormap(jet)
hold on
plot3(HY_X./1000,HY_Y./1000,[1000],'pr','linewidth',2.0,'markersize',28,'markerfacecolor',[1 0 0]);hold on;
xlabel('Along Strike (km)');
ylabel('Down Dip (km)');
caxis([min(min(elevi)) max(max(elevi)) ])
%     colorbar
colorbar('SouthOutside','fontsize',20)
xlim([0 L/1000]);
ylim([-W/1000 0]);
box on

fig_name = [char(path_f),['E',num2str(scen_num),'_S',num2str(slip_scen),'_H',...
    num2str(Hypocenter),'_Slip']];
print('-f2', '-djpeg','-r100',fig_name);

% $$$$$$$$$$$$$$$ %
if tag_crsp_rake == 1
    disp(' ')
    disp('*** Correlated random perturbation of rake angle ***')
    [rake_crsp] = crsp(L,W,ax,ay,dlon,dlat,elevi,mean_rake,cov_rake,eta_rake,type_psd_rake);
    
    figure(3)
    set(gca,'fontsize',20);
    surf(loni./1000,lati./1000,rake_crsp)
    view(2)
    shading flat
    xlabel('Along Strike (km)');
    %        ylabel('Down Dip (km)');
    %        title('\it{\lambda_f}');
    colorbar('southoutside'); colormap(jet);
    xlim([0 max(max(loni./1000))]);
    ylim([min(min(lati./1000)) 0]);
    hold on
    
    fig_name = [char(path_f),['E',num2str(scen_num),'_S',num2str(slip_scen),'_H',...
        num2str(Hypocenter),'_Rake_crsp']];
    print('-f3', '-djpeg','-r100',fig_name);
    
    Rd = interp2(loni,lati,rake_crsp,Xg,Yg);
    
end

if tag_crsp_rt == 1
    disp(' ')
    disp('*** Correlated random perturbation of rise time ***')
    [tau_crsp] = crsp(L,W,ax,ay,dlon,dlat,elevi,mean_rt,cov_rt,eta_rt,type_psd_rt);
    
    
    figure(4)
    set(gca,'fontsize',20);
    surf(loni./1000,lati./1000,tau_crsp)
    view(2)
    shading flat
    xlabel('Along Strike (km)');
    %        ylabel('Down Dip (km)');
    %        title('\it{\tau_R}');
    colorbar('southoutside'); colormap(jet);
    xlim([0 max(max(loni./1000))]);
    ylim([min(min(lati./1000)) 0]);
    hold on
    
    fig_name = [char(path_f),['E',num2str(scen_num),'_S',num2str(slip_scen),'_H',...
        num2str(Hypocenter),'_RiseTime_crsp']];
    print('-f4', '-djpeg','-r100',fig_name);
    
    RTd = interp2(loni,lati,tau_crsp,Xg,Yg);
    
end

if tag_crsp_vr == 1
    disp([' '])
    disp(['*** Correlated random perturbation of rupture velocity ***'])
    [vr_crsp] = crsp(L,W,ax,ay,dlon,dlat,elevi,mean_vr,cov_vr,eta_vr,type_psd_vr);
    
    figure(5)
    set(gca,'fontsize',20);
    surf(loni./1000,lati./1000,vr_crsp)
    view(2)
    shading flat
    xlabel('Along Strike (km)');
    %        ylabel('Down Dip (km)');
    %        title('\it{V_R}');
    colorbar
    xlim([0 max(max(loni./1000))]);
    ylim([min(min(lati./1000)) 0]);
    hold on
    
    fig_name = [char(path_f),['E',num2str(scen_num),'_S',num2str(slip_scen),'_H',...
        num2str(Hypocenter),'_Vr_crsp']];
    print('-f4', '-djpeg','-r100',fig_name);
    
    VRd = interp2(loni,lati,vr_crsp,Xg,Yg);
    
end

% $$$$$$$$$$$$$$$ %

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% FAULT COMPUTATIONS
tini=cputime;
disp([' '])
disp(['BEGIN - Computing area and so on... '])

x1 = grid_x(ctria_1(1:ctria)')';
y1 = grid_y(ctria_1(1:ctria)')';
z1 = grid_z(ctria_1(1:ctria)')';

x2 = grid_x(ctria_2(1:ctria)')';
y2 = grid_y(ctria_2(1:ctria)')';
z2 = grid_z(ctria_2(1:ctria)')';

x3 = grid_x(ctria_3(1:ctria)')';
y3 = grid_y(ctria_3(1:ctria)')';
z3 = grid_z(ctria_3(1:ctria)')';

one = ones(1,ctria);

area_ctria = zeros(1,ctria);
dist_bar_hypo = zeros(1,ctria);

x_average = zeros(1,ctria);
y_average = zeros(1,ctria);
z_average = zeros(1,ctria);
mu = zeros(1,ctria);
delta_u = zeros(1,ctria);

% http://mathworld.wolfram.com/TriangleArea.html
for i = 1:ctria
    A1 = [y1(i) z1(i) one(i);
        y2(i) z2(i) one(i);
        y3(i) z3(i) one(i);];
    
    A2 = [z1(i) x1(i) one(i);
        z2(i) x2(i) one(i);
        z3(i) x3(i) one(i);];
    
    A3 = [x1(i) y1(i) one(i);
        x2(i) y2(i) one(i);
        x3(i) y3(i) one(i);];
    
    area_ctria(i) = 1/2.* sqrt( (det(A1)).^2 + (det(A2)).^2 + (det(A3)).^2 );
    x_average(i) = (x1(i) + x2(i) + x3(i))/3;
    y_average(i) = (y1(i) + y2(i) + y3(i))/3;
    z_average(i) = (z1(i) + z2(i) + z3(i))/3;
    dist_bar_hypo(i) = sqrt((x_average(i)-XC)^2 + (y_average(i)-YC)^2 ...
        + (z_average(i)-ZC)^2);
    
    
    delta_u(i) = Sd(i); % from srcmod database
    
    for j = 1:num_tag_ctria
        if  ctria_tag(i)==j
            mu(i) = rho(j)*Vs(j)^2;
        end
        %ONLY FOR ISTANBUL NH
        if (tagcase == 16)
            if  (ctria_tag(i) == 1)
                Depth_nh = FO_z - z_average(i);
                if(Depth_nh < 0 ); Depth_nh = 0; end
                Vs_nh =1350 + 23.33*(Depth_nh)^(0.5);
                rho_nh = 2100 + 5.69*(Depth_nh)^(0.5);
                mu(i) = rho_nh*Vs_nh^2;
            end
        end
    end
end

% scale factor to match seismic moment value from Hanks &
% Kanamori (1979) relationship
m0 = area_ctria.*mu.*delta_u;
scale_factor = Magnitude0/sum(m0);

if iscale == 1
    m0 = m0.*scale_factor;
    % slip pattern after scale factor
    elevi = elevi.*scale_factor;
else
    m0 = m0;
    elevi = elevi;
end
%     sum(m0);

delta_u_fin = m0./(area_ctria.*mu);

triangles = [ctria_1 ctria_2 ctria_3];
figure(80)
values_deltau_pnt = zeros(length(grid_x),1);
for i = 1 : ctria
    ind1 = triangles(i,1);
    ind2 = triangles(i,2);
    ind3 = triangles(i,3);
    values_deltau_pnt(ind1) = delta_u_fin(i);
    values_deltau_pnt(ind2) = delta_u_fin(i);
    values_deltau_pnt(ind3) = delta_u_fin(i);
end


trisurf(triangles,grid_x,grid_y,grid_z,values_deltau_pnt)

fprintf('Magnitude Target %6.2f \n',  2/3*log10(Magnitude0) -6);
fprintf('Scale factor  %6.2f \n',  scale_factor);
fprintf('Computed Magnitude %6.2f \n',  2/3*log10(sum(m0)) -6)

figure(8);
set(gca,'fontsize',20);
surf(loni./1000,lati./1000,elevi); %shading flat; hold on;
view(2)
hold on
%     plot3(Xg/1000,Yg/1000,delta_u_fin,'.r','markersize',20);
shading flat;
colormap(jet)
hold on
plot3(HY_X./1000,HY_Y./1000,[1000],'pr','linewidth',2.0,'markersize',28,'markerfacecolor',[1 0 0]);hold on;
xlabel('Along Strike (km)');
ylabel('Down Dip (km)');
caxis([min(min(elevi)) max(max(elevi)) ])
colorbar('SouthOutside','fontsize',20)
xlim([0 max(max(loni./1000))]);
ylim([min(min(lati./1000)) 0]);
box on


fig_name = [char(path_f),[eqk_str,'_Slip']];
print('-f8','-djpeg','-r100',fig_name);


disp(['END - Computing area and so on... ',char(num2str(cputime-tini)),' sec.'])
disp([' '])


% ======================================================================= %
%  ================ writing output files ================================ %
% ======================================================================= %

%     if grid>0
%             fprintf(fid,'%i  %+13.7e  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
%     end
%     fprintf(fid,'End Coordinates\n');
%     fprintf(fid,'Elements\n');
%     fprintf(fid,'# element node_1 node_2 node_3 material_number\n');
%     if ctria>0
%             fprintf(fid,'%i    %i  %i  %i    %i\n',...
%                         [ctria_id(1:ctria)';...
%                         ctria_1(1:ctria)'; ctria_2(1:ctria)'; ctria_3(1:ctria)'; ctria_tag(1:ctria)']);
%     end
%     fprintf(fid,'End Elements\n');
%
%      fclose(fid);


% fid = fopen(nomefilew2,'w');
%     fprintf(fid,'MONITOR FILES\n');
%
%     if grid>0
%             fprintf(fid,'MONITOR  %+13.7e  %+13.7e  %+13.7e\n',[ grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
%     end
% fclose(fid);


% fid = fopen(nomefilew3,'w');
%     fprintf(fid,'%i \n',num_nodes);
%
%     if grid>0
%             fprintf(fid,'%i  %+13.7e  %+13.7e  %+13.7e\n',[grid_id(1:grid)'; grid_x(1:grid)'; grid_y(1:grid)'; grid_z(1:grid)']);
%     end
% fclose(fid);

% $$$$$$$$$$$$$$$ %

% ==================================================== %
% random perturbations of rupture velocity
% ==================================================== %
if tag_crsp_vr == 0
    vel_rup = [];
    for i = 1:ctria
        if (ctria_tag(i) == 1 && strcmp(utmzone,'35T'))
            Depth_nh = FO_z - z_average(i);
            if(Depth_nh < 0 ); Depth_nh = 0; end
            Vs_nh =1350 + 23.33*(Depth_nh)^(0.5);
            vel_rup(i) = 0.85*Vs_nh;
        else
            vel_rup(i) = 0.85*Vs(ctria_tag(i));
        end
        
        
    end
    
    vel_rnd_mod = vel_rup;
    %vel_rnd_mod = velocity.*one;
    
elseif tag_crsp_vr == 1
    for i = 1:ctria
        vel_rnd(i) = VRd(i)
        vel_rnd_mod(i) = vel_rnd(i);
        if vel_rnd_mod(i) > vel_threshold
            vel_rnd_mod(i) = vel_threshold;
        end
    end
end

figure(6)
set(gca,'fontsize',16);
plot(vel_rnd_mod,'k+-','Linewidth',3); hold on;
ylabel('Rupture Velocity (m/s)','fontsize',16);
xlabel('# ctria','fontsize',16);


% ==================================================== %
% rake parameter
% ==================================================== %
if tag_crsp_rake == 0
    % Slip vector (uses slip_vec.m)
    slip = slip_vec(STRIKE,RAKE,DIP);
    slip_x = slip_x_mu.*one;
    slip_y = slip_y_mu.*one;
    slip_z = slip_z_mu.*one;
    % Fault normal vector (uses norm_vec.m)
    norm_x = norm_x_mu.*one;
    norm_y = norm_y_mu.*one;
    norm_z = norm_z_mu.*one;
    % slip norm vector transformation
    [Sp] = fault_points(v1,v2,slip);
    S_X = Sp(1);
    S_X = S_X.*one;
    S_Y = Sp(2);
    S_Y = S_Y.*one;
    
elseif tag_crsp_rake == 1
    for i = 1:ctria
        % Slip vector (uses slip_vec.m)
        slip = slip_vec(STRIKE,Rd(i),DIP);
        slip_x(i) = slip(1);slip_y(i)= slip(2);slip_z(i) = slip(3);
        % Fault normal vector (uses norm_vec.m)
        normr = norm_vec(STRIKE,Rd(i),DIP);
        norm_x(i) = normr(1);norm_y(i) = normr(2);norm_z(i) = normr(3);
        % slip norm vector transformation
        [Sp] = fault_points(v1,v2,slip);
        S_X(i) = Sp(1);
        S_Y(i) = Sp(2);
    end
end
figure(7)
set(gca,'fontsize',20);
quiver(Xg./1000,Yg./1000,S_X,S_Y,'-r')
xlabel('Along Strike (km)','fontsize',20);
ylabel('Along Dip (km)','fontsize',20);
xlim([0 (L+0.01*L)/1000]);
ylim([(-W-0.01*W)/1000 100/1000]);

fig_name = [char(path_f),'E',num2str(scen_num),'_S',num2str(slip_scen),'_H',...
    num2str(Hypocenter),'_RakeVec'];
print('-f7','-djpeg','-r100',fig_name);

% ==================================================== %
% rise time parameter
% ==================================================== %
tag_nf = nfunc.*one;
if tag_crsp_rt == 1
    tau_all = RTd.*one;
else
    tau_all = tau.*one;
end
% $$$$$$$$$$$$$$$ %


fid = fopen(nomefilew4,'w');

% if ctria>0
%     fprintf(fid,'SISM  %d  2     %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e    %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e\n',...
%         [tag_nf;XC.*one; YC.*one; ZC.*one;...
%         x1; y1; z1;...
%         x2; y2; z2;...
%         x3; y3; z3;...
%         slip_x; slip_y; slip_z;...
%         norm_x; norm_y; norm_z;...
%         vel_rnd_mod; m0;tau_all]);
% end

for i = 1 : ctria
    
    if(m0(i) > 0)
        fprintf(fid,'SISM  %d  2  %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e  %+13.7e  %+13.7e  %+13.7e    %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e \n',...
            tag_nf(i), XC(1) , YC(1) , ZC(1) , x1(i), y1(i), z1(i), x2(i), y2(i), z2(i), x3(i), y3(i), z3(i), ...
            slip_x(i), slip_y(i), slip_z(i),  norm_x(i), norm_y(i), norm_z(i), dist_bar_hypo(i)/vel_rnd_mod(i), m0(i), tau_all(i));
        
    end
end


fclose(fid);

disp(['scale factor for seismic moment: ', num2str(scale_factor)])

% scalar seismic moment after scale factor
M0_final = sum(m0);
u_mean = mean(delta_u_fin);
u_max = max(delta_u_fin);



fid5 = fopen(nomefilew5,'w');
fprintf(fid5,'%s\n',['*********************** SCENARIO PARAMETERS *********************** ']);
fprintf(fid5,'%s\n',[utmzone,' : ',loc_str]);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['**** FAULT GEOMETRY ****']);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['Fault ID = ',fault_str]);
% fprintf(fid5,'%s\n',['Maximum moment magnitude Mw_max = ',num2str(Mw_max)]);
fprintf(fid5,'%s\n',['Fault Length = ',num2str(L./1000),' km']);
fprintf(fid5,'%s\n',['Fault Width = ',num2str(W./1000),' km']);
fprintf(fid5,'%s\n',['Fault Strike = ',num2str(STRIKE),'°']);
fprintf(fid5,'%s\n',['Fault Dip = ',num2str(DIP),'°']);
fprintf(fid5,'%s\n',['Fault Origin (FO): X_FO = ',num2str(FO_x),' m - Y_FO = ',num2str(FO_y),' m']);
fprintf(fid5,'%s\n',['Latitude of FO = ',num2str(Lat_FO),'° - Longiutde of FO = ',num2str(Lon_FO),'°'....
    ' - utmzone ',utmzone]);
fprintf(fid5,'%s\n',['Depth of FO = ',num2str(FO_z/1000),' km']);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['**** EQK SCENARIO   ****']);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['Scenario ID = ',eqk_str]);
fprintf(fid5,'%s\n',['Moment Magnitude = ',num2str(Mw)]);
fprintf(fid5,'%s\n',['Slip model #: ',slip_str,' : ', nomer2]);
fprintf(fid5,'%s\n',['Rupture Length L = ',num2str(Leff./1000),' km']);
fprintf(fid5,'%s\n',['Rupture Width W = ',num2str(Weff./1000),' km']);
fprintf(fid5,'%s\n',['Normalized distance L0 of rupture scenario along strike from FO (w.r.t Lmax) = ',num2str(L0./L)]);
fprintf(fid5,'%s\n',['Normalized distance W0 of rupture scenario down dip from FO (w.r.t Wmax) = ',num2str(W0./W)]);
fprintf(fid5,'%s\n',['Normalized distance Hypo_X of hypocenter along strike from FO (w.r.t Lmax) = ',num2str(HY_X./L)]);
fprintf(fid5,'%s\n',['Normalized distance Hypo_Y of hypocenter down dip from FO (w.r.t Wmax) = ',num2str(abs(HY_Y)./W)]);
fprintf(fid5,'%s\n',['X_HY = ',num2str(XC),' m - Y_HY = ',num2str(YC),' m']);
fprintf(fid5,'%s\n',['Latitude of HYPO = ',num2str(Lat_HY),'° - Longitude of HYPO = ',num2str(Lon_HY),'°',...
    ' - utmzone ',utmzone]);
fprintf(fid5,'%s\n',['Focal Depth = ',num2str(ZC/1000),' km']);
fprintf(fid5,'%s\n',['Rake angle = ',num2str(RAKE),' deg']);
fprintf(fid5,'%s\n',['Rupture Velocity = ',num2str(velocity/1000),' km/s']);
fprintf(fid5,'%s\n',['Rise time = ',num2str(tau),' s']);
fprintf(fid5,'%s\n',['Mean slip value: ',num2str(u_mean),' m']);
fprintf(fid5,'%s\n',['Maximum slip value: ',num2str(u_max),' m']);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['****   SLIP MODEL   ****']);
fprintf(fid5,'%s\n',['************************']);
fprintf(fid5,'%s\n',['Slip model ID: ', slip_str]);
fprintf(fid5,'%s\n',['Slip model reference: ',nomer2]);
fprintf(fid5,'%s\n',['Mean slip value (original): ',num2str(u_mean_orig),' m']);
fprintf(fid5,'%s\n',['Maximum slip value (original): ',num2str(u_max_orig),' m']);
% $$$$$$$$$$$$$$$ %
fprintf(fid5,'%s\n',['*************************************************']);
fprintf(fid5,'%s\n',['****   CORRELATED RANDOM SOURCE PARAMETERS   ****']);
fprintf(fid5,'%s\n',['*************************************************']);
if tag_crsp_rake == 0
    out = 'NO';
else
    out = 'YES';
end
fprintf(fid5,'%s\n',['Random perturbation of rake angle : ',out]);
fprintf(fid5,'%s\n',['(Mean) Rake angle : ',num2str(RAKE),'°']);
if tag_crsp_rake == 1
    fprintf(fid5,'%s\n',['Rake - coefficient of variation: ',num2str(cov_rake)]);
    fprintf(fid5,'%s\n',['Rake - correlation length in along-strike direction: ',num2str(ax),' m']);
    fprintf(fid5,'%s\n',['Rake - correlation length in down-dip direction: ',num2str(ay),' m']);
    fprintf(fid5,'%s\n',['Rake - correlation coefficient with slip distribution: ',num2str(eta_rake)]);
end
if tag_crsp_rt == 0
    out = 'NO';
else
    out = 'YES';
end
fprintf(fid5,'%s\n',['Random perturbation of rise time : ',out]);
fprintf(fid5,'%s\n',['(Mean) Rise Time : ',num2str(tau),' s']);
if tag_crsp_rt == 1
    fprintf(fid5,'%s\n',['Rise Time - coefficient of variation: ',num2str(cov_rt)]);
    fprintf(fid5,'%s\n',['Rise Time - correlation length in along-strike direction: ',num2str(ax),' m']);
    fprintf(fid5,'%s\n',['Rise Time - correlation length in down-dip direction: ',num2str(ay),' m']);
    fprintf(fid5,'%s\n',['Rise Time - correlation coefficient with slip distribution: ',num2str(eta_rt)]);
end
if tag_crsp_vr == 0
    out = 'NO';
else
    out = 'YES';
end
fprintf(fid5,'%s\n',['Random perturbation of rupture velocity : ',out]);
fprintf(fid5,'%s\n',['(Mean) Rupt. Vel. : ',num2str(velocity/1000),' km/s']);
if tag_crsp_vr == 1
    fprintf(fid5,'%s\n',['Rupt. Vel. - coefficient of variation: ',num2str(cov_vr)]);
    fprintf(fid5,'%s\n',['Rupt. Vel. - correlation length in along-strike direction: ',num2str(ax),' m']);
    fprintf(fid5,'%s\n',['Rupt. Vel. - correlation length in down-dip direction: ',num2str(ay),' m']);
    fprintf(fid5,'%s\n',['Rupt. Vel. - correlation coefficient with slip distribution: ',num2str(eta_vr)]);
    fprintf(fid5,'%s\n',['Rupt. Vel. - velocity threshold = ',num2str(vel_threshold/1000),' km/s']);
end
fclose(fid5);
% $$$$$$$$$$$$$$$ %




fid6 = fopen(nomefilew6,'a');
fprintf(fid6,'%s\n',['scen_ID loc_ID fault_ID Mw slip_ID L W L0 W0 Hypo_X Hypo_Y',...
    'Hypo_Lon Hypo_Lat Hypo_Depth Vr tau rake']);
fprintf(fid6,'%s\n',[eqk_str,' ', loc_str,' ', fault_str,' ', num2str(Mw),' ', ...
    slip_str,' ', num2str(L_scen./1000),' ',num2str(W_scen./1000),' ',...
    num2str(L0./L),' ', num2str(W0./W),' ',...
    num2str(HY_X./L), ' ',num2str(abs(HY_Y)./W),' ',...
    num2str(Lat_HY),' ',num2str(Lon_HY),' ',num2str(ZC/1000),' ',...
    num2str(velocity/1000),' ',num2str(tau),' ',num2str(RAKE)]);
fclose(fid6);

fid_mate_out = fopen(nomefile_mate,'a');

% if ctria>0
%     fprintf(fid_mate_out,'SISM  %d  2     %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e    %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e\n',...
%         [tag_nf;XC.*one; YC.*one; ZC.*one;...
%         x1; y1; z1;...
%         x2; y2; z2;...
%         x3; y3; z3;...
%         slip_x; slip_y; slip_z;...
%         norm_x; norm_y; norm_z;...
%         vel_rnd_mod; m0;tau_all]);
% end

for i = 1 : ctria
    
    if(m0(i) > 0)
        fprintf(fid,'SISM  %d  2  %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e    %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e      %+13.7e  %+13.7e  %+13.7e     %+13.7e    %+13.7e    %+13.7e \n',...
            tag_nf(i), XC(1) , YC(1) , ZC(1) , x1(i), y1(i), z1(i), x2(i), y2(i), z2(i), x3(i), y3(i), z3(i), ...
            slip_x(i), slip_y(i), slip_z(i),  norm_x(i), norm_y(i), norm_z(i), dist_bar_hypo(i)/vel_rnd_mod(i), m0(i), tau_all(i));
        
    end
end



fclose(fid_mate_out);

% copyfile(nomefile_mate_out,nomefile_mate)
% delete(nomefile_mate_out);

disp('---------------------------------------------------------------')
