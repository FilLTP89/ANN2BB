
clear all
close all
clc
warning off


path = ['C:\Users\chiara.smerzini\Documents\Thessaloniki2014-2015\',...
    'Work\STREST\Thessaloniki_Simulation\'];  

casename_soil = 'Volvi1978_SD3_k2_ALL_final_improved_NL';
casename_rock = 'Volvi1978_SD3_k2_rock'; 

% casename_soil = 'E00560_Anthemoutas_composite_SD3_k2_ALL_improved_NL'; 
% casename_rock = 'E00560_Anthemoutas_composite_SD3_k2_rock'; 

subpath_speed_soil = ['SPEED\run\',casename_soil,'\'];
subpath_speed_rock = ['SPEED\run\',casename_rock,'\'];
% subpath_map_soil = ['PGMAP\']; 
% subpath_map_rock = ['PGMAP\'];
% subpath_map_soil = ['PGMAP_1.5Hz\']; 
% subpath_map_rock = ['PGMAP_1.5Hz\'];
subpath_map_soil = ['PGMAP_2Hz\']; 
subpath_map_rock = ['PGMAP_2Hz\'];
nn = 16455; 

type = 'PGV'; 
% type = 'PSA_Tn2'; 

outfile = [path,subpath_speed_soil,subpath_map_soil,type,'_ratio.shp']; 

nomefile_soil = [path,subpath_speed_soil,subpath_map_soil,type,'.shp']; 
lett_soil = shaperead(nomefile_soil);


nomefile_rock = [path,subpath_speed_rock,subpath_map_rock,type,'.shp']; 
lett_rock = shaperead(nomefile_rock);


% check on coordinates
for ii = 1:nn
    Xs(ii) = lett_soil(ii).X; 
    Ys(ii) = lett_soil(ii).Y; 
    Xr(ii) = lett_rock(ii).X;
    Yr(ii) = lett_rock(ii).Y; 
end

Xcheck = Xs-Xr;
Ycheck = Ys-Yr; 
if (abs(sum(Xcheck))) > 1e-3 || (abs(sum(Ycheck))) > 1e-3
    disp('error in rock/soil coor'); 
    stop
end



[PGMAP(1:nn).Geometry] = deal('Point');
for ii = 1:nn
    [PGMAP(ii).X] = Xs(ii); 
    [PGMAP(ii).Y] = Ys(ii);
    if type(1:3) == 'PGV'
        [PGMAP(ii).Rx] = lett_soil(ii).PGVx/lett_rock(ii).PGVx;
        [PGMAP(ii).Ry] = lett_soil(ii).PGVy/lett_rock(ii).PGVy;
        [PGMAP(ii).Rz] = lett_soil(ii).PGVz/lett_rock(ii).PGVz;
        [PGMAP(ii).Rp] = lett_soil(ii).PGVp/lett_rock(ii).PGVp;
        [PGMAP(ii).Ro] = lett_soil(ii).PGVo/lett_rock(ii).PGVo;
        [PGMAP(ii).Rgmh] = lett_soil(ii).PGVgmh/lett_rock(ii).PGVgmh;
    elseif type(1:3) == 'PGD'
        [PGMAP(ii).Rx] = lett_soil(ii).PGDx/lett_rock(ii).PGDx;
        [PGMAP(ii).Ry] = lett_soil(ii).PGDy/lett_rock(ii).PGDy;
        [PGMAP(ii).Rz] = lett_soil(ii).PGDz/lett_rock(ii).PGDz;
        [PGMAP(ii).Rp] = lett_soil(ii).PGDp/lett_rock(ii).PGDp;
        [PGMAP(ii).Ro] = lett_soil(ii).PGDo/lett_rock(ii).PGDo;
        [PGMAP(ii).Rgmh] = lett_soil(ii).PGDgmh/lett_rock(ii).PGDgmh;  
    elseif type(1:3) == 'SD_'
        [PGMAP(ii).Rx] = lett_soil(ii).SDx/lett_rock(ii).SDx;
        [PGMAP(ii).Ry] = lett_soil(ii).SDy/lett_rock(ii).SDy;
        [PGMAP(ii).Rz] = lett_soil(ii).SDz/lett_rock(ii).SDz;
        [PGMAP(ii).Rp] = lett_soil(ii).SDp/lett_rock(ii).SDp;
        [PGMAP(ii).Ro] = lett_soil(ii).SDo/lett_rock(ii).SDo;
        [PGMAP(ii).Rgmh] = lett_soil(ii).SDgmh/lett_rock(ii).SDgmh;  
    elseif type(1:3) == 'PSA'
        [PGMAP(ii).Rx] = lett_soil(ii).PSAx/lett_rock(ii).PSAx;
        [PGMAP(ii).Ry] = lett_soil(ii).PSAy/lett_rock(ii).PSAy;
        [PGMAP(ii).Rz] = lett_soil(ii).PSAz/lett_rock(ii).PSAz;
        [PGMAP(ii).Rp] = lett_soil(ii).PSAp/lett_rock(ii).PSAp;
        [PGMAP(ii).Ro] = lett_soil(ii).PSAo/lett_rock(ii).PSAo;
        [PGMAP(ii).Rgmh] = lett_soil(ii).PSAgmh/lett_rock(ii).PSAgmh;  
    end
end


shapewrite(PGMAP,outfile); 





