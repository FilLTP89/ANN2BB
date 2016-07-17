ccc

z = 2.25:.9:5.85;
x = 0.75.*ones(numel(z),1)';
y = 0.58.*ones(numel(z),1)';

v = 120;
f = 380;
ti = 1/f;
tf = 0.2;
t = ti+(z-z(1))./v;

TRAVERSE = 1;

% return
%%
fid = fopen('input.spec', 'w+');

fprintf(fid,'# -*- mode: perl -*-\n');
fprintf(fid,'run_name = "Sm_h";\n');
fprintf(fid,'\n');

fprintf(fid,'# duration of the run\n');
fprintf(fid,'sim_time = %11.4f;\n',tf);
fprintf(fid,'mesh_file = "mesh4spec"; # input mesh file\n');
fprintf(fid,'mat_file = "material.input";\n');
fprintf(fid,'dim=3;\n');
fprintf(fid,'\n');

fprintf(fid,'snapshots {\n');
fprintf(fid,'    save_snap = true;\n');
fprintf(fid,'    snap_interval = 0.004;\n');
fprintf(fid,'};\n');
fprintf(fid,'\n');

fprintf(fid,'# Description des capteurs\n');
fprintf(fid,'save_traces = true;\n');
fprintf(fid,'station_file = "captures.dat";\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'# Fichier protection reprise\n');
fprintf(fid,'prorep=false;\n');
fprintf(fid,'prorep_iter=1000;\n');
fprintf(fid,'restart_iter=370;\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

count=0;

for i = 1:numel(t)
    
    fprintf(fid,'# introduce a source %i\n',i);
    fprintf(fid,'source {\n');
    fprintf(fid,'# coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)\n');
    fprintf(fid,'coords = %11.2f. %11.2f. %11.2f.;\n',x(i),y(i),z(i));
    fprintf(fid,'# the numbers before the labels are here to help convert from previous input.spec format\n');
    fprintf(fid,'# Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)\n');
    fprintf(fid,'type = impulse;\n');
    fprintf(fid,'# Direction 0.x,1.y ou 2.z (only for Impulse)\n');
    fprintf(fid,'dir = 0. 1. 0.;\n');
    fprintf(fid,'# Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus\n');
    fprintf(fid,'func = ricker;\n');
    fprintf(fid,'tau = %11.10f;\n',t(i));
    fprintf(fid,'freq = %11.2f;   # source main frequency / cutoff frequency\n',f);
    fprintf(fid,'};\n');
    
    fprintf(fid,'\n');
    
    count = count+1;
    
end

if TRAVERSE
    count=0;
    for i = 1:numel(t)
        
        fprintf(fid,'# introduce a source %i\n',i);
        fprintf(fid,'source {\n');
        fprintf(fid,'# coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)\n');
        fprintf(fid,'coords = %11.2f. %11.2f. %11.2f.;\n',-x(i),y(i),z(i));
        fprintf(fid,'# the numbers before the labels are here to help convert from previous input.spec format\n');
        fprintf(fid,'# Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)\n');
        fprintf(fid,'type = impulse;\n');
        fprintf(fid,'# Direction 0.x,1.y ou 2.z (only for Impulse)\n');
        fprintf(fid,'dir = 0. 1. 0.;\n');
        fprintf(fid,'# Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus\n');
        fprintf(fid,'func = ricker;\n');
        fprintf(fid,'tau = %11.10f;\n',t(i));
        fprintf(fid,'freq = %11.2f;   # source main frequency / cutoff frequency\n',f);
        fprintf(fid,'};\n');
        
        fprintf(fid,'\n');
        
        count = count+1;
        
    end
end
fprintf(fid,'time_scheme {\n');
fprintf(fid,'    accel_scheme = false;  # Acceleration scheme for Newmark\n');
fprintf(fid,'    veloc_scheme = true;   # Velocity scheme for Newmark\n');
fprintf(fid,'    alpha = 0.5;           # alpha (Newmark parameter)\n');
fprintf(fid,'    beta = -0.5;           # beta (Newmark parameter)\n');
fprintf(fid,'    gamma = 1;             # gamma (Newmark parameter)\n');
fprintf(fid,'    courant=0.2;\n');
fprintf(fid,'};\n');
fprintf(fid,'\n');

fprintf(fid,'amortissement {\n');
fprintf(fid,'    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)\n');
fprintf(fid,'    atn_band = 10  0.05;   # attenuation period band\n');
fprintf(fid,'    atn_period = 0.2;      # model period\n');
fprintf(fid,'};\n');

fclose(fid);