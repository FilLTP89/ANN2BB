% ================
% rewrite_monitor_format
% Editor: Ilario Mazzieri & Filippo Gatti
% Politecnico di Milano MOX/DICA
% Ecole Centrale Paris - Laboratoire MSSMat
% Copyright 2015
% NOTES
% rewrite_monitor_format: function to rewrite monitor format from speed
% output
% N.B.: files MONITORXXXXX.D ==> monitorXXXXX.d
% INPUT: OUT_OPT (output options (as written in SPEED.input file))
%        path_start (directory where MONITOR files are stored)
%        path_end (directroy where monitor files will be saved)
% REFERENCES:
% http://speed.mox.polimi.it/DOWNLOAD/DOXYGEN_NEW/html/page2.html
% ================

function rewrite_monitor_format(varargin)
    
    OUT_OPT = varargin{1};
    if isunix
        pslash='/';
    elseif ispc
        pslash='\';
    end
    if nargin==1
        path_start=cd;
        path_end=cd;
    elseif nargin==2
        path_start=varargin{2};
        path_end=cd;
    elseif nargin==3
        path_start=varargin{2};
        path_end=varargin{3};
    end
    path_start=sprintf('%s%s',path_start,pslash);
    path_end=sprintf('%s%s',path_end,pslash);
    
    INFO = load([path_start,'MONITOR.INFO']);
    
    t0 = 0;
    T = INFO(1,1);                             %final time
    dt_s = INFO(2,1);                          %deltat simulation
    ndt_monit = INFO(3,1);
    dt = ndt_monit*dt_s;                       %deltat monitor
    
    MPI_num_proc = INFO(4,1);                  %number of mpi proc
    MPI_mnt_id = INFO([5:5+MPI_num_proc-1],1); %id mpi for monitors
    
    
    % DISPLACEMENT
    
    if(OUT_OPT(1) == 1)
        
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Displacement...']);
            filename1 = 'MONITORXXXXX.D';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.D'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.D'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.D'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.D'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.D'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:4+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.d'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.d'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.d'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.d'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.d'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e \n', time(h), values(h,1), values(h,2), values(h,3));
                    end
                    fclose(file_id);
                    
                    k = k + 3;
                end
            end
            disp('Done');
            
        end
        
    end
    
    % VEOLCITY
    
    if(OUT_OPT(2) == 1)
        
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Velocity...']);
            filename1 = 'MONITORXXXXX.V';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.V'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.V'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.V'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.V'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.V'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:4+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.v'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.v'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.v'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.v'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.v'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e \n', time(h), values(h,1), values(h,2), values(h,3));
                    end
                    fclose(file_id);
                    
                    k = k + 3;
                end
            end
            disp('Done');
        end
        
        
    end
    
    % ACCELERATION
    
    if(OUT_OPT(3) == 1)
        
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Acceleration...']);
            filename1 = 'MONITORXXXXX.A';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.A'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.A'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.A'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.A'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.A'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:4+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.a'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.a'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.a'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.a'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.a'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e \n', time(h), values(h,1), values(h,2), values(h,3));
                    end
                    fclose(file_id);
                    
                    k = k + 3;
                end
            end
            disp('Done');
        end
        
    end
    
    
    % STRAIN
    
    if(OUT_OPT(4) == 1)
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Strain tensor...']);
            filename1 = 'MONITORXXXXX.E';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.E'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.E'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.E'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.E'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.E'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:7+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.e'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.e'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.e'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.e'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.e'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e %10.8e   %10.8e   %10.8e \n', ...
                            time(h), values(h,1), values(h,2), values(h,3), values(h,4), values(h,5), values(h,6));
                    end
                    fclose(file_id);
                    
                    k = k + 6;
                end
            end
            disp('Done');
        end
        
    end
    
    
    % STRESS
    
    if(OUT_OPT(5) == 1)
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Stress tensor...']);
            filename1 = 'MONITORXXXXX.S';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.S'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.S'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.S'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.S'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.S'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:7+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.s'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.s'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.s'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.s'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.s'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e %10.8e   %10.8e   %10.8e \n', ...
                            time(h), values(h,1), values(h,2), values(h,3), values(h,4), values(h,5), values(h,6));
                    end
                    fclose(file_id);
                    
                    k = k + 6;
                end
            end
            disp('Done');
        end
        
    end
    
    
    
    % OMEGA
    
    if(OUT_OPT(6) == 1)
        
        
        Num_of_tot_mon = 0;
        
        for i = 1 : MPI_num_proc
            disp(['Processing MONITOR ',num2str(i),' for Rotations...']);
            filename1 = 'MONITORXXXXX.O';
            filename2 = 'MONITORXXXXX.INFO';
            
            if(MPI_mnt_id(i) ~= 0 && i <= 10)
                filename1 = ['MONITOR0000',num2str(i-1),'.O'];
                filename2 = ['MONITOR0000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100)
                filename1 = ['MONITOR000',num2str(i-1),'.O'];
                filename2 = ['MONITOR000',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 1000)
                filename1 = ['MONITOR00',num2str(i-1),'.O'];
                filename2 = ['MONITOR00',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 10000)
                filename1 = ['MONITOR0',num2str(i-1),'.O'];
                filename2 = ['MONITOR0',num2str(i-1),'.INFO'];
            elseif(MPI_mnt_id(i) ~= 0 && i <= 100000)
                filename1 = ['MONITOR',num2str(i-1),'.O'];
                filename2 = ['MONITOR',num2str(i-1),'.INFO'];
            end
            
            fid = fopen([path_start,filename2]);
            
            if(fid ~= -1)
                
                ST = fclose(fid);
                INFO_MONITOR = load([path_start,filename2]);
                Num_of_mon = INFO_MONITOR(1);
                Id_of_mon = INFO_MONITOR([2:2+Num_of_mon-1]);
                Num_of_tot_mon = Num_of_tot_mon + Num_of_mon;
                
                
                VAL_MONITOR = load([path_start,filename1]);
                k = 0;
                for j = 1: Num_of_mon
                    time = VAL_MONITOR(:,1);
                    values = VAL_MONITOR(:,[2+k:4+k]);
                    
                    
                    if(Id_of_mon(j) < 10)
                        datafilename = ['monitor0000',num2str(Id_of_mon(j)),'.o'];
                    elseif(Id_of_mon(j) < 100)
                        datafilename = ['monitor000',num2str(Id_of_mon(j)),'.o'];
                    elseif(Id_of_mon(j) < 1000)
                        datafilename = ['monitor00',num2str(Id_of_mon(j)),'.o'];
                    elseif(Id_of_mon(j) < 10000)
                        datafilename = ['monitor0',num2str(Id_of_mon(j)),'.o'];
                    elseif(Id_of_mon(j) < 100000)
                        datafilename = ['monitor',num2str(Id_of_mon(j)),'.o'];
                    end
                    
                    file_id = fopen([path_end,datafilename], 'w');
                    for h = 1: length(time)
                        fprintf(file_id, '%10.8e   %10.8e   %10.8e  %10.8e \n', time(h), values(h,1), values(h,2), values(h,3));
                    end
                    fclose(file_id);
                    
                    k = k + 3;
                end
            end
            disp('Done');
        end
        
        
    end
    return
end
