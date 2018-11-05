%% *Parse outcomes of numerical analysis*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_pbs_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _mon (monitor structure)_
%% OUTPUT:
% * _mon (monitor structure)_
% * _pbs (structure of numerical simulations)_
%% N.B.
% Need for _speed_monitor_name.m,hisada_monitor_name.m,vel2acc.m,
% dis2acc.m,PGAVD_eval.m_
function [varargout] = syn2ann_pbs_parser(varargin)
    global flag_map
    
    I1 = 0.05;
    %% SET-UP
    % _monitor structure_
    mon = varargin{1};
    bhr = varargin{2};
    %%
    % _parsing metadata_
    
    mtd = importdata(mon.fnm);
    idm = find(strcmpi('Monitor ID',mtd.textdata(1,:))==1);
    idd = find(strcmpi('Repi [km]',mtd.textdata(1,:))==1);
    ids = find(strcmpi('Station code',mtd.textdata(1,:))==1);
    idg.eutm = find(strcmpi('E_UTM [m]',mtd.textdata(1,:))==1);
    idg.nutm = find(strcmpi('N_UTM [m]',mtd.textdata(1,:))==1);
    
    for i_ = 1:mon.na
        idn = find(round(mtd.data(:,idm-1))==mon.id(i_)==1,1,'first');
        mon.dep(i_) = mtd.data(idn,idd-1);
        mon.st(i_) = mtd.textdata(idn+1,ids);
        mon.eutm(i_) = mtd.data(idn,idg.eutm-1);
        mon.nutm(i_) = mtd.data(idn,idg.nutm-1);
        
        if ~flag_map
            if strcmpi(mon.st(i_),bhr.nm(i_))
                fprintf('monitor %u matched to record %s!\n',mon.id(i_),bhr.nm{i_});
            else
                warning('monitor and record do not match');
                keyboard
            end
        end
    end
    pbs.mon = mon;
    %%
    % _parsing records_
    switch(lower(mon.typ))
        case 'speed'
            %%
            % _speed simulations_
            
            for i_ = 1:mon.na % number of monitors
                fprintf('monitor: \n');
                disp(mon.id(i_));
                % file name
                str = speed_monitor_name(mon.id(i_),mon.rc{1},mon.pt);
                % read monitor file
                str = importdata(str);
                % time-step
                pbs.mon.dtm(i_) = mean(diff(str(:,1)));
                % parse acceleration components
                fprintf('components SPEED: \n');
                for j_ = 1:mon.nc
                    cpp = mon.cp{j_};
                    disp(cpp);
                    pbs.syn{i_}.(strcat('th',mon.rc{1})).(cpp) = str(:,mon.ci(j_)+1);
                    switch mon.rc{1}
                        case 'v'
                            % [TODO]
                        case 'd'
                            [pbs.syn{i_}.tha.(cpp),...
                                pbs.syn{i_}.thv.(cpp),...
                                pbs.syn{i_}.thd.(cpp)] = ...
                                bpf_thd_speed(pbs.mon.dtm(i_),...
                                pbs.syn{i_}.thd.(cpp),mon.lfr,mon.hfr);
                    end
                end
                % time-step number
                pbs.mon.ntm(i_) = numel(pbs.syn{i_}.tha.(cpp));
                % time-vector
                pbs.mon.vtm{i_} = pbs.mon.dtm(i_)*(0:pbs.mon.ntm(i_)-1);
                
            end
        case 'sem3d'
            %%
            % _sem3d simulations_
            
            %% *PARSING*
            sem3d = parse_sem_results('pfn', mon.pt);
            
            %% *RESAMPLING*
            sem3d = sem_rsmpl(sem3d,'dtt',0.005);
            
            %% *BASELINE CORRECTION + FILTERING*
            sem3d = sem_bpf(sem3d,'hfr',[],'lfr',0.05,'avd','idc');
            rc_sem3d = {'Accel';'Veloc';'Displ'};
            cpn_sem3d = {'ew','x';'ns','y';'ud','z'};
            
            for i_ = 1:mon.na % number of monitors
                fprintf('monitor: \n');
                disp(mon.id(i_));
                
                % time-step
                pbs.mon.dtm(i_) = sem3d.dTime;
                
                % parse acceleration components
                fprintf('components SEM3D: \n');
                for k_=1:numel(rc_sem3d)
                    for j_ = 1:mon.nc
                        cpp = mon.cp{j_};
                        idx_sem3d = find(strcmpi(cpn_sem3d(:,1),cpp)==1);
                        disp(cpp);
                        pbs.syn{i_}.(strcat('th',lower(rc_sem3d{k_}(1)))).(cpp) = ...
                            sem3d.(rc_sem3d{k_}).(cpn_sem3d{idx_sem3d,2})(:,mon.id(i_));
                        
                    end
                end
                % time-step number
                pbs.mon.ntm(i_) = numel(pbs.syn{i_}.tha.(cpp));
                % time-vector
                pbs.mon.vtm{i_} = pbs.mon.dtm(i_)*(0:pbs.mon.ntm(i_)-1);
            end
    end
    pbs.mon.rc  = {'a';'v';'d'};
    pbs.mon.nr  = numel(pbs.mon.rc);
    %% OUTPUT
    varargout{1} = mon;
    varargout{2} = pbs;
    return
end
