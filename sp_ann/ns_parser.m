%% *Parse outcomes of numerical analysis*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ns_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _mon    = monitor structure_
% * mon.pt = path to monitor files    (string)
% * mon.fn = monitor metadata filename(string)
% * mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * mon.id = monitor identity         (integer)
% * mon.na = number of monitors       (integer)
% * mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * mon.nr = number of records        (integer)
% * mon.cp = motion component         (integer: 1,2,3)
% * mon.nc = number of components     (integer)
%% OUTPUT:
% * _mon    = monitor structure_
% * mon.pt = path to monitor files    (string)
% * mon.fn = monitor metadata filename(string)
% * mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * mon.id = monitor identity         (integer)
% * mon.dep = epicentral distance     (real vector)
% * mon.stn = monitor names           (string vector)
% * mon.na = number of monitors       (integer)
% * mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * mon.nr = number of records        (integer)
% * mon.cp = motion component         (integer: 1,2,3)
% * mon.nc = number of components     (integer)
%%
% * _nss    = structure of numerical simulations_
%%
% * _nss.mon    = structure of monitor data_
% * nss.mon.pt = path to monitor files        (string)
% * nss.mon.fn = monitor metadata filename    (string)
% * nss.mon.tp = type of monitor              (string: 'S'(speed),'H'(hisada))
% * nss.mon.id = monitor identity             (integer vector)
% * nss.mon.dep = epicentral distance         (real vector)
% * nss.mon.stn = monitor names               (string vector)
% * nss.mon.na = number of monitors           (integer)
% * nss.mon.rc = monitor record               (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * nss.mon.nr = number of records            (integer)
% * nss.mon.cp = motion component             (integer: 1,2,3)
% * nss.mon.nc = number of components         (integer)
% * nss.mon.dtm(mon.na,1) = time-steps        (real vector)
% * nss.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
% * nss.mon.ntm(mon.na,1) = time step number  (real vector)
%%
% * _nss.syn(mon.na,1)     = structure vectors (cell vector)_
% * nss.syn{i}.tha.x      = x-acceleration    (real vector)
% * nss.syn{i}.tha.y      = y-acceleration    (real vector)
% * nss.syn{i}.tha.z      = z-acceleration    (real vector)
% * nss.syn{i}.thv.x      = x-velocity        (real vector)
% * nss.syn{i}.thv.y      = y-velocity        (real vector)
% * nss.syn{i}.thv.z      = z-velocity        (real vector)
% * nss.syn{i}.thd.x      = x-displacement    (real vector)
% * nss.syn{i}.thd.y      = y-displacement    (real vector)
% * nss.syn{i}.thd.z      = z-displacement    (real vector)
% * nss.syn{i}.pga.x(1)   = x-time-pga        (real vector)
% * nss.syn{i}.pga.x(2)   = x-pga             (real vector)
% * nss.syn{i}.pga.y(1)   = y-time-pga        (real vector)
% * nss.syn{i}.pga.y(2)   = y-pga             (real vector)
% * nss.syn{i}.pga.z(1)   = z-time-pga        (real vector)
% * nss.syn{i}.pga.z(2)   = z-pga             (real vector)
% * nss.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
% * nss.syn{i}.pgv.x(2)   = x-pgv             (real vector)
% * nss.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
% * nss.syn{i}.pgv.y(2)   = y-pgv             (real vector)
% * nss.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
% * nss.syn{i}.pgv.z(2)   = z-pgv             (real vector)
% * nss.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
% * nss.syn{i}.pgd.x(2)   = x-pgd             (real vector)
% * nss.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
% * nss.syn{i}.pgd.y(2)   = y-pgd             (real vector)
% * nss.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
% * nss.syn{i}.pgd.z(2)   = z-pgd             (real vector)
%% N.B.
% Need for _speed_monitor_name.m,hisada_monitor_name.m,vel2acc.m,
% dis2acc.m,PGAVD_eval.m_
function [varargout] = ns_parser(varargin)
    
    %% SET-UP
    %%
    % _monitor structure_
    mon = varargin{1};
    %%
    % _parsing metadata_
    mtd = importdata(mon.fn);
    idm = find(strcmpi('Monitor ID',mtd.textdata(1,:))==1)-2;
    idd = find(strcmpi('Repi (km)',mtd.textdata(1,:))==1)-2;
    ids = find(strcmpi('Station code',mtd.textdata(1,:))==1);
    for i_ = 1:mon.na
        idn = round(mtd.data(:,idm))==mon.id(i_);
        mon.dep(i_) = mtd.data(idn,idd);
        mon.stn(i_) = mtd.textdata(idn,ids);
    end
    nss.mon = mon;
    %%
    % _parsing records_
    switch(lower(mon.tp))
        case 's'
            %%
            % _speed simulations_
            for i_ = 1:mon.na % number of monitors
                % file name
                str = speed_monitor_name(mon.id(i_),mon.rc{1},mon.pt);
                % read monitor file
                str = importdata(str);
                % time-vector
                nss.mon.vtm(i_) = {str(:,1)};
                % time-step
                nss.mon.dtm(i_) = diff(nss.mon.vtm{i_}(1:2));
                % time-step number
                nss.mon.ntm(i_) = numel(nss.mon.vtm{i_});
                for k_ = 1:mon.nc % number of components
                    eval(sprintf('nss.syn{i_}.th%s.%s = str(:,%u);',...
                        mon.rc{1},mon.cp{k_},k_+1));
                    if strcmpi(mon.rc{1},'v')
                        eval(sprintf(['nss.syn{i_}.tha.%s,nss.syn{i_}.thv.%s,',...
                            'nss.syn{i_}.thd.%s] = vel2acc(nss.mon.dtm(i_),',...
                            'nss.syn{i_}.thv.%s,.1,3);'],...
                            mon.cp{k_},mon.cp{k_},mon.cp{k_},mon.cp{k_}));
                    elseif strcmpi(mon.rc{1},'d')
                        eval(sprintf(['[nss.syn{i_}.tha.%s,nss.syn{i_}.thv.%s,',...
                            'nss.syn{i_}.thd.%s] = dis2acc(nss.mon.dtm(i_),',...
                            'nss.syn{i_}.thd.%s,.1,3);'],...
                            mon.cp{k_},mon.cp{k_},mon.cp{k_},mon.cp{k_}));
                    end
                end
                
            end
        case 'h'
            %%
            % _hisada simulations_
            for i_ = 1:mon.na % number of monitors
                flag = 1;
                for k_ = 1:mon.nc % number of components
                    % file name
                    str = hisada_monitor_name(mon.id(i_),mon.cp{k_},mon.pt);
                    % read monitor file
                    str = importdata(str);
                    if flag    % time-vector
                        nss.mon.vtm(i_) = {str(:,1)};
                        % time-step
                        nss.mon.dtm(i_) = diff(nss.mon.vtm{i_}(1:2));
                        % time-step number
                        nss.mon.ntm(i_) = numel(nss.mon.vtm{i_});
                        flag = 0;
                    end
                    eval(sprintf('nss.syn{i_}.thv.%s = str(:,2);',...
                        mon.cp{k_}));
                    eval(sprintf(['[nss.syn{i_}.tha.%s,nss.syn{i_}.thv.%s,',...
                        'nss.syn{i_}.thd.%s] = vel2acc(nss.mon.dtm(i_),',...
                        'nss.syn{i_}.thv.%s,.1,2);'],...
                        mon.cp{k_},mon.cp{k_},mon.cp{k_},mon.cp{k_}));
                end
                
            end
    end
    nss.mon.rc  = {'a';'v';'d'};
    nss.mon.nr  = numel(nss.mon.rc);
    %% COMPUTING MAXIMUM VALUES
    %%
    % _compute PGA/PGV/PGD from numerical simulations
    for i_ = 1:mon.na % number of monitors
        for j_ = 1:mon.nc % number of components
            eval(sprintf(['[nss.syn{i_}.pga.%s(1),nss.syn{i_}.pga.%s(2),'...
                'nss.syn{i_}.pgv.%s(1),nss.syn{i_}.pgv.%s(2),',...
                'nss.syn{i_}.pgd.%s(1),nss.syn{i_}.pgd.%s(2)] = ',...
                'PGAVD_eval(nss.mon.dtm(i_),nss.syn{i_}.tha.%s,',...
                'nss.syn{i_}.thv.%s,nss.syn{i_}.thd.%s);'],...
                mon.cp{j_},mon.cp{j_},mon.cp{j_},mon.cp{j_},mon.cp{j_},...
                mon.cp{j_},mon.cp{j_},mon.cp{j_},mon.cp{j_}));
        end
    end
    %% OUTPUT
    varargout{1} = mon;
    varargout{2} = nss;
    return
end
