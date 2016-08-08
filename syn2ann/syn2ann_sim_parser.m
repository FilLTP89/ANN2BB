%% *Parse outcomes of numerical analysis*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_sim_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _mon (monitor structure)_
%% OUTPUT:
% * _mon (monitor structure)_
% * _nss (structure of numerical simulations)_
%% N.B.
% Need for _speed_monitor_name.m,hisada_monitor_name.m,vel2acc.m,
% dis2acc.m,PGAVD_eval.m_
function [varargout] = syn2ann_sim_parser(varargin)
    I1 = 0.05;
    %% SET-UP
    % _monitor structure_
    mon = varargin{1};
    bhr = varargin{2};
    %%
    % _parsing metadata_
    mtd = importdata(mon.fn);
    idm = find(strcmpi('Monitor ID',mtd.textdata(1,:))==1);
    idd = find(strcmpi('Repi (km)',mtd.textdata(1,:))==1);
    ids = strcmpi('Station code',mtd.textdata(1,:));
    idg.eutm = find(strcmpi('E_UTM [m]',mtd.textdata(1,:))==1);
    idg.nutm = find(strcmpi('N_UTM [m]',mtd.textdata(1,:))==1);
    
    for i_ = 1:mon.na
        idn = find(round(mtd.data(:,idm-1))==mon.id(i_)==1);
        mon.dep(i_) = mtd.data(idn,idd-1);
        mon.st(i_) = mtd.textdata(idn+1,ids);
        mon.eutm(i_) = mtd.data(idn,idg.eutm-1);
        mon.nutm(i_) = mtd.data(idn,idg.nutm-1);
        if strcmpi(mon.st(i_),bhr.nm(i_))
            fprintf('monitor %u matched to record %s!\n',mon.id(i_),bhr.nm{i_});
        else
            warning('monitor and record do not match');
            keyboard
        end
    end
    nss.mon = mon;
    %%
    % _parsing records_
    switch(lower(mon.tp))
        case 's'
            %%
            % _speed simulations_
            
            for i_ = 1:mon.na % number of monitors
                fprintf('monitor: \n');
                disp(mon.id(i_));
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
                % name of the stations
                fprintf('components: \n');
                for j_ = 1:mon.nc
                    cpp = mon.cp{j_};
                    disp(cpp);
                    nss.syn{i_}.(strcat('th',mon.rc{1})).(cpp) = str(:,mon.ci(j_)+1);
                    switch mon.rc{1}
                        case 'v'
                            [nss.syn{i_}.tha.(cpp),...
                                nss.syn{i_}.thv.(cpp),...
                                nss.syn{i_}.thd.(cpp)] = ...
                                vel2acc(nss.mon.dtm(i_),...
                                nss.syn{i_}.thv.(cpp),mon.lfr,mon.hfr);
                        case 'd'
                            [nss.syn{i_}.tha.(cpp),...
                                nss.syn{i_}.thv.(cpp),...
                                nss.syn{i_}.thd.(cpp)] = ...
                                dis2acc_speed(nss.mon.dtm(i_),...
                                nss.syn{i_}.thd.(cpp),mon.lfr,mon.hfr);
                    end
                end
            end
        case 'h'
            %%
            % _hisada simulations_
            for i_ = 1:mon.na % number of monitors
                flag = 1;
                for j_ = 1:mon.nc % number of components
                    cpp = mon.cp{j_};
                    % file name
                    str = hisada_monitor_name(mon.id(i_),cpp,mon.pt);
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
                    nss.syn{i_}.(strcat('th',mon.rc{1})).(cpp) = ...
                        str(:,2);
                    [nss.syn{i_}.tha.(cpp),...
                        nss.syn{i_}.thv.(cpp),...
                        nss.syn{i_}.thd.(cpp)] = ...
                        vel2acc(nss.mon.dtm(i_),...
                        nss.syn{i_}.thv.(cpp),mon.lfr,mon.hfr);
                end
            end
    end
    nss.mon.rc  = {'a';'v';'d'};
    nss.mon.nr  = numel(nss.mon.rc);
    %% OUTPUT
    varargout{1} = mon;
    varargout{2} = nss;
    return
end
