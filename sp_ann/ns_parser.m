%% *Parse outcomes of numerical analysis*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ns_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _mon (monitor structure)_
%% OUTPUT:
% * _mon (monitor structure)_
% * _nss (structure of numerical simulations)_
%% N.B.
% Need for _speed_monitor_name.m,hisada_monitor_name.m,vel2acc.m,
% dis2acc.m,PGAVD_eval.m_
function [varargout] = ns_parser(varargin)
    I1 = 0.05;
    %% SET-UP
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
        idn = idn(2:end-1);
        mon.dep(i_) = mtd.data(idn,idd);
        mon.st(i_) = mtd.textdata(idn,ids);
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
                
                for j_ = 1:mon.nc
                    nss.syn{i_}.(strcat('th',mon.rc{1})).(mon.cp{j_}) = str(:,mon.ci(j_)+1);
                    switch mon.rc{1}
                        case 'v'
                        [nss.syn{i_}.tha.(mon.cp{j_}),...
                            nss.syn{i_}.thv.(mon.cp{j_}),...
                            nss.syn{i_}.thd.(mon.cp{j_})] = ...
                            vel2acc(nss.mon.dtm(i_),...
                            nss.syn{i_}.thv.(mon.cp{j_}),mon.lfr,mon.hfr);
                    case 'd'
                       [nss.syn{i_}.tha.(mon.cp{j_}),...
                            nss.syn{i_}.thv.(mon.cp{j_}),...
                            nss.syn{i_}.thd.(mon.cp{j_})] = ...
                            dis2acc(nss.mon.dtm(i_),...
                            nss.syn{i_}.thd.(mon.cp{j_}),mon.lfr,mon.hfr);
                    end
%                     [nss.syn{i_}.AT5.(nss.mon.cp{j_}),....
%                         nss.syn{i_}.AI5.(nss.mon.cp{j_}),...
%                         nss.syn{i_}.Ain.(nss.mon.cp{j_})] = ...
%                         arias_intensity(nss.syn{i_}.tha.(nss.mon.cp{j_}),...
%                         nss.mon.dtm(i_),I1);
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
                    nss.syn{i_}.(strcat('th',mon.rc{1})).(mon.cp{k_}) = ...
                        str(:,2);
                    [nss.syn{i_}.tha.(mon.cp{k_}),nss.syn{i_}.thv.(mon.cp{k_}),...
                        nss.syn{i_}.thd.(mon.cp{k_})] = ...
                        vel2acc(nss.mon.dtm(i_),...
                        nss.syn{i_}.thv.(mon.cp{k_}),mon.lfr,mon.hfr);
%                     [nss.syn{i_}.AT5.(nss.mon.cp{j_}),....
%                         nss.syn{i_}.AI5.(nss.mon.cp{j_}),...
%                         nss.syn{i_}.Ain.(nss.mon.cp{j_})] = ...
%                         arias_intensity(nss.syn{i_}.tha.(nss.mon.cp{j_}),...
%                         nss.mon.dtm(i_),I1);
                end
            end
    end
    nss.mon.rc  = {'a';'v';'d'};
    nss.mon.nr  = numel(nss.mon.rc);
%     %%
%     % _compute PGA/PGV/PGD from numerical simulations
%     for i_ = 1:mon.na % number of monitors
%         for j_ = 1:mon.nc % number of components
%             [nss.syn{i_}.pga.(mon.cp{j_})(1),nss.syn{i_}.pga.(mon.cp{j_})(2),...
%                 nss.syn{i_}.pgv.(mon.cp{j_})(1),nss.syn{i_}.pgv.(mon.cp{j_})(2),...
%                 nss.syn{i_}.pgd.(mon.cp{j_})(1),nss.syn{i_}.pgd.(mon.cp{j_})(2)] = ...
%                 PGAVD_eval(nss.mon.dtm(i_),nss.syn{i_}.tha.(mon.cp{j_}),...
%                 nss.syn{i_}.thv.(mon.cp{j_}),nss.syn{i_}.thd.(mon.cp{j_}));
%         end
%     end
    %% OUTPUT
    varargout{1} = mon;
    varargout{2} = nss;
    return
end
