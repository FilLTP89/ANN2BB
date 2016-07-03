%% *Parse recorded motions at KKNPP*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _rc_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _bhr (borehole data structure)_
%% OUTPUT:
% * _bhr (borehole data structure)_
%% N.B.
% Need for _itaca_monitor_name.m_, _parse_itaca_file.m_,
% _band_pass_filter.m_, _vel2acc.m_, _dis2acc.m_, _syn2ann_thp.m_

function [varargout] = rc_parser(varargin)
    %% *SET-UP*
    % _borehole data structure_
    bhr = varargin{1};
    bhr.na = bhr.ns*bhr.nd;
    %% *PARSING RECORDS*
    count_nm = 0;
    rec.mon.na = bhr.na;
    switch bhr.tp{1}
        case 'itaca'
            % stations
            for i_ = 1:bhr.ns
                % devices
                for j_ = 1:bhr.nd(i_)
                    count_nm = count_nm+1;
                    % _monitor identity_
                    rec.mon.nc=0;
                    rec.mon.nr=0;
                    % _directions_
                    for ii_ = 1:bhr.nc
                        rec.mon.cp{ii_}=bhr.cp{ii_};
                        rec.mon.nc = rec.mon.nc+1;
                        % _components_
                        for jj_ = 1:bhr.nr
                            rec.mon.rc{jj_} = bhr.rc{jj_};
                            rec.mon.nr = rec.mon.nr+1;
                            bhr.nm{count_nm} = itaca_monitor_name...
                                (bhr.st{i_}.id{1},bhr.ev{1},bhr.st{i_}.dv{j_},...
                                bhr.cp{ii_},bhr.rc{jj_},bhr.st{i_}.ni{1},bhr.pt);
                            fprintf('filename: %s\n',bhr.nm{count_nm});
                            % _parsing_
                            [rec.mon.dtm(count_nm),...
                                rec.mon.ntm(count_nm),...
                                rec.mon.vtm{count_nm},...
                                rec.syn{count_nm}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})] = ...
                                parse_itaca_file(bhr.nm{count_nm});
                        end
                    end
                end
            end
        case 'kknpp'
            % [TODO]
    end
    %% *ACCELERATION/VELOCITY/DISPLACEMENT*
    switch strcat(cell2mat(bhr.rc(:)'))
        case 'a'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        band_pass_filter(rec.mon.dtm(i_),...
                        rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
        case 'v'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        vel2acc(rec.mon.dtm(i_),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
        case 'd'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        dis2acc(rec.mon.dtm(i_),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
        case 'av'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        band_pass_filter(rec.mon.dtm(i_),...
                        rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
        case 'ad'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        band_pass_filter(rec.mon.dtm(i_),...
                        rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
        case 'vd'
            for i_ = 1:rec.mon.na
                for j_ = 1:rec.mon.nc
                    [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                        dis2acc(rec.mon.dtm(i_),...
                        rec.syn{i_}.thd.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
                end
            end
    end
    rec.bhr = bhr;
    %% *OUTPUT*
    varargout{1} = bhr;
    varargout{2} = rec;
    return
end
