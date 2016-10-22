%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_rec_parser_: function to parse synthetics from speed/hisada
%% INPUT:
% * _bhr (borehole data structure)_
%% OUTPUT:
% * _bhr (borehole data structure)_
%% N.B.
% Need for _itaca_monitor_name.m_, _parse_itaca_file.m_,
% _band_pass_filter.m_, _vel2acc.m_, _dis2acc.m_, _syn2ann_thp.m_

function [varargout] = syn2ann_rec_parser(varargin)
    %% *SET-UP*
    % _borehole data structure_
    bhr = varargin{1};
    cor = varargin{2};
    %% *PARSING RECORDS*
    count_na = 0;
    count_fn = 0;
    rec.mon.na = bhr.na;
    switch bhr.tp{1}
        case 'itaca'
            % stations
            for i_ = 1:bhr.ns
                % devices
                for j_ = 1:bhr.nd(i_)
                    count_na = count_na+1;
                    % _monitor identity_
                    rec.mon.nc=0;
                    rec.mon.nr=0;
                    bhr.nm{count_na} = sprintf('%s%s',bhr.st{i_}.id{1},bhr.st{i_}.dv{j_});
                    % _directions_
                    for ii_ = 1:bhr.nc
                        rec.mon.cp{ii_} = bhr.cp{ii_};
                        rec.mon.nc = rec.mon.nc+1;
                        % _components_
                        for jj_ = 1:bhr.nr
                            count_fn = count_fn+1;
                            rec.mon.rc{jj_} = bhr.rc{jj_};
                            rec.mon.nr = rec.mon.nr+1;
                            % _file-name_
                            bhr.fn{count_fn} = itaca_monitor_name...
                                (bhr.st{i_}.id{1},bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},...
                                bhr.cp{ii_},bhr.rc{jj_},bhr.st{i_}.ni,bhr.pt);
                            fprintf('filename: %s\n',bhr.fn{count_fn});
                            switch strcat(bhr.st{i_}.id{1},bhr.st{i_}.dv{j_})
                                case {'MRN','MIR08'}
                                rec.mon.dtm(count_na)=0.005;
                                rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})=...
                                    cor.(bhr.nm{count_na}).(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_});
                                rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})= ...
                                    rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})...
                                    (~isnan(rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})));
                                rec.mon.ntm(count_na)=numel(rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_}));
                                rec.mon.vtm{count_na} = (0:rec.mon.ntm(count_na)-1)*rec.mon.dtm(count_na);
                            case {'AQK','AQU'}
                                % _parsing_
                                [rec.mon.dtm(count_na),...
                                    rec.mon.ntm(count_na),...
                                    rec.mon.vtm{count_na},...
                                    rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})] = ...
                                    parse_itaca_file_new(bhr.fn{count_fn});
                            end
                        end
                    end
                end
            end
        case 'kknpp'
            % stations
            for i_ = 1:bhr.ns
                % devices
                for j_ = 1:bhr.nd(i_)
                    count_na = count_na+1;
                    % _monitor identity_
                    rec.mon.nc=0;
                    rec.mon.nr=0;
                    bhr.nm{count_na} = sprintf('%s%s',bhr.st{i_}.id{1},bhr.st{i_}.dv{j_});
                    % _components_
                    for jj_ = 1:bhr.nr
                        count_fn = count_fn+1;
                        rec.mon.rc{jj_} = bhr.rc{jj_};
                        rec.mon.nr = rec.mon.nr+1;
                        % _file-name_
                        bhr.fn{count_fn} = kknpp_monitor_name(bhr.st{i_}.id{1},...
                            bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},bhr.pt);
                        fprintf('filename: %s\n',bhr.fn{count_fn});
                        % _parsing_
                        rec.mon.cp = bhr.cp;
                        rec.mon.nc = bhr.nc;
                        [rec.mon.dtm(count_na),...
                            rec.mon.ntm(count_na),...
                            rec.mon.vtm{count_na},...
                            rec.syn{count_na}.(strcat('th',bhr.rc{jj_}))] = ...
                            parse_kknpp_file(bhr.fn{count_fn},bhr.cp);
                    end
                end
            end
    end
    
    %% *ACCELERATION/VELOCITY/DISPLACEMENT*
    if strcmpi(strcat(cell2mat(bhr.rc(:)')),'a')
        for i_ = 1:rec.mon.na
            for j_ = 1:rec.mon.nc
                keyboard
                [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                    band_pass_filter(rec.mon.dtm(i_),...
                    rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
            end
        end
    elseif strcmpi(strcat(cell2mat(bhr.rc(:)')),'v')
        for i_ = 1:rec.mon.na
            for j_ = 1:rec.mon.nc
                [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                    vel2acc(rec.mon.dtm(i_),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
            end
        end
    elseif strcmpi(strcat(cell2mat(bhr.rc(:)')),'d')
        for i_ = 1:rec.mon.na
            for j_ = 1:rec.mon.nc
                [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                    dis2acc(rec.mon.dtm(i_),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
            end
        end
    elseif strcmpi(strcat(cell2mat(bhr.rc(:)')),'av')
        for i_ = 1:rec.mon.na
            for j_ = 1:rec.mon.nc
                [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                    band_pass_filter(rec.mon.dtm(i_),...
                    rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
            end
        end
    elseif strcmpi(strcat(cell2mat(bhr.rc(:)')),'ad')
        for i_ = 1:rec.mon.na
            for j_ = 1:rec.mon.nc
                [rec.syn{i_}.tha.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thv.(rec.mon.cp{j_}),...
                    rec.syn{i_}.thd.(rec.mon.cp{j_})] = ...
                    band_pass_filter(rec.mon.dtm(i_),...
                    rec.syn{i_}.tha.(rec.mon.cp{j_}),bhr.lfr,bhr.hfr);
            end
        end
    elseif strcmpi(strcat(cell2mat(bhr.rc(:)')),'vd')
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
