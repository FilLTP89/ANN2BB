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

function [varargout] = trann_rec_parser(varargin)
    
    %% *SET-UP*
    % _borehole data structure_
    bhr = varargin{1};
    %% *PARSING RECORDS*
    count_na = 0;
    count_fn = 0;
    rec.mon.na = bhr.na;
    % stations
    
    for i_ = 1:bhr.ns
        switch lower(bhr.st{i_}.tp{1})
            case 'avd'
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
                            bhr.fn{count_fn} = fullfile(bhr.pt,strcat(bhr.st{i_}.id{1},...
                                bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},'_',...
                                bhr.cp{ii_},bhr.st{i_}.ni{2}));
                            fprintf('filename: %s\n',bhr.fn{count_fn});
                            
                            % _parsing_
                            [rec.mon.dtm(count_na),...
                                rec.mon.ntm(count_na),...
                                rec.mon.vtm{count_na},...
                                rec.syn{count_na}.tha.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thv.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thd.(bhr.cp{ii_})] = ...
                                parse_avd_file(bhr.fn{count_fn},0.01);
                            
                        end
                    end
                end
            case {'itaca'}
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
                        cpn = bhr.cp{ii_}(1);
                        if strcmpi(cpn,'u')
                            cpn = 'z';
                        end
                        % _components_
                        for jj_ = 1:bhr.nr
                            count_fn = count_fn+1;
                            rec.mon.rc{jj_} = bhr.rc{jj_};
                            rec.mon.nr = rec.mon.nr+1;
                            % _file-name_
                            bhr.fnm{count_fn} = itaca_monitor_name...
                                (bhr.st{i_}.id{1},bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},...
                                cpn,bhr.rc{jj_},bhr.st{i_}.ni,bhr.pt);
                            fprintf('filename: %s\n',bhr.fnm{count_fn});
                            % _parsing_
                            [rec.mon.dtm(count_na),...
                                rec.mon.ntm(count_na),...
                                rec.mon.vtm{count_na},...
                                rec.syn{count_na}.(strcat('th',bhr.rc{jj_})).(bhr.cp{ii_})] = ...
                                parse_itaca_file_new(bhr.fnm{count_fn});
                            rec.syn{count_na}.tha.(bhr.cp{ii_}) = ...
                                rec.syn{count_na}.tha.(bhr.cp{ii_});
                            [rec.syn{count_na}.tha.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thv.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thd.(bhr.cp{ii_})] = ...
                                idc_tha(rec.mon.dtm(count_na),...
                                rec.syn{count_na}.tha.(bhr.cp{ii_}));
                        end
                    end
                end
            case {'knet','kiknet'}
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
                            bhr.fn{count_fn} = fullfile(bhr.pt,strcat(bhr.st{i_}.id{1},...
                                bhr.st{i_}.ni{1},bhr.st{i_}.ev{1},...
                                '.',upper(bhr.cp{ii_}),bhr.st{i_}.dv{j_}));
                            fprintf('filename: %s\n',bhr.fn{count_fn});
                            
                            % _parsing_
                            [rec.mon.dtm(count_na),...
                                rec.mon.ntm(count_na),...
                                rec.mon.vtm{count_na},...
                                rec.syn{count_na}.tha.(bhr.cp{ii_}),~] = ...
                                parse_kiknet_file(bhr.fn{count_fn});
                            rec.syn{count_na}.tha.(bhr.cp{ii_}) = ...
                                rec.syn{count_na}.tha.(bhr.cp{ii_})/100;
                            [rec.syn{count_na}.tha.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thv.(bhr.cp{ii_}),...
                                rec.syn{count_na}.thd.(bhr.cp{ii_}),...
                                rec.mon.vtm{count_na},rec.mon.pads{count_na}] = ...
                                bpf_tha(rec.mon.dtm(count_na),...
                                rec.syn{count_na}.tha.(bhr.cp{ii_}),0.05,[]);
                            
                            rec.mon.vtm{count_na} = rec.mon.vtm{count_na}-...
                                rec.mon.vtm{count_na}(rec.mon.pads{count_na}(1)+1);
                        end
                    end
                end
        end
        bhr.st{i_}.ev{1}(strfind(bhr.st{i_}.ev{1},'.')) = '_';
    end
    rec.bhr = bhr;
    %% *OUTPUT*
    varargout{1} = bhr;
    varargout{2} = rec;
    return
end
