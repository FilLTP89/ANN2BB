%% *Parse recorded motions at KKNPP*
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

function [varargout] = syn2ann_bhr_parser(varargin)
    %% *SET-UP*
    % _borehole data structure_
    bhr = varargin{1};
    %% *PARSING RECORDS*
    count_na = 0;
    count_fn = 0;
    switch bhr.tp{1}
        case 'itaca'
            % stations
            for i_ = 1:bhr.ns
                % devices
                for j_ = 1:bhr.nd(i_)
                    count_na = count_na+1;
                    % _monitor identity_
                    bhr.nm{count_na} = sprintf('%s%s',bhr.st{i_}.id{1},bhr.st{i_}.dv{j_});
                    % _directions_
                    for ii_ = 1:bhr.nc
                        % _components_
                        for jj_ = 1:bhr.nr
                            count_fn = count_fn+1;
                            % _file-name_
                            bhr.fn{count_fn} = itaca_monitor_name...
                                (bhr.st{i_}.id{1},bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},...
                                bhr.cp{ii_},bhr.rc{jj_},bhr.st{i_}.ni,bhr.pt);
                            fprintf('filename: %s\n',bhr.fn{count_fn});
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
                    bhr.nm{count_na} = sprintf('%s%s',bhr.st{i_}.id{1},bhr.st{i_}.dv{j_});
                    % _components_
                    for jj_ = 1:bhr.nr
                        count_fn = count_fn+1;
                        % _file-name_
                        bhr.fn{count_fn} = kknpp_monitor_name(bhr.st{i_}.id{1},...
                            bhr.st{i_}.ev{1},bhr.st{i_}.dv{j_},bhr.pt);
                        fprintf('filename: %s\n',bhr.fn{count_fn});
                    end
                end
            end
    end
    %% *OUTPUT*
    varargout{1} = bhr;
    return
end
