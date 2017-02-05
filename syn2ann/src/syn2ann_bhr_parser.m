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
    for i_ = 1:bhr.ns
        % devices
        for j_ = 1:bhr.nd(i_)
            count_na = count_na+1;
            % _monitor identity_
            bhr.nm{count_na} = '';
            % _directions_
            for ii_ = 1:bhr.nc
                % _components_
                for jj_ = 1:bhr.nr
                    count_fn = count_fn+1;
                    % _file-name_
                    
                end
            end
        end
    end
    %% *OUTPUT*
    varargout{1} = bhr;
    return
end
