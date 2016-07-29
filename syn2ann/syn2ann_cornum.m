%% *Process time-histories*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_cornum_: function to process numerical time-histories and
% correct them
%% INPUT:
% * _sas (syn2ann structure)_
%% OUTPUT:
% * _sas (syn2ann structure)_
%% N.B.
% Need for _PGAVD_eval.m_, _arias_intensity.m_
function [varargout] = syn2ann_cornum(varargin)
    %% *SET-UP*
    sas = varargin{1};
    for i_ = 1:sas.mon.na
        for j_ = 1:sas.mon.nc
            %% ARIAS INTENSITY
            [sas.syn{i_}.tha.(sas.mon.cp{j_}),....
                sas.syn{i_}.thv.(sas.mon.cp{j_}),...
                sas.syn{i_}.thd.(sas.mon.cp{j_})] = ...
                band_pass_filter(sas.mon.dtm(i_),sas.syn{i_}.tha.(sas.mon.cp{j_}),...
                [],[]);
        end
    end
    %% *OUTPUT*
    varargout{1} = sas;
    return
end
