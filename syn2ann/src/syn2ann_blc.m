%% *Process time-histories*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_thp_: function to process time-histories
%% INPUT:
% * _sas (syn2ann structure)_
%% OUTPUT:
% * _sas (syn2ann structure)_
%% N.B.
% Need for _PGAVD_eval.m_, _arias_intensity.m_
function [varargout] = syn2ann_blc(varargin)
    %% *SET-UP*
    sas = varargin{1};
    for i_ = 1:sas.mon.na
        for j_ = 1:sas.mon.nc
            %% *BASELINE CORRECTION*
            [sas.syn{i_}.tha.(sas.mon.cp{j_}),....
                sas.syn{i_}.thv.(sas.mon.cp{j_}),...
                sas.syn{i_}.thd.(sas.mon.cp{j_}),...
                sas.mon.vtm{i_},...
                sas.mon.npd(:,i_)] = ...
                bpf_tha(sas.mon.dtm(i_),...
                sas.syn{i_}.tha.(sas.mon.cp{j_}),0.05,40);
            sas.mon.ntm(i_) = numel(sas.syn{i_}.tha.(sas.mon.cp{j_}));
            sas.mon.vtm{i_} = sas.mon.vtm{i_}-sas.mon.vtm{i_}(sas.mon.npd(1,i_)+1);
        end
    end
    %% *OUTPUT*
    varargout{1} = sas;
    return
end
