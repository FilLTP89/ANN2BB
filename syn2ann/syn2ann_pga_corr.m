%% *Process time-histories*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_pga_corr_: function to process time-histories
%% INPUT:
% * _sas (syn2ann structure)_
%% OUTPUT:
% * _sas (syn2ann structure)_
%% N.B.
% Need for _PGAVD_eval.m_, _arias_intensity.m_
function [varargout] = syn2ann_pga_corr(varargin)
    %% *SET-UP*
    sas_inp = varargin{1};
    sas_tar = varargin{2};
    
    for i_ = 1:sas_inp.mon.na
        for j_ = 1:sas_inp.mon.nc
            cp = sas_inp.mon.cp{j_};
            [sas_inp.syn{i_}.tha.(cp)] = ...
                adjust_pga(sas_inp.mon.dtm(i_),...
                sas_inp.syn{i_}.tha.(cp),...
                sas_tar.syn{i_}.psa.(cp)(1),1);
            sas_inp.syn{i_}.psa.(cp)(1) = max(abs(sas_inp.syn{i_}.tha.(cp)));
        end
    end
    %% *OUTPUT*
    varargout{1} = sas_inp;
    return
end
