%% *Train ANN network onto synthetics*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _apply_ann2hbs_justPSA_: function to apply trained ANN onto hybrid LF/HF
% synthetics, upon PSA(LF) values
%% INPUT:
% * _hbs (hybrid synthetics structure)_
% * _ann (trained Artificial Neural Network (ANN) structure)_
%% OUTPUT:
% * _trs (trained/simulated ann structure)_
%% *N.B.*
% Need for:
% _sim.m,trann_check_vTn.m_

function [varargout] = apply_ann2hbs_justPSA(varargin)
    %% *SET-UP*
    hbs = varargin{1};
    ann = varargin{2};
    %
    % _check input/target natural periods with hbs_
    %
    inn = cell(hbs.mon.nc,1);
    for j_ = 1:hbs.mon.nc
        %
        cpp.hbs = (hbs.mon.cp{j_});
        %
        [trs.(cpp.hbs).iid,trs.(cpp.hbs).tid] = ...
            trann_check_vTn(ann.(cpp.hbs).inp,ann.(cpp.hbs).tar,hbs.mon,1e-8);
        %
        psa = -999*ones(hbs.mon.na,hbs.mon.nT);
        for i_ = 1:hbs.mon.na
            psa(i_,:) = hbs.syn{i_}.psa.(cpp.hbs)(:)';
        end
        %
        inn{j_} = psa(:,trs.(cpp.hbs).iid);
    end
    
    %% *ANN SIMULATION*
    trs.mon.cp = hbs.mon.cp;
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_}.*100);
        trs.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp;
        trs.(hbs.mon.cp{j_}).mon.vTn = ...
            [ann.(hbs.mon.cp{j_}).tar.vTn;ann.(hbs.mon.cp{j_}).inp.vTn];
        trs.(hbs.mon.cp{j_}).mon.nT = numel(trs.(hbs.mon.cp{j_}).mon.vTn);
        trs.(hbs.mon.cp{j_}).mon.tol = ann.(hbs.mon.cp{j_}).tol;
        trs.(hbs.mon.cp{j_}).mon.nit = ann.(hbs.mon.cp{j_}).nit;
        for i_ = 1:hbs.mon.na
            out = 10.^(sim(ann.(hbs.mon.cp{j_}).net,inp(i_,:)'));
            out = out./100;
            trs.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}) = ...
                [out(:);inn{j_}(i_,:)'];
        end
    end
    
    %% *OUTPUT*
    varargout{1} = trs;
    return
end
