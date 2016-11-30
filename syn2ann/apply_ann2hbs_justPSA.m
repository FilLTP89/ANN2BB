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
    for j_ = 1:hbs.mon.nc
        cpp = (hbs.mon.cp{j_});
        [trs.(cpp).iid,trs.(cpp).tid] = ...
            trann_check_vTn(ann.(cpp).inp,ann.(cpp).tar,hbs.mon,1e-8);
    end
    
    inn = cell(hbs.mon.nc,1);
    
    for j_ = 1:hbs.mon.nc
        psa = -999*ones(hbs.mon.na,hbs.mon.nT);
        cpp = hbs.mon.cp{j_};
        for i_ = 1:hbs.mon.na
            psa(i_,:) = hbs.syn{i_}.psa.(cpp)(:)';
        end
        inn{j_} = psa(:,trs.(cpp).iid);
    end
    
    %% *ANN SIMULATION*
    trs.mon.cp = hbs.mon.cp;
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_}.*100);
        trs.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp;
        trs.(hbs.mon.cp{j_}).mon.vTn = ...
            [ann.(hbs.mon.cp{j_}).tar.vTn;ann.(hbs.mon.cp{j_}).inp.vTn];
        trs.(hbs.mon.cp{j_}).mon.nT = numel(trs.(hbs.mon.cp{j_}).mon.vTn);
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
