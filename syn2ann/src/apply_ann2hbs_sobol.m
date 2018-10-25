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

function [out] = apply_ann2hbs_sobol(x)
    global no contr epsilon1 
    global ann hbs dsx s_
    for i_ = 1:hbs.clc
        out(:,i_) = 10.^(sim(ann.(hbs.mon.cp{dsx}).net,x(i_,:).'));
    end
    out = out(s_,:)./100;
    return
end
