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

function [varargout] = prepare_ann2hbs_sobol(varargin)
    %% *SET-UP*
    global srt
    hbs = varargin{1};
    ann = varargin{2};
    pdf = varargin{3};
    cv  = varargin{4};
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

    VVarEntree = cell(hbs.mon.na,hbs.mon.nc);

    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_}.*100);
        trs.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp;
        trs.(hbs.mon.cp{j_}).mon.vTn = ...
            [ann.(hbs.mon.cp{j_}).tar.vTn;ann.(hbs.mon.cp{j_}).inp.vTn];
        trs.(hbs.mon.cp{j_}).mon.nT = numel(trs.(hbs.mon.cp{j_}).mon.vTn);
        trs.(hbs.mon.cp{j_}).mon.tol = ann.(hbs.mon.cp{j_}).tol;
        trs.(hbs.mon.cp{j_}).mon.nit = ann.(hbs.mon.cp{j_}).nit;
        trs.(hbs.mon.cp{j_}).mon.TnC = ann.(hbs.mon.cp{j_}).TnC;
        %
        for i_ = 1:hbs.mon.na
            param = inp(i_,:)';
            VVarEntree{i_,j_} = ones(numel(param),4);
            VVarEntree{i_,j_}(:,1) =  (1:numel(param)).';
            VVarEntree{i_,j_}(:,2) = pdf*ones(numel(param),1);
            VVarEntree{i_,j_}(:,3) = param(:).*(1.0-cv);
            VVarEntree{i_,j_}(:,4) = param(:).*(1.0+cv);
        end
    end
    
    %% *OUTPUT*
    varargout{1} = trs;
    varargout{2} = VVarEntree;

    return
end
