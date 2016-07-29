%% *Train ANN network onto synthetics*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ann2hbs_train_: function to train ANN network onto hybrid LF/HF
% synthetics, trained PSA(LF) values
%% INPUT:
% * _hbs (hybrid synthetics structure)_
% * _ann (trained Artificial Neural Network (ANN) structure)_
%% OUTPUT:
% * _trs (trained/simulated ann structure)_
%% N.B. 
% Need for _sim.m_
function [varargout] = ann2hbs_train_justPSA(varargin)
    %% *SET-UP*
    hbs = varargin{1};
    ann = varargin{2};
    
    % _input periods_
    for j_ = 1:hbs.mon.nc
        cpp = (hbs.mon.cp{j_});
        for k_ = 1:ann.(cpp).inp.nT
            idx = find(abs(hbs.mon.vTn-ann.(cpp).inp.vTn(k_))<1e-8);
            if ~isempty(idx)
                trs.(cpp).iid(k_) = idx;
            end
        end
    end
    
    % _target natural periods_
    for j_ = 1:hbs.mon.nc
        cpp = (hbs.mon.cp{j_});
        for k_ = 1:ann.(cpp).tar.nT
            idx = find(abs(hbs.mon.vTn-ann.(cpp).tar.vTn(k_))<1e-8);
            if ~isempty(idx)
                trs.(cpp).tid(k_) = idx;
            end
        end
    end
    
    inn = cell(hbs.mon.nc,1);
    
    for j_ = 1:hbs.mon.nc
        psa = -ones(hbs.mon.na,hbs.mon.nT);
        pgv = -ones(hbs.mon.na,1);
        pgd = -ones(hbs.mon.na,1);
        cpp = hbs.mon.cp{j_};
        for i_ = 1:hbs.mon.na
            psa(i_,:) = hbs.syn{i_}.psa.(cpp)(:)';
        end
        inn{j_} = [psa(:,trs.(cpp).iid)];
    end
    
    %% TRAIN NETWORK
    trs.mon.cp = hbs.mon.cp;
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_}.*100);
        trs.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp;
        trs.(hbs.mon.cp{j_}).mon.vTn = ...
            [ann.(hbs.mon.cp{j_}).tar.vTn;ann.(hbs.mon.cp{j_}).inp.vTn];
        trs.(hbs.mon.cp{j_}).mon.nT = numel(trs.(hbs.mon.cp{j_}).mon.vTn);
        for i_ = 1:hbs.mon.na
            try
                out = 10.^(sim(ann.(hbs.mon.cp{j_}).net,inp(i_,:)'));
                out = out./100;
                trs.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}) = [out(:);inn{j_}(i_,1:end-1)'];
            catch
                keyboard
            end
        end
    end
    
    %% OUTPUT
    varargout{1} = trs;
    return
end
