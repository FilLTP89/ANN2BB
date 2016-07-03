%% *Train ANN network onto synthetics*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ann2hbs_train_: function to train ANN network onto hybrid LF/HF
% synthetics
%% INPUT:
% * _hbs (hybrid synthetics structure)_
% * _ann (trained Artificial Neural Network (ANN) structure)_
%% OUTPUT:
% * _trs (trained/simulated ann structure)_
%% N.B. 
% Need for _sim.m_
function [varargout] = ann2hbs_train(varargin)
    %% *SET-UP*
    hbs = varargin{1};
    ann = varargin{2};
    % _input periods_
    for j_ = 1:hbs.mon.nc
        for k_ = 1:ann.(hbs.mon.cp{j_}).inp.nT
            idx = find(abs(hbs.mon.vTn-ann.(hbs.mon.cp{j_}).inp.vTn(k_))<1e-8);
            if ~isempty(idx)
                trs.(hbs.mon.cp{j_}).iid(k_) = idx;
            end
        end
    end
    %%
    % _target natural periods_
    for j_ = 1:hbs.mon.nc
        for k_ = 1:ann.(hbs.mon.cp{j_}).tar.nT
            idx = find(abs(hbs.mon.vTn-ann.(hbs.mon.cp{j_}).tar.vTn(k_))<1e-8);
            if ~isempty(idx)
                trs.(hbs.mon.cp{j_}).tid(k_) = idx;
            end
        end
    end
    
    inn = cell(hbs.mon.nc,1);
    for j_ = 1:hbs.mon.nc
        psa = reshape(cell2mat(cellfun(@(v) v.psa.(hbs.mon.cp{j_}),hbs.syn,...
            'UniformOutput',0)),hbs.mon.na,hbs.mon.nT);
        pgv = reshape(cell2mat(cellfun(@(v) v.pgv.(hbs.mon.cp{j_})(2),hbs.syn,...
            'UniformOutput',0)),hbs.mon.na,1);
        pgd = reshape(cell2mat(cellfun(@(v) v.pgd.(hbs.mon.cp{j_})(2),hbs.syn,...
            'UniformOutput',0)),hbs.mon.na,1);
        inn{j_} = [psa(:,trs.(hbs.mon.cp{j_}).iid),abs(pgv),abs(pgd)];
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
                trs.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}) = [out(:);inn{j_}(i_,1:end-2)'];
            catch
                keyboard
            end
        end
    end
    
    %% OUTPUT
    varargout{1} = trs;
    return
end