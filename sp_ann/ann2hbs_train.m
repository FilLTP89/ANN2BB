function [varargout] = ann2hbs_train(varargin)
    %===============
    % Train ANN network onto synthetics
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % ann2hbs_train: function to train ANN network onto hybrid LF/HF
    % synthetics
    % INPUT:  ann (trained artificial neural network structure)
    % OUTPUT: onn (output structure)
    % N.B. Need for
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    ann = varargin{1};
    hbs = varargin{2};
    inn =  cell(hbs.mon.nc,ann.tar.nT);
    onn =  cell(hbs.mon.nc,ann.tar.nT);
    %======================================================================
    %======================================================================
    % ANN TRAINING INPUT
    %======================================================================
    %----------------------------------------------------------------------
    % input periods
    %----------------------------------------------------------------------
    for k_ = 1:ann.inp.nT
        idx = find(abs(hbs.mon.vTn-ann.inp.vTn(k_))<1e-8);
        if ~isempty(idx)
            iid(k_) = idx;
        end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % target natural periods
    %----------------------------------------------------------------------
    for k_ = 1:ann.tar.nT
        idx = find(abs(hbs.mon.vTn-ann.tar.vTn(k_))<1e-8);
        if ~isempty(idx)
            tid(k_) = idx;
        end
    end
    %----------------------------------------------------------------------
    for j_ = 1:hbs.mon.nc
        eval(sprintf(['psa = reshape(cell2mat(arrayfun(@(v) v{:}.psa.%s,hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.nT,hbs.mon.na)'';'],hbs.mon.cp{j_}));
        eval(sprintf(['pgv = reshape(cell2mat(arrayfun(@(v) v{:}.pgv.%s(2),hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.na,1);'],hbs.mon.cp{j_}));
        eval(sprintf(['pgd = reshape(cell2mat(arrayfun(@(v) v{:}.pgd.%s(2),hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.na,1);'],hbs.mon.cp{j_}));
        inn{j_} = [psa(:,iid),pgv,pgd];
    end
    %======================================================================
    %======================================================================
    % TRAIN NETWORK
    %======================================================================
    
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_});
        for i_ = 1:hbs.mon.na
            try
                out = 10.^(sim(ann.net,inp(i_,:)'));
                onn{j_,i_}.vTn = [out(:);ann.tar.vTn(:)];
                onn{j_,i_}.nTn = numel(onn{j_,i_}.vTn);
                onn{j_,i_}.psa = [out(:);inn{j_}(i_,1:end-2)'];
            catch
                keyboard
            end
        end
    end
    
    %======================================================================
    varargout{1} = onn;
    return
end
