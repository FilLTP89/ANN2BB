%% *Train ANN network onto synthetics*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ann2hbs_train_: function to train ANN network onto hybrid LF/HF
% synthetics
%% INPUT:
%%
% * _hbs (hybrid synthetics structure)_
%%
% * _hbs.mon = monitor structure_
% * hbs.mon.pt = path to monitor files    (string)
% * hbs.mon.fn = monitor metadata filename(string)
% * hbs.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * hbs.mon.id = monitor identity         (integer)
% * hbs.mon.dep = epicentral distance     (real vector)
% * hbs.mon.stn = monitor names           (string vector)
% * hbs.mon.na = number of monitors       (integer)
% * hbs.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * hbs.mon.nr = number of records        (integer)
% * hbs.mon.cp = motion component         (integer: 1,2,3)
% * hbs.mon.nc = number of components     (integer)
% * hbs.mon.na  = number of accelerograms to be generated
% * hbs.mon.nc  = number of motion components (integer)
% * hbs.mon.cp(mon.nc,1)  = motion components (string vector)
% * hbs.mon.dtm(mon.na,1) = time-steps        (real vector)
% * hbs.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
% * hbs.mon.ntm(mon.na,1) = time step number  (real vector)
% * hbs.mon.nT            = natural period number (real)
% * hbs.mon.vTn(mon.nT,1) = natural period vector (real vector)
% * hbs.mon.zeta          = damping               (real)
% * hbs.mon.dfr(mon.na,1) = frequency step   (real vector)
% * hbs.mon.nfr(mon.na,1) = frequency number (real vector)
% * hbs.mon.nNy(mon.na,1) = Nyquist index    (real vector)
% * hbs.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%%
% * _hbs.syn(mtd.na,1)    = structure vectors (cell vector)_
% * hbs.syn{i}.tha.x      = x-acceleration    (real vector)
% * hbs.syn{i}.tha.y      = y-acceleration    (real vector)
% * hbs.syn{i}.tha.z      = z-acceleration    (real vector)
% * hbs.syn{i}.rsd.x = x-displacement spectrum
% * hbs.syn{i}.rsd.y = y-displacement spectrum
% * hbs.syn{i}.rsd.z = z-displacement spectrum
% * hbs.syn{i}.psa.x = x-pseudo-acceleration spectrum
% * hbs.syn{i}.psa.y = y-pseudo-acceleration spectrum
% * hbs.syn{i}.psa.z = z-pseudo-acceleration spectrum
% * hbs.syn{i}.fsa.x = x-fourier spectrum
% * hbs.syn{i}.fsa.y = y-fourier spectrum
% * hbs.syn{i}.fsa.z = z-fourier spectrum
%%
% * _ann (trained Artificial Neural Network (ANN) structure)_
% * ann.hrz = trained ANN for horizontal components
% * ann.vrt = trained ANN for vertical   components
%% OUTPUT:
% * _trs (output structure)_

function [varargout] = ann2hbs_train(varargin)
    
    %% SET-UP
    hbs = varargin{1};
    ann = varargin{2};
    trs.mon = hbs.mon;
    
    %%
    % _input periods_
    
    for j_ = 1:hbs.mon.nc
        for k_ = 1:ann.(hbs.mon.cp{j_}).inp.nT
            idx = find(abs(hbs.mon.vTn-ann.(hbs.mon.cp{j_}).inp.vTn(k_))<1e-8);
            if ~isempty(idx)
                iid(k_) = idx;
            end
        end
    end
    %%
    % _target natural periods_
    for j_ = 1:hbs.mon.nc
        for k_ = 1:ann.(hbs.mon.cp{j_}).tar.nT
            idx = find(abs(hbs.mon.vTn-ann.(hbs.mon.cp{j_}).tar.vTn(k_))<1e-8);
            if ~isempty(idx)
                tid(k_) = idx;
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
        inn{j_} = [psa(:,iid),pgv,pgd];
    end
    
    %% TRAIN NETWORK
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_});
        for i_ = 1:hbs.mon.na
            try
                out = 10.^(sim(ann.(hbs.mon.cp{j_}).net,inp(i_,:)'));
                trs.syn{i_}.vTn.(hbs.mon.cp{j_}) = ...
                    [ann.(hbs.mon.cp{j_}).tar.vTn;ann.(hbs.mon.cp{j_}).inp.vTn];
                trs.syn{i_}.nTn.(hbs.mon.cp{j_}) = numel(trs.syn{i_}.vTn.(hbs.mon.cp{j_}));
                trs.syn{i_}.psa.(hbs.mon.cp{j_}) = [out(:);inn{j_}(i_,1:end-2)'];
            catch
                keyboard
            end
        end
    end

    %% OUTPUT
    varargout{1} = trs;
    return
end
