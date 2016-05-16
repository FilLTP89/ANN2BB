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
% * _hbs.syn(mtd.na,1)     = structure vectors (cell vector)_
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
% * _ann (trained artificial neural network structure)_
%% OUTPUT:
% * _trs (output structure)_

function [varargout] = ann2hbs_train(varargin)
    
    %% SET-UP
    hbs = varargin{1};
    ann = varargin{2};
    inn =  cell(hbs.mon.nc,ann.tar.nT);
    trs =  cell(hbs.mon.nc,ann.tar.nT);
    %%
    % _input periods_
    for k_ = 1:ann.inp.nT
        idx = find(abs(hbs.mon.vTn-ann.inp.vTn(k_))<1e-8);
        if ~isempty(idx)
            iid(k_) = idx;
        end
    end
    %%
    % _target natural periods_
    for k_ = 1:ann.tar.nT
        idx = find(abs(hbs.mon.vTn-ann.tar.vTn(k_))<1e-8);
        if ~isempty(idx)
            tid(k_) = idx;
        end
    end
    
    for j_ = 1:hbs.mon.nc
        eval(sprintf(['psa = reshape(cell2mat(arrayfun(@(v) v{:}.psa.%s,hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.nT,hbs.mon.na)'';'],hbs.mon.cp{j_}));
        eval(sprintf(['pgv = reshape(cell2mat(arrayfun(@(v) v{:}.pgv.%s(2),hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.na,1);'],hbs.mon.cp{j_}));
        eval(sprintf(['pgd = reshape(cell2mat(arrayfun(@(v) v{:}.pgd.%s(2),hbs.syn,',...
            '''UniformOutput'',0)),hbs.mon.na,1);'],hbs.mon.cp{j_}));
        inn{j_} = [psa(:,iid),pgv,pgd];
    end
    %% TRAIN NETWORK
    
    for j_ = 1:hbs.mon.nc
        inp = log10(inn{j_});
        for i_ = 1:hbs.mon.na
            try
                out = 10.^(sim(ann.net,inp(i_,:)'));
                trs{j_,i_}.vTn = [ann.tar.vTn;ann.inp.vTn];
                trs{j_,i_}.nTn = numel(trs{j_,i_}.vTn);
                trs{j_,i_}.psa = [out(:);inn{j_}(i_,1:end-2)'];
            catch
                keyboard
            end
        end
    end
    %% OUTPUT
    varargout{1} = trs;
    return
end
