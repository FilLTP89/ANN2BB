%% *Spectral Matching*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _synthetics2ann_spectral_matching_: function to match hybrid synthetic
% spectra with trained ANN results
%% INPUT:
%         _hbs    = structure of hybrid synthetics_
%%
%         _hbs.mon = monitor structure_
%         hbs.mon.pt = path to monitor files    (string)
%         hbs.mon.fn = monitor metadata filename(string)
%         hbs.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
%         hbs.mon.id = monitor identity         (integer)
%         hbs.mon.dep = epicentral distance     (real vector)
%         hbs.mon.stn = monitor names           (string vector)
%         hbs.mon.na = number of monitors       (integer)
%         hbs.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
%         hbs.mon.nr = number of records        (integer)
%         hbs.mon.cp = motion component         (integer: 1,2,3)
%         hbs.mon.nc = number of components     (integer)
%         hbs.mon.na  = number of accelerograms to be generated
%         hbs.mon.nc  = number of motion components (integer)
%         hbs.mon.cp(mon.nc,1)  = motion components (string vector)
%         hbs.mon.dtm(mon.na,1) = time-steps        (real vector)
%         hbs.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
%         hbs.mon.ntm(mon.na,1) = time step number  (real vector)
%         hbs.mon.nT            = natural period number (real)
%         hbs.mon.vTn(mon.nT,1) = natural period vector (real vector)
%         hbs.mon.zeta          = damping               (real)
%         hbs.mon.dfr(mon.na,1) = frequency step   (real vector)
%         hbs.mon.nfr(mon.na,1) = frequency number (real vector)
%         hbs.mon.nNy(mon.na,1) = Nyquist index    (real vector)
%         hbs.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%%
%         _hbs.syn(mtd.na,1)     = structure vectors (cell vector)_
%         hbs.syn{i}.tha.x      = x-acceleration    (real vector)
%         hbs.syn{i}.tha.y      = y-acceleration    (real vector)
%         hbs.syn{i}.tha.z      = z-acceleration    (real vector)
%         hbs.syn{i}.rsd.x = x-displacement spectrum
%         hbs.syn{i}.rsd.y = y-displacement spectrum
%         hbs.syn{i}.rsd.z = z-displacement spectrum
%         hbs.syn{i}.psa.x = x-pseudo-acceleration spectrum
%         hbs.syn{i}.psa.y = y-pseudo-acceleration spectrum
%         hbs.syn{i}.psa.z = z-pseudo-acceleration spectrum
%         hbs.syn{i}.fsa.x = x-fourier spectrum
%         hbs.syn{i}.fsa.y = y-fourier spectrum
%         hbs.syn{i}.fsa.z = z-fourier spectrum
%         _ann_
%% OUTPUT:
%
%% N.B.
function [varargout] = synthetics2ann_spectral_matching(varargin)
    
    %% SET-UP
    hbs = varargin{1};
    trs = varargin{2};
    
    %% SPECTRAL SCALING
    out.mon = hbs.mon;
    for j_ = 1:hbs.mon.nc
        for i_ = 1:hbs.mon.na
            eval(sprintf(['[out.mon.dtm(i_),out.syn{i_}.tha.%s,out.syn{i_}.thv.%s,'...
                'out.syn{i_}.thd.%s,out.mon.dfr(i_),out.syn{i_}.psa.%s] = ',...
                'spectral_scaling(hbs.mon.dtm(i_),hbs.syn{i_}.tha.%s,',...
                'trs{i_,j_}.psa,trs{i_,j_}.vTn);'],hbs.mon.cp{j_},hbs.mon.cp{j_},...
                hbs.mon.cp{j_},hbs.mon.cp{j_},hbs.mon.cp{j_}));
        end
    end
    %% OUTPUT
    varargout{1} = out;
    return
end