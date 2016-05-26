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
%         trs.mon.pt = path to monitor files    (string)
%         trs.mon.fn = monitor metadata filename(string)
%         trs.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
%         trs.mon.id = monitor identity         (integer)
%         trs.mon.dep = epicentral distance     (real vector)
%         trs.mon.stn = monitor names           (string vector)
%         trs.mon.na = number of monitors       (integer)
%         trs.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
%         trs.mon.nr = number of records        (integer)
%         trs.mon.cp = motion component         (integer: 1,2,3)
%         trs.mon.nc = number of components     (integer)
%         trs.mon.na  = number of accelerograms to be generated
%         trs.mon.nc  = number of motion components (integer)
%         trs.mon.cp(mon.nc,1)  = motion components (string vector)
%         trs.mon.dtm(mon.na,1) = time-steps        (real vector)
%         trs.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
%         trs.mon.ntm(mon.na,1) = time step number  (real vector)
%         trs.mon.nT            = natural period number (real)
%         trs.mon.vTn(mon.nT,1) = natural period vector (real vector)
%         trs.mon.zeta          = damping               (real)
%         trs.mon.dfr(mon.na,1) = frequency step   (real vector)
%         trs.mon.nfr(mon.na,1) = frequency number (real vector)
%         trs.mon.nNy(mon.na,1) = Nyquist index    (real vector)
%         trs.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%%
%         _hbs.syn(mtd.na,1)     = structure vectors (cell vector)_
%         trs.syn{i}.tha.x      = x-acceleration    (real vector)
%         trs.syn{i}.tha.y      = y-acceleration    (real vector)
%         trs.syn{i}.tha.z      = z-acceleration    (real vector)
%         trs.syn{i}.rsd.x = x-displacement spectrum
%         trs.syn{i}.rsd.y = y-displacement spectrum
%         trs.syn{i}.rsd.z = z-displacement spectrum
%         trs.syn{i}.psa.x = x-pseudo-acceleration spectrum
%         trs.syn{i}.psa.y = y-pseudo-acceleration spectrum
%         trs.syn{i}.psa.z = z-pseudo-acceleration spectrum
%         trs.syn{i}.fsa.x = x-fourier spectrum
%         trs.syn{i}.fsa.y = y-fourier spectrum
%         trs.syn{i}.fsa.z = z-fourier spectrum
%         _ann_
%% OUTPUT:
%
%% N.B.
function [varargout] = synthetics2ann_spectral_matching(varargin)
    
    %% SET-UP
    hbs = varargin{1};
    trs = varargin{2};
    
    %% SPECTRAL SCALING
    out.mon = trs.mon;
    
    for j_ = 1:trs.mon.nc
        [out.mon.dtm,tha.(hbs.mon.cp{j_}),thv.(hbs.mon.cp{j_}),thd.(hbs.mon.cp{j_})] =...
            arrayfun(@(A,B,C) spectral_scaling(A,B{:}.tha.(hbs.mon.cp{j_}),...
            C{:}.psa.(hbs.mon.cp{j_}),C{:}.vTn.(hbs.mon.cp{j_})),...
            hbs.mon.dtm(:),hbs.syn(:),trs.syn(:),'UniformOutput',0);
        out.syn{:}.tha=tha;
        out.syn{:}.thv=thv;
        out.syn{:}.thd=thd;
    end
    %% OUTPUT
    varargout{1} = out;
    return
end