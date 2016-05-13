function [varargout] = synthetics2ann_spectral_matching(varargin)
    %===============
    % Spectral Matching
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % synthetics2ann_spectral_matching: function to match psa spectra
    % between synthetics and Artificial Neural Network
    % OUTPUT: hbs    = structure of hybrid synthetics
    %         hbs.mon = monitor structure
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
    %         hbs.syn(mtd.na,1)     = structure vectors (cell vector)
    %         hbs.syn{i}.tha.x      = x-acceleration    (real vector)
    %         hbs.syn{i}.tha.y      = y-acceleration    (real vector)
    %         hbs.syn{i}.tha.z      = z-acceleration    (real vector)
    %         ann
    % OUTPUT: hbs    = structure of hybrid synthetics
    %         hbs.mon = monitor structure
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
    %         hbs.syn(mtd.na,1)     = structure vectors (cell vector)
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
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    hbs = varargin{1};
    ann = varargin{2};
    
    varargout{1} = hbs;
    return
end