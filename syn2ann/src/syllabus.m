%% *BOREHOLE*
%
% _bhr (borehole data structure)_
%
% * bhr.pt = path to record files (string)
%
% * bhr.st = station id      (string)
%
% * bhr.ns = numel(bhr.st) (integer)
%
% * bhr.ev = recorded events (string vector)
%
% * bhr.ne = numel(bhr.ev) (integer)
%
% * bhr.tp = type of record file (string vector)
%
% * bhr.lb = record's labels (string vector)
%
% * bhr.dv = borehole device list (string vector)
%
% * bhr.cd = borehole control devices (integer vector)
%
% * bhr.nd = numel(bhr.dv) (integer)
%
% * bhr.cp = motion directions (string vector)
%
% * bhr.nc = numel(bhr.cp) (integer)
%
% * bhr.rc = motion component (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
%
% * bhr.nr = numel(bhr.rc) (integer)
%
% * bhe.id = recording device identity (string)
%
% * bhr.nm = record filenames (string vector)

%% *MONITOR*
%
% _mon (monitor structure)_
%
% * mon.pt = path to monitor files    (string)
%
% * mon.fn = monitor metadata filename(string)
%
% * mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
%
% * mon.id = monitor identity         (integer)
%
% * mon.na = number of monitors       (integer)
%
% * mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
%
% * mon.nr = number of records        (integer)
%
% * mon.cp = motion component         (string vector)
%
% * mon.nc = number of components     (integer)
%
% * mon.dep = epicentral distance     (real vector)
%
% * mon.st = monitor names            (string vector)

%% *RECORDED DATA*
%
% * _rec    = structure of numerical simulations_
%
% * _rec.ev(bhr.ne,1)    = array of structure of borehole data_
%
% * rec.ev.bhr.id = station id      (string)
%
% * rec.ev.bhr.ev = recorded events (string vector)
%
% * rec.ev.bhr.ne  = number of recorded events(integer)
%
% * rec.ev.bhr.tp  = type of record file (string vector)
%
% * rec.ev.bhr.lb = record's labels (string vector)
%
% * rec.ev.bhr.dv = borehole device list (string vector)
%
% * _rec.syn(bhr.na,1)     = structure vectors (cell vector)_
%
% * rec.syn{i}.tha.x      = x-acceleration    (real vector)
% * rec.syn{i}.tha.y      = y-acceleration    (real vector)
% * rec.syn{i}.tha.z      = z-acceleration    (real vector)
% * rec.syn{i}.thv.x      = x-velocity        (real vector)
% * rec.syn{i}.thv.y      = y-velocity        (real vector)
% * rec.syn{i}.thv.z      = z-velocity        (real vector)
% * rec.syn{i}.thd.x      = x-displacement    (real vector)
% * rec.syn{i}.thd.y      = y-displacement    (real vector)
% * rec.syn{i}.thd.z      = z-displacement    (real vector)
%
% * rec.syn{i}.pga.x(1)   = x-time-pga        (real vector)
% * rec.syn{i}.pga.x(2)   = x-pga             (real vector)
% * rec.syn{i}.pga.y(1)   = y-time-pga        (real vector)
% * rec.syn{i}.pga.y(2)   = y-pga             (real vector)
% * rec.syn{i}.pga.z(1)   = z-time-pga        (real vector)
% * rec.syn{i}.pga.z(2)   = z-pga             (real vector)
% * rec.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
% * rec.syn{i}.pgv.x(2)   = x-pgv             (real vector)
% * rec.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
% * rec.syn{i}.pgv.y(2)   = y-pgv             (real vector)
% * rec.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
% * rec.syn{i}.pgv.z(2)   = z-pgv             (real vector)
% * rec.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
% * rec.syn{i}.pgd.x(2)   = x-pgd             (real vector)
% * rec.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
% * rec.syn{i}.pgd.y(2)   = y-pgd             (real vector)
% * rec.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
% * rec.syn{i}.pgd.z(2)   = z-pgd             (real vector)
%
%% NUMERICAL SIMULATION
% _nss (structure of numerical simulations)_
%
% * _nss.mon    = structure of monitor data_
%
% * nss.mon.pt = path to monitor files    (string)
% * nss.mon.fn = monitor metadata filename(string)
% * nss.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * nss.mon.id = monitor identity         (integer)
% * nss.mon.na = number of monitors       (integer)
% * nss.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * nss.mon.nr = number of records        (integer)
% * nss.mon.cp = motion component         (string vector)
% * nss.mon.nc = number of components     (integer)
% * nss.mon.dep = epicentral distance     (real vector)
% * nss.mon.st = monitor names            (string vector)
% * nss.mon.dtm(mon.na,1) = time-steps        (real vector)
% * nss.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
% * nss.mon.ntm(mon.na,1) = time step number  (real vector)
% * nss.mon.nT            = natural period number (real)
% * nss.mon.vTn(mon.nT,1) = natural period vector (real vector)
% * nss.mon.zeta          = damping               (real)
% * nss.mon.dfr(mon.na,1) = frequency-step   (real vector)
% * nss.mon.nfr(mon.na,1) = frequency number (real vector)
% * nss.mon.nNy(mon.na,1) = Nyquist index    (real vector)
% * nss.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%..........................................................................
% * _nss.syn(mon.na,1)     = structure vectors (cell vector)_
%..........................................................................
% * nss.syn{i}.tha.x      = x-acceleration    (real vector)
% * nss.syn{i}.tha.y      = y-acceleration    (real vector)
% * nss.syn{i}.tha.z      = z-acceleration    (real vector)
% * nss.syn{i}.thv.x      = x-velocity        (real vector)
% * nss.syn{i}.thv.y      = y-velocity        (real vector)
% * nss.syn{i}.thv.z      = z-velocity        (real vector)
% * nss.syn{i}.thd.x      = x-displacement    (real vector)
% * nss.syn{i}.thd.y      = y-displacement    (real vector)
% * nss.syn{i}.thd.z      = z-displacement    (real vector)
%..........................................................................
% * nss.syn{i}.pga.x(1)   = x-time-pga        (real vector)
% * nss.syn{i}.pga.x(2)   = x-pga             (real vector)
% * nss.syn{i}.pga.y(1)   = y-time-pga        (real vector)
% * nss.syn{i}.pga.y(2)   = y-pga             (real vector)
% * nss.syn{i}.pga.z(1)   = z-time-pga        (real vector)
% * nss.syn{i}.pga.z(2)   = z-pga             (real vector)
% * nss.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
% * nss.syn{i}.pgv.x(2)   = x-pgv             (real vector)
% * nss.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
% * nss.syn{i}.pgv.y(2)   = y-pgv             (real vector)
% * nss.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
% * nss.syn{i}.pgv.z(2)   = z-pgv             (real vector)
% * nss.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
% * nss.syn{i}.pgd.x(2)   = x-pgd             (real vector)
% * nss.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
% * nss.syn{i}.pgd.y(2)   = y-pgd             (real vector)
% * nss.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
% * nss.syn{i}.pgd.z(2)   = z-pgd             (real vector)
%..........................................................................
% * nss.syn{i}.rsd.x      = x-displacement response spectrum (real vectors)
% * nss.syn{i}.rsd.y      = y-displacement response spectrum (real vectors)
% * nss.syn{i}.rsd.z      = z-displacement response spectrum (real vectors)
%..........................................................................
% * nss.syn{i}.psa.x      = x-pseudo-acceleration spectrum
% * nss.syn{i}.psa.y      = y-pseudo-acceleration spectrum
% * nss.syn{i}.psa.z      = z-pseudo-acceleration spectrum
%..........................................................................
% * nss.syn{i}.fsa.x      = x-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsa.y      = y-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsa.z      = z-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsv.x      = x-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsv.y      = y-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsv.z      = z-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsd.x      = x-displacement fourier spectrum (real vector)
% * nss.syn{i}.fsd.y      = y-displacement fourier spectrum (real vector)
% * nss.syn{i}.fsd.z      = z-displacement fourier spectrum (real vector)
%==========================================================================

%% SYNTHETICS SABETTA&PUGLIESE 1996
% * _sps (structure of sp synthetics)_
%..........................................................................
% * _sps.mtd     = metadata structure_
%..........................................................................
% * sps.mtd.mw  = moment magnitude (real)
% * sps.mtd.dtm = time step        (real)
% * sps.mtd.scc = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)
% * sps.mtd.sst = st. dev. of GMPE (0=median value,1=84th percentile)
% * sps.mtd.scl = scale factor [1=cm/s/s]
% * sps.mtd.ivd = flag for output in velocity and displacement
%..........................................................................
% * _sps.mon = monitor structure_
%..........................................................................
% * sps.mon.pt = path to monitor files    (string)
% * sps.mon.fn = monitor metadata filename(string)
% * sps.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * sps.mon.id = monitor identity         (integer)
% * sps.mon.na = number of monitors       (integer)
% * sps.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * sps.mon.nr = number of records        (integer)
% * sps.mon.cp = motion component         (integer: 1,2,3)
% * sps.mon.nc = number of components     (integer)
% * sps.mon.na  = number of accelerograms to be generated
% * sps.mon.nc  = number of motion components (integer)
% * sps.mon.dep = epicentral distance     (real vector)
% * sps.mon.st = monitor names            (string vector)
% * sps.mon.cp(mon.nc,1)  = motion components (string vector)
% * sps.mon.dtm(mon.na,1) = time-steps        (real vector)
% * sps.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
% * sps.mon.ntm(mon.na,1) = time step number  (real vector)
% * sps.mon.nT            = natural period number (real)
% * sps.mon.vTn(mon.nT,1) = natural period vector (real vector)
% * sps.mon.zeta          = damping               (real)
% * sps.mon.dfr(mon.na,1) = frequency step   (real vector)
% * sps.mon.nfr(mon.na,1) = frequency number (real vector)
% * sps.mon.nNy(mon.na,1) = Nyquist index    (real vector)
% * sps.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%..........................................................................
% * _sps.syn(mtd.na,1)     = structure vectors (cell vector)_
%..........................................................................
% * sps.syn{i}.tha.x      = x-acceleration    (real vector)
% * sps.syn{i}.tha.y      = y-acceleration    (real vector)
% * sps.syn{i}.tha.z      = z-acceleration    (real vector)
% * sps.syn{i}.thv.x      = x-velocity        (real vector)
% * sps.syn{i}.thv.y      = y-velocity        (real vector)
% * sps.syn{i}.thv.z      = z-velocity        (real vector)
% * sps.syn{i}.thd.x      = x-displacement    (real vector)
% * sps.syn{i}.thd.y      = y-displacement    (real vector)
% * sps.syn{i}.thd.z      = z-displacement    (real vector)
%..........................................................................
% * sps.syn{i}.pga.x(1)   = x-time-pga        (real vector)
% * sps.syn{i}.pga.x(2)   = x-pga             (real vector)
% * sps.syn{i}.pga.y(1)   = y-time-pga        (real vector)
% * sps.syn{i}.pga.y(2)   = y-pga             (real vector)
% * sps.syn{i}.pga.z(1)   = z-time-pga        (real vector)
% * sps.syn{i}.pga.z(2)   = z-pga             (real vector)
% * sps.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
% * sps.syn{i}.pgv.x(2)   = x-pgv             (real vector)
% * sps.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
% * sps.syn{i}.pgv.y(2)   = y-pgv             (real vector)
% * sps.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
% * sps.syn{i}.pgv.z(2)   = z-pgv             (real vector)
% * sps.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
% * sps.syn{i}.pgd.x(2)   = x-pgd             (real vector)
% * sps.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
% * sps.syn{i}.pgd.y(2)   = y-pgd             (real vector)
% * sps.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
% * sps.syn{i}.pgd.z(2)   = z-pgd             (real vector)
%..........................................................................
% * sps.syn{i}.rsd.x = x-displacement spectrum
% * sps.syn{i}.rsd.y = y-displacement spectrum
% * sps.syn{i}.rsd.z = z-displacement spectrum
%..........................................................................
% * sps.syn{i}.psa.x = x-pseudo-acceleration spectrum
% * sps.syn{i}.psa.y = y-pseudo-acceleration spectrum
% * sps.syn{i}.psa.z = z-pseudo-acceleration spectrum
%..........................................................................
% * sps.syn{i}.fsa.x = x-fourier spectrum
% * sps.syn{i}.fsa.y = y-fourier spectrum
% * sps.syn{i}.fsa.z = z-fourier spectrum
%==========================================================================

%% HYBRID SYNTHETICS
% _hbs (hybrid synthetics structure)_
%..........................................................................
% * _hbs.mon = monitor structure_
%..........................................................................
% * hbs.mon.pt = path to monitor files    (string)
% * hbs.mon.fn = monitor metadata filename(string)
% * hbs.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * hbs.mon.id = monitor identity         (integer)
% * hbs.mon.dep = epicentral distance     (real vector)
% * hbs.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * hbs.mon.nr = number of records        (integer)
% * hbs.mon.cp = motion component         (integer: 1,2,3)
% * hbs.mon.nc = number of components     (integer)
% * hbs.mon.na  = number of accelerograms to be generated
% * hbs.mon.nc  = number of motion components (integer)
% * hbs.mon.cp(mon.nc,1)  = motion components (string vector)
% * hbs.mon.st = monitor names           (string vector)
% * hbs.mon.na = number of monitors       (integer)
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
%..........................................................................
% * _hbs.syn(mtd.na,1)     = structure vectors (cell vector)_
%..........................................................................
% * hbs.syn{i}.tha.x      = x-acceleration    (real vector)
% * hbs.syn{i}.tha.y      = y-acceleration    (real vector)
% * hbs.syn{i}.tha.z      = z-acceleration    (real vector)
% * hbs.syn{i}.thv.x      = x-velocity        (real vector)
% * hbs.syn{i}.thv.y      = y-velocity        (real vector)
% * hbs.syn{i}.thv.z      = z-velocity        (real vector)
% * hbs.syn{i}.thd.x      = x-displacement    (real vector)
% * hbs.syn{i}.thd.y      = y-displacement    (real vector)
% * hbs.syn{i}.thd.z      = z-displacement    (real vector)
%..........................................................................
% * hbs.syn{i}.pga.x(1)   = x-time-pga        (real vector)
% * hbs.syn{i}.pga.x(2)   = x-pga             (real vector)
% * hbs.syn{i}.pga.y(1)   = y-time-pga        (real vector)
% * hbs.syn{i}.pga.y(2)   = y-pga             (real vector)
% * hbs.syn{i}.pga.z(1)   = z-time-pga        (real vector)
% * hbs.syn{i}.pga.z(2)   = z-pga             (real vector)
% * hbs.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
% * hbs.syn{i}.pgv.x(2)   = x-pgv             (real vector)
% * hbs.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
% * hbs.syn{i}.pgv.y(2)   = y-pgv             (real vector)
% * hbs.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
% * hbs.syn{i}.pgv.z(2)   = z-pgv             (real vector)
% * hbs.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
% * hbs.syn{i}.pgd.x(2)   = x-pgd             (real vector)
% * hbs.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
% * hbs.syn{i}.pgd.y(2)   = y-pgd             (real vector)
% * hbs.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
% * hbs.syn{i}.pgd.z(2)   = z-pgd             (real vector)
% * hbs.syn{i}.rsd.x = x-displacement spectrum
% * hbs.syn{i}.rsd.y = y-displacement spectrum
% * hbs.syn{i}.rsd.z = z-displacement spectrum
% * hbs.syn{i}.psa.x = x-pseudo-acceleration spectrum
% * hbs.syn{i}.psa.y = y-pseudo-acceleration spectrum
% * hbs.syn{i}.psa.z = z-pseudo-acceleration spectrum
% * hbs.syn{i}.fsa.x = x-fourier spectrum
% * hbs.syn{i}.fsa.y = y-fourier spectrum
% * hbs.syn{i}.fsa.z = z-fourier spectrum
%..........................................................................

