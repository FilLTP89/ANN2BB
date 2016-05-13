function [varargout] = sp_spectra(varargin)
    %===============
    % Compute response/fourier spectra
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % sp_spectra: function to compute response/fourier spectra for a set
    % of synthetics records generated with Sabetta and Pugliese method
    % INPUT:  sps    = structure of sp synthetics
    %         sps.mtd     = metadata structure
    %         sps.mtd.mw  = moment magnitude (real)
    %         sps.mtd.scc = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)
    %         sps.mtd.sst = st. dev. of GMPE (0=median value,1=84th percentile)
    %         sps.mtd.scl = scale factor [1=cm/s/s]
    %         sps.mtd.ivd = flag for output in velocity and displacement
    %         sps.mon = monitor structure
    %         sps.mon.pt = path to monitor files    (string)
    %         sps.mon.fn = monitor metadata filename(string)
    %         sps.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
    %         sps.mon.id = monitor identity         (integer)
    %         sps.mon.dep = epicentral distance     (real vector)
    %         sps.mon.stn = monitor names           (string vector)
    %         sps.mon.na = number of monitors       (integer)
    %         sps.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
    %         sps.mon.nr = number of records        (integer)
    %         sps.mon.cp = motion component         (integer: 1,2,3)
    %         sps.mon.nc = number of components     (integer)
    %         sps.mon.na  = number of accelerograms to be generated
    %         sps.mon.nc  = number of motion components (integer)
    %         sps.mon.cp(mon.nc,1)  = motion components (string vector)
    %         sps.mon.dtm(mon.na,1) = time-steps        (real vector)
    %         sps.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
    %         sps.mon.ntm(mon.na,1) = time step number  (real vector)
    %         sps.syn(mon.na,1)     = structure vectors (cell vector)
    %         sps.syn{i}.tha.x      = x-acceleration    (real vector)
    %         sps.syn{i}.tha.y      = y-acceleration    (real vector)
    %         sps.syn{i}.tha.z      = z-acceleration    (real vector)
    %         sps.syn{i}.pga.x(1)   = x-time-pga        (real vector)
    %         sps.syn{i}.pga.x(2)   = x-pga             (real vector)
    %         sps.syn{i}.pga.y(1)   = y-time-pga        (real vector)
    %         sps.syn{i}.pga.y(2)   = y-pga             (real vector)
    %         sps.syn{i}.pga.z(1)   = z-time-pga        (real vector)
    %         sps.syn{i}.pga.z(2)   = z-pga             (real vector)
    %         sps.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
    %         sps.syn{i}.pgv.x(2)   = x-pgv             (real vector)
    %         sps.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
    %         sps.syn{i}.pgv.y(2)   = y-pgv             (real vector)
    %         sps.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
    %         sps.syn{i}.pgv.z(2)   = z-pgv             (real vector)
    %         sps.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
    %         sps.syn{i}.pgd.x(2)   = x-pgd             (real vector)
    %         sps.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
    %         sps.syn{i}.pgd.y(2)   = y-pgd             (real vector)
    %         sps.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
    %         sps.syn{i}.pgd.z(2)   = z-pgd             (real vector)     
    % OUTPUT: sps    = structure of sp synthetics
    %         sps.mtd     = metadata structure
    %         sps.mtd.mw  = moment magnitude (real)
    %         sps.mtd.scc = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)
    %         sps.mtd.sst = st. dev. of GMPE (0=median value,1=84th percentile)
    %         sps.mtd.scl = scale factor [1=cm/s/s]
    %         sps.mtd.ivd = flag for output in velocity and displacement
    %         sps.mon = monitor structure
    %         sps.mon.pt = path to monitor files    (string)
    %         sps.mon.fn = monitor metadata filename(string)
    %         sps.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
    %         sps.mon.id = monitor identity         (integer)
    %         sps.mon.dep = epicentral distance     (real vector)
    %         sps.mon.stn = monitor names           (string vector)
    %         sps.mon.na = number of monitors       (integer)
    %         sps.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
    %         sps.mon.nr = number of records        (integer)
    %         sps.mon.cp = motion component         (integer: 1,2,3)
    %         sps.mon.nc = number of components     (integer)
    %         sps.mon.na  = number of accelerograms to be generated
    %         sps.mon.nc  = number of motion components (integer)
    %         sps.mon.cp(mon.nc,1)  = motion components (string vector)
    %         sps.mon.dtm(mon.na,1) = time-steps        (real vector)
    %         sps.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
    %         sps.mon.ntm(mon.na,1) = time step number  (real vector)
    %         sps.mon.nT            = natural period number (real)
    %         sps.mon.vTn(mon.nT,1) = natural period vector (real vector)
    %         sps.mon.zeta          = damping               (real)
    %         sps.mon.dfr(mon.na,1) = frequency step   (real vector)    
    %         sps.mon.nfr(mon.na,1) = frequency number (real vector)
    %         sps.mon.nNy(mon.na,1) = Nyquist index    (real vector)
    %         sps.mon.vfr(mon.na,1) = frequency vector (real cell vector)
    %         sps.syn(mtd.na,1)     = structure vectors (cell vector)
    %         sps.syn{i}.tha.x      = x-acceleration    (real vector)
    %         sps.syn{i}.tha.y      = y-acceleration    (real vector)
    %         sps.syn{i}.tha.z      = z-acceleration    (real vector)
    %         sps.syn{i}.pga.x(1)   = x-time-pga        (real vector)
    %         sps.syn{i}.pga.x(2)   = x-pga             (real vector)
    %         sps.syn{i}.pga.y(1)   = y-time-pga        (real vector)
    %         sps.syn{i}.pga.y(2)   = y-pga             (real vector)
    %         sps.syn{i}.pga.z(1)   = z-time-pga        (real vector)
    %         sps.syn{i}.pga.z(2)   = z-pga             (real vector)
    %         sps.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
    %         sps.syn{i}.pgv.x(2)   = x-pgv             (real vector)
    %         sps.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
    %         sps.syn{i}.pgv.y(2)   = y-pgv             (real vector)
    %         sps.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
    %         sps.syn{i}.pgv.z(2)   = z-pgv             (real vector)
    %         sps.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
    %         sps.syn{i}.pgd.x(2)   = x-pgd             (real vector)
    %         sps.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
    %         sps.syn{i}.pgd.y(2)   = y-pgd             (real vector)
    %         sps.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
    %         sps.syn{i}.pgd.z(2)   = z-pgd             (real vector)     
    %         sps.syn{i}.rsd.x = x-displacement spectrum
    %         sps.syn{i}.rsd.y = y-displacement spectrum
    %         sps.syn{i}.rsd.z = z-displacement spectrum
    %         sps.syn{i}.psa.x = x-pseudo-acceleration spectrum
    %         sps.syn{i}.psa.y = y-pseudo-acceleration spectrum
    %         sps.syn{i}.psa.z = z-pseudo-acceleration spectrum
    %         sps.syn{i}.fsa.x = x-fourier spectrum
    %         sps.syn{i}.fsa.y = y-fourier spectrum
    %         sps.syn{i}.fsa.z = z-fourier spectrum
    % N.B. Need for newmark_sd.m
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    sps = varargin{1};
    %----------------------------------------------------------------------
    % response spectra parameters 
    %----------------------------------------------------------------------
    % minimum natural period
    Tn_min   = 0;
    % maximum natural period
    Tn_max   = 5;
    % natural period step
    dTn      = 0.05;
    sps.mon.vTn  = (Tn_min:dTn:Tn_max)';
    sps.mon.nT   = numel(sps.mon.vTn);
    sps.mon.zeta = 0.05;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % fourier spectra parameters 
    %----------------------------------------------------------------------
    for i_ = 1:sps.mon.na
        % Nyquist frequency
        fr_max      = 1/2/sps.mon.dtm(i_);
        % fft points
        sps.mon.nfr(i_) = 2^nextpow2(sps.mon.ntm(i_))+1;
        % frequency period step
        sps.mon.dfr(i_) = 1/sps.mon.dtm(i_)/(sps.mon.nfr(i_)-1);
        % Nyquist frequency index
        sps.mon.nNy(i_) = floor(fr_max/sps.mon.dfr(i_))+1;
        % frequency vector
        sps.mon.vfr(i_) = {sps.mon.dfr(i_)*(0:sps.mon.nfr(i_)-1)'};
    end
    %======================================================================
    %======================================================================
    % SD and PSA SPECTRA
    %======================================================================
    for i_ = 1:sps.mon.na
        for j_ = 1:sps.mon.nc
            eval(sprintf(['[sps.syn{i_}.rsd.%s,~,~,sps.syn{i_}.psa.%s,~] ='...
                'newmark_sd(sps.syn{i_}.tha.%s,sps.mon.dtm(i_),sps.mon.vTn,sps.mon.zeta);'],...
                sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_}));
        end
    end
    %======================================================================
    % FOURIER SPECTRUM
    %======================================================================
    for i_ = 1:sps.mon.na
        for j_ = 1:sps.mon.nc
            eval(sprintf(['[~,sps.syn{i_}.fsa.%s,~,~,~,~] ='...
                'super_fft(sps.mon.dtm(i_),sps.syn{i_}.tha.%s,0);'],...
                sps.mon.cp{j_},sps.mon.cp{j_}));
        end
    end
    %======================================================================
    varargout{1} = sps;
    return
end
