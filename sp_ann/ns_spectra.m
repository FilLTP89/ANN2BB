function [varargout] = ns_spectra(varargin)
    %===============
    % Compute response/fourier spectra
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % ns_spectra: function to compute response/fourier spectra for a set
    % of synthetics records outcomes of numerical simulations
    % INPUT:  nss    = structure of numerical simulations
    %         nss.mon    = structure of monitor data
    %         nss.mon.pt = path to monitor files        (string)
    %         nss.mon.tp = type of monitor              (string: 'S'(speed),'H'(hisada))
    %         nss.mon.id = monitor identity             (integer)
    %         nss.mon.na = number of monitors           (integer)
    %         nss.mon.rc = monitor record               (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
    %         nss.mon.nr = number of records            (integer)
    %         nss.mon.cp = motion component             (integer: 1,2,3)
    %         nss.mon.nc = number of components         (integer)
    %         nss.mon.dtm(mon.na,1) = time-steps        (real vector)
    %         nss.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
    %         nss.mon.ntm(mon.na,1) = time step number  (real vector)
    %         nss.syn(mon.na,1)     = structure vector  (cell vector)
    %         nss.syn{i}.tha.x      = x-acceleration    (real vector)
    %         nss.syn{i}.tha.y      = y-acceleration    (real vector)
    %         nss.syn{i}.tha.z      = z-acceleration    (real vector)
    %         nss.syn{i}.thv.x      = x-velocity        (real vector)
    %         nss.syn{i}.thv.y      = y-velocity        (real vector)
    %         nss.syn{i}.thv.z      = z-velocity        (real vector)
    %         nss.syn{i}.thd.x      = x-displacement    (real vector)
    %         nss.syn{i}.thd.y      = y-displacement    (real vector)
    %         nss.syn{i}.thd.z      = z-displacement    (real vector)
    %         cp (motion component)
    % OUTPUT: nss (vector simulation structure)
    %         nss.mon    = structure of monitor data
    %         nss.mon.pt = path to monitor files        (string)
    %         nss.mon.tp = type of monitor              (string: 'S'(speed),'H'(hisada))
    %         nss.mon.id = monitor identity             (integer)
    %         nss.mon.na = number of monitors           (integer)
    %         nss.mon.rc = monitor record               (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
    %         nss.mon.nr = number of records            (integer)
    %         nss.mon.cp = motion component             (integer: 1,2,3)
    %         nss.mon.nc = number of components         (integer)
    %         nss.mon.dtm(mon.na,1) = time-steps        (real vector)
    %         nss.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
    %         nss.mon.ntm(mon.na,1) = time step number  (real vector)
    %         nss.mon  = response-spectra structure
    %         nss.mon.nT            = natural period number (real)
    %         nss.mon.vTn(mon.nT,1) = natural period vector (real vector)
    %         nss.mon.zeta          = damping               (real)
    %         nss.mon.dfr(mon.na,1) = frequency-step   (real vector)
    %         nss.mon.nfr(mon.na,1) = frequency number (real vector)
    %         nss.mon.nNy(mon.na,1) = Nyquist index    (real vector)
    %         nss.mon.vfr(mon.na,1) = frequency vector (real cell vector)
    %         nss.syn(mon.na,1)     = structure vector (cell vector)
    %         nss.syn{i}.tha.x      = x-acceleration    (real vector)
    %         nss.syn{i}.tha.y      = y-acceleration    (real vector)
    %         nss.syn{i}.tha.z      = z-acceleration    (real vector)
    %         nss.syn{i}.thv.x      = x-velocity        (real vector)
    %         nss.syn{i}.thv.y      = y-velocity        (real vector)
    %         nss.syn{i}.thv.z      = z-velocity        (real vector)
    %         nss.syn{i}.thd.x      = x-displacement    (real vector)
    %         nss.syn{i}.thd.y      = y-displacement    (real vector)
    %         nss.syn{i}.thd.z      = z-displacement    (real vector)
    %         nss.syn{i}.pga.x(1)   = x-time-pga        (real vector)
    %         nss.syn{i}.pga.x(2)   = x-pga             (real vector)
    %         nss.syn{i}.pga.y(1)   = y-time-pga        (real vector)
    %         nss.syn{i}.pga.y(2)   = y-pga             (real vector)
    %         nss.syn{i}.pga.z(1)   = z-time-pga        (real vector)
    %         nss.syn{i}.pga.z(2)   = z-pga             (real vector)
    %         nss.syn{i}.pgv.x(1)   = x-time-pgv        (real vector)
    %         nss.syn{i}.pgv.x(2)   = x-pgv             (real vector)
    %         nss.syn{i}.pgv.y(1)   = y-time-pgv        (real vector)
    %         nss.syn{i}.pgv.y(2)   = y-pgv             (real vector)
    %         nss.syn{i}.pgv.z(1)   = z-time-pgv        (real vector)
    %         nss.syn{i}.pgv.z(2)   = z-pgv             (real vector)
    %         nss.syn{i}.pgd.x(1)   = x-time-pgd        (real vector)
    %         nss.syn{i}.pgd.x(2)   = x-pgd             (real vector)
    %         nss.syn{i}.pgd.y(1)   = y-time-pgd        (real vector)
    %         nss.syn{i}.pgd.y(2)   = y-pgd             (real vector)
    %         nss.syn{i}.pgd.z(1)   = z-time-pgd        (real vector)
    %         nss.syn{i}.pgd.z(2)   = z-pgd             (real vector)        
    %         nss.syn{i}.rsd.x      = x-displacement response spectrum (real vectors)
    %         nss.syn{i}.rsd.y      = y-displacement response spectrum (real vectors)
    %         nss.syn{i}.rsd.z      = z-displacement response spectrum (real vectors)
    %         nss.syn{i}.psa.x      = x-pseudo-acceleration spectrum
    %         nss.syn{i}.psa.y      = y-pseudo-acceleration spectrum
    %         nss.syn{i}.psa.z      = z-pseudo-acceleration spectrum
    %         nss.syn{i}.fsa.x      = x-acceleration fourier spectrum (real vector)
    %         nss.syn{i}.fsa.y      = y-acceleration fourier spectrum (real vector)
    %         nss.syn{i}.fsa.z      = z-acceleration fourier spectrum (real vector)
    %         nss.syn{i}.fsv.x      = x-velocity fourier spectrum (real vector)
    %         nss.syn{i}.fsv.y      = y-velocity fourier spectrum (real vector)
    %         nss.syn{i}.fsv.z      = z-velocity fourier spectrum (real vector)
    %         nss.syn{i}.fsd.x      = x-displacement fourier spectrum (real vector)
    %         nss.syn{i}.fsd.y      = y-displacement fourier spectrum (real vector)
    %         nss.syn{i}.fsd.z      = z-displacement fourier spectrum (real vector)
    % N.B. Need for newmark_sd.m
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    nss = varargin{1};
    %----------------------------------------------------------------------
    % response spectra parameters
    %----------------------------------------------------------------------
    % minimum natural period
    Tn_min   = 0;
    % maximum natural period
    Tn_max   = 5;
    % natural period step
    dTn      = 0.05;
    nss.mon.vTn  = (Tn_min:dTn:Tn_max)';
    nss.mon.nT   = numel(nss.mon.vTn);
    nss.mon.zeta = 0.05;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % fourier spectra parameters
    %----------------------------------------------------------------------
    for i_ = 1:nss.mon.na
        % Nyquist frequency
        fr_max      = 1/2/nss.mon.dtm(i_);
        % fft points
        nss.mon.nfr(i_) = 2^nextpow2(nss.mon.ntm(i_))+1;
        % frequency period step
        nss.mon.dfr(i_) = 1/nss.mon.dtm(i_)/(nss.mon.nfr(i_)-1);
        % Nyquist frequency index
        nss.mon.nNy(i_) = floor(fr_max/nss.mon.dfr(i_))+1;
        % frequency vector
        nss.mon.vfr(i_) = {nss.mon.dfr(i_)*(0:nss.mon.nfr(i_)-1)'};
    end
    %======================================================================
    %======================================================================
    % SD and PSA SPECTRA
    %======================================================================
    for i_ = 1:nss.mon.na
        for j_ = 1:nss.mon.nc
            eval(sprintf(['[nss.syn{i_}.rsd.%s,~,~,nss.syn{i_}.psa.%s,~] ='...
                'newmark_sd(nss.syn{i_}.tha.%s,nss.mon.dtm(i_),nss.mon.vTn,',...
                'nss.mon.zeta);'],...
                nss.mon.cp{j_},nss.mon.cp{j_},nss.mon.cp{j_}));
        end
    end
    %======================================================================
    %======================================================================
    % FOURIER SPECTRUM
    %======================================================================
    for i_ = 1:nss.mon.na
        for k_ = 1:nss.mon.nc
            eval(sprintf(['[~,nss.syn{i_}.fsa.%s,~,~,~,~] ='...
                'super_fft(nss.mon.dtm(i_),nss.syn{i_}.tha.%s,0);'],...
                nss.mon.cp{k_},nss.mon.cp{k_}));
        end
    end
    %======================================================================
    varargout{1} = nss;
    return
end
