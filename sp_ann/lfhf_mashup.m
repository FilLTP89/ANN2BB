%% *Compute response spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _compute_spectra_: function to compute response spectra for a set of records
%% INPUT:
%%
% * _nss (Numeric simulation structure)_
% * nss.mon    = structure of monitor data
% * nss.mon.pt = path to monitor files        (string)
% * nss.mon.tp = type of monitor              (string: 'S'(speed),'H'(hisada))
% * nss.mon.id = monitor identity             (integer)
% * nss.mon.na = number of monitors           (integer)
% * nss.mon.rc = monitor record               (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * nss.mon.nr = number of records            (integer)
% * nss.mon.cp = motion component             (integer: 1,2,3)
% * nss.mon.nc = number of components         (integer)
% * nss.mon.dtm(mon.na,1) = time-steps        (real vector)
% * nss.mon.vtm(mon.na,1) = time-vectors      (real cell vector)
% * nss.mon.ntm(mon.na,1) = time step number  (real vector)
% * nss.mon  = response-spectra structure
% * nss.mon.nT            = natural period number (real)
% * nss.mon.vTn(mon.nT,1) = natural period vector (real vector)
% * nss.mon.zeta          = damping               (real)
% * nss.mon.dfr(mon.na,1) = frequency-step   (real vector)
% * nss.mon.nfr(mon.na,1) = frequency number (real vector)
% * nss.mon.nNy(mon.na,1) = Nyquist index    (real vector)
% * nss.mon.vfr(mon.na,1) = frequency vector (real cell vector)
%%
% * _nss.syn(mon.na,1)     = structure vector (cell vector)_
% * nss.syn{i}.tha.x      = x-acceleration    (real vector)
% * nss.syn{i}.tha.y      = y-acceleration    (real vector)
% * nss.syn{i}.tha.z      = z-acceleration    (real vector)
% * nss.syn{i}.thv.x      = x-velocity        (real vector)
% * nss.syn{i}.thv.y      = y-velocity        (real vector)
% * nss.syn{i}.thv.z      = z-velocity        (real vector)
% * nss.syn{i}.thd.x      = x-displacement    (real vector)
% * nss.syn{i}.thd.y      = y-displacement    (real vector)
% * nss.syn{i}.thd.z      = z-displacement    (real vector)
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
% * nss.syn{i}.rsd.x      = x-displacement response spectrum (real vectors)
% * nss.syn{i}.rsd.y      = y-displacement response spectrum (real vectors)
% * nss.syn{i}.rsd.z      = z-displacement response spectrum (real vectors)
% * nss.syn{i}.psa.x      = x-pseudo-acceleration spectrum
% * nss.syn{i}.psa.y      = y-pseudo-acceleration spectrum
% * nss.syn{i}.psa.z      = z-pseudo-acceleration spectrum
% * nss.syn{i}.fsa.x      = x-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsa.y      = y-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsa.z      = z-acceleration fourier spectrum (real vector)
% * nss.syn{i}.fsv.x      = x-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsv.y      = y-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsv.z      = z-velocity fourier spectrum (real vector)
% * nss.syn{i}.fsd.x      = x-displacement fourier spectrum (real vector)
% * nss.syn{i}.fsd.y      = y-displacement fourier spectrum (real vector)
% * nss.syn{i}.fsd.z      = z-displacement fourier spectrum (real vector)
%%
% * _sps  (Sabetta&Pugliese structure of synthetics)_
%%
% * _sps.mtd     = metadata structure_
% * sps.mtd.mw  = moment magnitude (real)
% * sps.mtd.dtm = time step        (real)
% * sps.mtd.scc = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)
% * sps.mtd.sst = st. dev. of GMPE (0=median value,1=84th percentile)
% * sps.mtd.scl = scale factor [1=cm/s/s]
% * sps.mtd.ivd = flag for output in velocity and displacement
%%
% * _sps.mon = monitor structure_
% * sps.mon.pt = path to monitor files    (string)
% * sps.mon.fn = monitor metadata filename(string)
% * sps.mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
% * sps.mon.id = monitor identity         (integer)
% * sps.mon.dep = epicentral distance     (real vector)
% * sps.mon.stn = monitor names           (string vector)
% * sps.mon.na = number of monitors       (integer)
% * sps.mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
% * sps.mon.nr = number of records        (integer)
% * sps.mon.cp = motion component         (integer: 1,2,3)
% * sps.mon.nc = number of components     (integer)
% * sps.mon.na  = number of accelerograms to be generated
% * sps.mon.nc  = number of motion components (integer)
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
%%
% * _sps.syn(mtd.na,1)     = structure vectors (cell vector)_
% * sps.syn{i}.tha.x      = x-acceleration    (real vector)
% * sps.syn{i}.tha.y      = y-acceleration    (real vector)
% * sps.syn{i}.tha.z      = z-acceleration    (real vector)
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
% * sps.syn{i}.rsd.x = x-displacement spectrum
% * sps.syn{i}.rsd.y = y-displacement spectrum
% * sps.syn{i}.rsd.z = z-displacement spectrum
% * sps.syn{i}.psa.x = x-pseudo-acceleration spectrum
% * sps.syn{i}.psa.y = y-pseudo-acceleration spectrum
% * sps.syn{i}.psa.z = z-pseudo-acceleration spectrum
% * sps.syn{i}.fsa.x = x-fourier spectrum
% * sps.syn{i}.fsa.y = y-fourier spectrum
% * sps.syn{i}.fsa.z = z-fourier spectrum
%% OUTPUT:
% * _hbs (hybrid synthetics structure)_
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
function [varargout] = lfhf_mashup(varargin)
    %% SET-UP
    nss = varargin{1};
    sps = varargin{2};
    %% HYBRIDIZATION
    %%
    % _resampling and padding records_
    [slf,shf,hbs] = lfhf_padding(nss,sps);
    %%
    % _arias intensity_
    I1 = 0.05;
    [slf,shf]= lfhf_arias(I1,slf,shf);
    %%
    % _shift of hf signal according to T05 of arias intensity_
    shf = lfhf_shift(slf,shf);
    %%
    % _hybrid signals_
    hbs = lfhf_hybridator(slf,shf,hbs);
    %% OUTPUT
    varargout{1} = hbs;
    return
end
%% *Resampling records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_padding_: function that resambles lf and hf records to have
% the same dt and number of points.
%% INPUT:
% * _nss (numerical simulation structure)_
% * _sps (Sabetta/Pugliese synthetics structure)_
%% OUTPUT:
% * _dtm (common time step)_
% * _ntm (common number of time-points)_
% * _slf (lf record's structure)_
% * _shf (hf record's structure)_
% * _hbs (hybrid record's structure)_
function [varargout] = lfhf_padding(varargin)
    %% SET-UP
    nss = varargin{1};
    sps = varargin{2};
    %%
    % _copy structures in LF/HF_
    slf.mon = nss.mon;
    shf.mon = sps.mon;
    hbs.mon = nss.mon;
    slf.syn = nss.syn;
    shf.syn = sps.syn;
    hbs.syn = cell(nss.mon.na,1);
    
    % frequency band
    
    if strcmpi(hbs.mon.tp,'h')
        hbs.mon.fa = 1; % in Hz
        hbs.mon.fb = 1.5; % in Hz
    else
        hbs.mon.fa = 1; % in Hz
        hbs.mon.fb = 2; % in Hz
    end
    %% RESAMPLING
    %  _check time-step_
    [status,dtm,idx1,idx2] = check_dt(slf.mon.dtm,shf.mon.dtm);
    hbs.mon.dtm(:)   = dtm;
    if status
        % resampling
        for i_ = 1:numel(idx1)
            j_ = idx1(i_);
            k_ = idx2(i_);
            if logical(k_-1)
                fac = shf.mon.dtm(j_)/dtm(j_);
                [new_dtm,shf.syn{j_}.tha,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(shf.mon.dtm(j_),v,...
                    fac),shf.syn{j_}.tha,'UniformOutput',0);
                
                shf.mon.ntm(j_) = new_ntm.x;
                shf.mon.dtm(j_) = new_dtm.x;
                shf.mon.vtm{j_} = new_vtm.x;
            else
                fac = slf.mon.dtm(j_)/dtm(j_);
                [new_dtm,slf.syn{j_}.tha,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(slf.mon.dtm(j_),v,...
                    fac),slf.syn{j_}.tha,'UniformOutput',0);
                
                slf.mon.dtm(j_) = new_dtm.x;
                slf.mon.ntm(j_) = new_ntm.x;
                slf.mon.vtm{j_} = new_vtm.x;
            end
        end
    end
    %% PADDING RECORDS
    for i_ = 1:slf.mon.na
        [ntm,idx] = max([slf.mon.ntm(i_),shf.mon.ntm(i_)]);
        if idx==2     % padding low frequency
            for j_ = 1:slf.mon.nc
                eval(sprintf(['slf.syn{i_}.tha.%s(slf.mon.ntm(i_)+1:ntm) =',...
                    'zeros(ntm-slf.mon.ntm(i_),1);'],slf.mon.cp{j_},slf.mon.cp{j_}));
            end
        elseif idx==1 % padding high frequency
            for j_ = 1:shf.mon.nc
                eval(sprintf(['shf.syn{i_}.tha.%s(shf.mon.ntm(i_)+1:ntm) = ',...
                    'zeros(ntm-shf.mon.ntm(i_),1);'],shf.mon.cp{j_},shf.mon.cp{j_}));
            end
        end
        slf.mon.ntm(i_) = ntm;
        shf.mon.ntm(i_) = ntm;
        hbs.mon.ntm(i_) = ntm;
        slf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:ntm-1)'};
        shf.mon.vtm(i_) = {shf.mon.dtm(i_)*(0:ntm-1)'};
        hbs.mon.vtm(i_) = {hbs.mon.dtm(i_)*(0:ntm-1)'};
    end
    %% OUTPUT
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hbs;
    return
end

%% *Check time-steps and resample records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _check_dt_: function that check if time-step are similar and resamples
% the lf hf records accordingly.
%% INPUT:
% * _dt1 (time-step vector)_
% * _dt2 (time-step vector)_
%% OUTPUT:
% * _dtm (common time step)_
% * _ntm (common number of time-points)_
function [varargout] = check_dt(varargin)
    %% SET-UP
    dt1 = varargin{1}(:);
    dt2 = varargin{2}(:);
    %% CHECK AND CORRECT DT
    [dtm,~ ] = min([dt1 dt2],[],2);
    [~,idx2] = max([dt1 dt2],[],2);
    idx1 = find(abs(dt1-dt2)>1e-10);
    
    if any(idx1)
        status  = 1;
    end
    %% OUTPUT
    varargout{1} = status;
    varargout{2} = dtm;
    varargout{3} = idx1;
    varargout{4} = idx2;
    return
end

%% *Compute Arias Intensity*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_arias_: compute arias intensity and T5% for lf and hf records
%% INPUT:
% * _I1 (percentage of Arias Intensity)_
% * _slf(low frequency records)_
% * _shf(high frequency records)_
%% OUTPUT:
% * _slf(low frequency records)_
% * _shf(high frequency records)_
function [varargout] = lfhf_arias(varargin)
    %% SET-UP
    I1  = varargin{1};
    slf = varargin{2};
    shf = varargin{3};
    %% ARIAS INTENSITY
    for i_ = 1:slf.mon.na
        for j_ = 1:slf.mon.nc
            % low frequency-arias intensity
            eval(sprintf(['[slf.syn{i_}.AT5.%s,slf.syn{i_}.AI5.%s,',...
                'slf.syn{i_}.Ain.%s] = arias_intensity(',...
                'slf.syn{i_}.tha.%s,slf.mon.dtm(i_),I1);'],...
                slf.mon.cp{j_},slf.mon.cp{j_},slf.mon.cp{j_},slf.mon.cp{j_}));
            % high frequency-arias intensity
            eval(sprintf(['[shf.syn{i_}.AT5.%s,shf.syn{i_}.AI5.%s,',...
                'shf.syn{i_}.Ain.%s] = arias_intensity(',...
                'shf.syn{i_}.tha.%s,shf.mon.dtm(i_),I1);'],...
                shf.mon.cp{j_},shf.mon.cp{j_},shf.mon.cp{j_},shf.mon.cp{j_}));
        end
    end
    %% OUTPUT
    varargout{1} = slf;
    varargout{2} = shf;
    return
end

%% *Shift time records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_shift_: shift hf records according to the difference in T5%
% between lf and hf records
%% INPUT:
% * _slf(low frequency records)_
% * _shf(high frequency records)_
%% OUTPUT:
% * _shf(high frequency records)_
function [varargout] = lfhf_shift(varargin)
    %% SET-UP
    slf = varargin{1};
    shf = varargin{2};
    %% TIME SHIFTING
    
    for i_ = 1:slf.mon.na
        for j_ = 1:slf.mon.nc
            eval(sprintf('hf = zeros(size(shf.syn{i_}.tha.%s));',slf.mon.cp{j_}));
            eval(sprintf(['is = shf.syn{i_}.AI5.%s-',...
                'slf.syn{i_}.AI5.%s;'],slf.mon.cp{j_},slf.mon.cp{j_}));
            ss = (is>=0)+1;
            is = abs(is)+1;
            if ss == 2 % positive is
                idx00 = is:slf.mon.ntm;
                idx01 = idx00-is+1;
                idx11 = idx01(end)+1:idx01(end)+is-1;
                eval(sprintf('hf(idx01) = shf.syn{i_}.tha.%s(idx00);',...
                    slf.mon.cp{j_}));
                hf(idx11) = 0.;
            else % negative is
                idx0  = 1:is-1;
                idx00 = is:slf.mon.ntm;
                idx01 = idx00-is+1;
                hf(idx0) = 0.;
                eval(sprintf('hf(idx00) = shf.syn{i_}.tha.%s(idx01);',...
                    slf.mon.cp{j_}));
            end
            eval(sprintf('shf.syn{i_}.tha.%s = hf;',shf.mon.cp{j_}));
        end
    end
    %% OUTPUT
    varargout{1} = shf;
    return
end

    %% *Fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_hybridator_: function to hybridize lf and hf records
%% INPUT: 
% * _slf(low frequency records)_
% * _shf(high frequency records)_
% * _hbs(hybrid records)_
%% OUTPUT: 
% * _slf(low frequency records)_
% * _shf(high frequency records)_
% * _hbs(hybrid records)_
function [varargout] = lfhf_hybridator(varargin)
    
    %% SET-UP
    slf = varargin{1};
    shf = varargin{2};
    hbs = varargin{3};
    %% BROAD-BAND SIGNALS
    for i_ = 1:slf.mon.na
        nfr = 2^nextpow2(hbs.mon.ntm(i_))+1;
        dfr = 1/hbs.mon.dtm(i_)/(nfr-1);
        vfr = dfr*(0:nfr-1)';
        fr_max = 0.5/hbs.mon.dtm(i_);
        nfa = round(hbs.mon.fa/dfr);
        nfb = round(hbs.mon.fb/dfr);
        fac = pi./(dfr*(nfb-nfa-1));
        hbs.mon.nfr(i_) = nfr;
        hbs.mon.dfr(i_) = dfr;
        hbs.mon.vfr(i_) = {vfr};
        hbs.mon.nNy(i_) = floor(fr_max/dfr)+1;
        
        for j_ = 1:slf.mon.nc            
            %%
            % _Fourier spectra_
            %             eval(sprintf('slf.syn{i_}.tha.%s(hbs.mon.ntm(i_)+1:nfr) = 0.;',...
            %                 slf.mon.cp{j_}));
            %             eval(sprintf('shf.syn{i_}.tha.%s(hbs.mon.ntm(i_)+1:nfr) = 0.;',...
            %                 slf.mon.cp{j_}));
            eval(sprintf(['slf.syn{i_}.fsa.%s = ',...
                'slf.mon.dtm(i_)*fft(slf.syn{i_}.tha.%s,nfr);'],...
                slf.mon.cp{j_},slf.mon.cp{j_}));
            eval(sprintf(['shf.syn{i_}.fsa.%s = ',...
                'shf.mon.dtm(i_)*fft(shf.syn{i_}.tha.%s,nfr);'],...
                slf.mon.cp{j_},slf.mon.cp{j_}));
            %% 
            % _wighting lf and hf_
            
            for k_ = 1:hbs.mon.nfr
                % weight LF part
                WLF = 1*(k_<=nfa)+...
                    (0.5+0.5*cos(fac*hbs.mon.dfr(i_)*(k_-nfa-1)))*(k_>nfa && k_ <= nfb);
                % weight HF part
                WHF = 1-WLF;
                eval(sprintf(['hbs.syn{i_}.fsa.%s(k_) = ',...
                    'WLF.*slf.syn{i_}.fsa.%s(k_)+',...
                    'WHF.*shf.syn{i_}.fsa.%s(k_);'],...
                    slf.mon.cp{j_},slf.mon.cp{j_},slf.mon.cp{j_}));
                
            end
            %%
            % _Inverse Fourier Transform_            
            eval(sprintf(['hbs.syn{i_}.tha.%s = ',...
                'real(ifft(hbs.syn{i_}.fsa.%s,hbs.mon.ntm(i_),''symmetric''))./',...
                'hbs.mon.dtm(i_);'],...
                slf.mon.cp{j_},slf.mon.cp{j_}));
            eval(sprintf(['[hbs.syn{i_}.pga.%s(1),hbs.syn{i_}.pga.%s(2),'...
                'hbs.syn{i_}.pgv.%s(1),hbs.syn{i_}.pgv.%s(2),',...
                'hbs.syn{i_}.pgd.%s(1),hbs.syn{i_}.pgd.%s(2)] = ',...
                'PGAVD_eval(hbs.mon.dtm(i_),hbs.syn{i_}.tha.%s);'],...
                hbs.mon.cp{j_},hbs.mon.cp{j_},hbs.mon.cp{j_},hbs.mon.cp{j_},...
                hbs.mon.cp{j_},hbs.mon.cp{j_},hbs.mon.cp{j_}));
            eval(sprintf(['[hbs.syn{i_}.rsd.%s,~,~,hbs.syn{i_}.psa.%s,~] ='...
                'newmark_sd(hbs.syn{i_}.tha.%s,hbs.mon.dtm(i_),hbs.mon.vTn,',...
                'hbs.mon.zeta);'],...
                hbs.mon.cp{j_},hbs.mon.cp{j_},hbs.mon.cp{j_}));
            
            %%
            % _Filtering_
            %
            %             eval(sprintf(['[hbs.syn{i_}.tha.%s,hbs.syn{i_}.thv.%s,',...
            %                 'hbs.syn{i_}.thd.%s] = band_pass_filter(hbs.mon.dtm(i_),',...
            %                 'hbs.syn{i_}.tha.%s);'],hbs.mon.cp{j_},hbs.mon.cp{j_},...
            %                 hbs.mon.cp{j_},hbs.mon.cp{j_}));
            
            % _Final Fourier Spectra_
            %             eval(sprintf(['[~,hbs.syn{i_}.fsa.%s,~,~,~,~] = ',...
            %                 'super_fft(hbs.mon.dtm(i_),hbs.syn{i_}.tha.%s,0,hbs.mon.ntm);'],...
            %                 hbs.mon.cp{j_},hbs.mon.cp{j_}));
            
        end
    end
    %% OUTPUT
    varargout{1} = hbs;
    return
end
