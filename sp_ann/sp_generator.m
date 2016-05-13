function [varargout] = sp_generator(varargin)
    %===============
    % Generate synthetics - Sabetta & Pugliese approach
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % sp_generator: function to produce non-stationary accelerograms
    % according to the approach by Sabetta & Pugliese, 1996.
    % INPUT:  mon    = monitor structure
    %         mon.pt = path to monitor files    (string)
    %         mon.fn = monitor metadata filename(string)
    %         mon.tp = type of monitor          (string: 'S'(speed),'H'(hisada))
    %         mon.id = monitor identity         (integer)
    %         mon.dep = epicentral distance     (real vector)
    %         mon.stn = monitor names           (string vector)
    %         mon.na = number of monitors       (integer)
    %         mon.rc = monitor record           (string: 'a'(acceleration),'v'(velocity),'d'(displacement))
    %         mon.nr = number of records        (integer)
    %         mon.cp = motion component         (integer: 1,2,3)
    %         mon.nc = number of components     (integer)
    %         str    = extra metadata file      (string)
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
    % N.B.: it has been verified that with ssc = 2 the method of sabetta
    % gives as output acceleration THs with an average response spectrum
    % compliant with the (median+isig/2) spectrum given by the GMPE
    % N.B. Need for sabetta.m
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    sps.mon    = varargin{1};
    sps.mon.fn = varargin{2};
    %======================================================================
    %======================================================================
    % PARSING EXTRA METADATA
    %======================================================================
    str = importdata(sps.mon.fn);
    % earthquake magnitude
    sps.mtd.mw  = str(1,1);
    % site condition factor
    sps.mtd.scc = str(2,1);
    % site condition std
    sps.mtd.sst = str(3,1);
    % scale factor
    sps.mtd.scl = str(4,1);
    % output flag
    sps.mtd.ivd = str(5,1);
    %======================================================================
    %======================================================================
    % SABETTA & PUGLIESE SYNTHETICS
    %======================================================================
    sps.syn = cell(sps.mon.na,1);
    for i_ = 1:sps.mon.na % number of components
        [vtm,tha]=sabetta(sps.mtd.mw,sps.mon.dep(i_),sps.mtd.scc,...
            sps.mtd.sst,sps.mon.dtm(i_),1,sps.mtd.scl);
        for j_ = 1:sps.mon.nc
            eval(sprintf(['[sps.syn{i_}.tha.%s,sps.syn{i_}.thv.%s,',...
                'sps.syn{i_}.thd.%s] = band_pass_filter(sps.mon.dtm(i_),',...
                'tha);'],sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_}));
            eval(sprintf(['[sps.syn{i_}.pga.%s(1),sps.syn{i_}.pga.%s(2),'...
                'sps.syn{i_}.pgv.%s(1),sps.syn{i_}.pgv.%s(2),',...
                'sps.syn{i_}.pgd.%s(1),sps.syn{i_}.pgd.%s(2)] = ',...
                'PGAVD_eval(sps.mon.dtm(i_),sps.syn{i_}.tha.%s,',...
                'sps.syn{i_}.thv.%s,sps.syn{i_}.thd.%s);'],...
                sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_},...
                sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_},sps.mon.cp{j_},...
                sps.mon.cp{j_}));
        end
        sps.mon.vtm(i_) = {vtm};
        sps.mon.ntm(i_) = numel(vtm);
    end
    %======================================================================
    varargout{1} = sps;
    return
end
