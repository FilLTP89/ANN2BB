function [varargout] = lfhf_mashup(varargin)
    %===============
    % Compute response spectra
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % compute_spectra: function to compute response spectra for a set of records
    % INPUT:  nss (numerical simulation structure)
    %         sps (Sabetta/Pugliese synthetics structure)
    % OUTPUT: hbs (hybrid synthetics structure)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    nss = varargin{1};
    sps = varargin{2};
    %----------------------------------------------------------------------
    % parameters
    %----------------------------------------------------------------------
    icorr = 0;
    tpad = 15;
    %======================================================================
    %======================================================================
    % HYBRIDIZATION
    %======================================================================
    %----------------------------------------------------------------------
    % adjust vectors
    %----------------------------------------------------------------------
    [slf,shf,hbs] = lfhf_padding(nss,sps);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % arias intensity
    %----------------------------------------------------------------------
    I1 = 0.05;
    [slf,shf]= lfhf_arias(I1,slf,shf);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % shift of hf signal according to T05 of arias intensity
    %----------------------------------------------------------------------
    shf = lfhf_shift(slf,shf);
    %----------------------------------------------------------------------
    % hybrid signals
    %----------------------------------------------------------------------
    hbs = lfhf_hybridator(slf,shf,hbs);
    %======================================================================
    varargout{1} = hbs;
    return
end
%
function [varargout] = lfhf_padding(varargin)
    %===============
    % Resampling records
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % lfhf_padding: function that resambles lf and hf records to have
    % the same dt and number of points.
    % INPUT:  nss (numerical simulation structure)
    %         sps (Sabetta/Pugliese synthetics structure)
    % OUTPUT: dtm (common time step)
    %         ntm (common number of time-points)
    %         slf (lf record's structure)
    %         shf (hf record's structure)
    %         hbs (hybrid record's structure)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    nss = varargin{1};
    sps = varargin{2};
    
    slf.mon = nss.mon;
    shf.mon = sps.mon;
    hbs.mon = nss.mon;
    slf.syn = nss.syn;
    shf.syn = sps.syn;
    hbs.syn = cell(nss.mon.na,1);
    %----------------------------------------------------------------------
    % frequency band
    %----------------------------------------------------------------------
    hbs.mon.fa = 1.5; % in Hz
    hbs.mon.fb = 2.0; % in Hz
    %======================================================================
    %======================================================================
    % CHECK DT
    %======================================================================
    [status,dtm] = check_dt(slf.mon.dtm,shf.mon.dtm);
    
    if status
        hbs.mon.dtm(:) = dtm;
    else
        error('different time steps');
    end
    %======================================================================
    %======================================================================
    % PADDING RECORDS
    %======================================================================
    for i_ = 1:nss.mon.na
        
        [ntm,idx] = max([slf.mon.ntm(i_),shf.mon.ntm(i_)]);
        
        if idx==2     % padding low frequency
            for j_ = 1:slf.mon.nc
                eval(sprintf(['slf.syn{i_}.tha.%s(slf.mon.ntm(i_)+1:ntm) =',...
                    'zeros(ntm-slf.mon.ntm(i_),1);'],slf.mon.cp{j_},slf.mon.cp{j_}));
            end
        elseif idx==1 % padding high frequency
            for j_ = 1:slf.mon.nc
                eval(sprintf(['shf.syn{i_}.tha.%s(shf.mon.ntm(i_)+1:ntm) = ',...
                    'zeros(ntm-shf.mon.ntm(i_),1);'],...
                    slf.mon.cp{j_}));
            end
        end
        slf.mon.ntm(i_) = ntm;
        shf.mon.ntm(i_) = ntm;
        hbs.mon.ntm(i_) = ntm;
        slf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:ntm-1)'};
        shf.mon.vtm(i_) = {shf.mon.dtm(i_)*(0:ntm-1)'};
        hbs.mon.vtm(i_) = {hbs.mon.dtm(i_)*(0:ntm-1)'};
    end
    %======================================================================
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hbs;
    return
end
%
function [varargout] = check_dt(varargin)
    %===============
    % check and uniform time-steps
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % check_dt: function that check if time-step are similar and resamples
    % the lf hf records accordingly.
    % INPUT:  dt1 (time-step vector)
    %         dt2 (time-step vector)
    % OUTPUT: dtm (common time step)
    %         ntm (common number of time-points)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    dt1 = varargin{1};
    dt2 = varargin{2};
    status = 1;
    flag   = 0;
    dtm    = dt1;
    %======================================================================
    %======================================================================
    % CHECK AND CORRECT DT
    %======================================================================
    for i_ = 1:numel(dt1)
        for j_ = 1:numel(dt2)
            if abs(dt1(i_)-dt2(j_))>1e-10
                status=0;
                flag=1;
                keyboard
                break;
            end
        end
        if flag
            break;
        end
    end
    %======================================================================
    varargout{1} = status;
    varargout{2} = dtm;
    return
end
%==========================================================================
%==========================================================================
function [varargout] = lfhf_arias(varargin)
    %===============
    % Compute Arias Intensity
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % lfhf_arias: compute arias intensity and T5% for lf and hf records
    % INPUT:  I1 (percentage of Arias Intensity)
    %         slf(low frequency records)
    %         shf(high frequency records)
    % OUTPUT: slf(low frequency records)
    %         shf(high frequency records)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    I1  = varargin{1};
    slf = varargin{2};
    shf = varargin{3};
    %======================================================================
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
    %======================================================================
    varargout{1} = slf;
    varargout{2} = shf;
    return
end
%==========================================================================
%==========================================================================
function [varargout] = lfhf_shift(varargin)
    %===============
    % Shift time records
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % lfhf_shift: shift hf records according to the difference in T5%
    % between lf and hf records
    % INPUT:  slf(low frequency records)
    %         shf(high frequency records)
    % OUTPUT: shf(high frequency records)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    slf = varargin{1};
    shf = varargin{2};
    %======================================================================
    %======================================================================
    % TIME SHIFTING
    %======================================================================
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
    %======================================================================
    varargout{1} = shf;
    return
end
%==========================================================================
%==========================================================================
function [varargout] = lfhf_hybridator(varargin)
    %===============
    % Fourier spectra
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % lfhf_hybridator: function to hybridize lf and hf records
    % INPUT:  slf(low frequency records)
    %         shf(high frequency records)
    %         hbs(hybrid records)
    % OUTPUT: slf(low frequency records)
    %         shf(high frequency records)
    %         hbs(hybrid records)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    slf = varargin{1};
    shf = varargin{2};
    hbs = varargin{3};
    %======================================================================
    %======================================================================
    % BROAD-BAND SIGNALS
    %======================================================================
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
            %--------------------------------------------------------------
            % Fourier spectra
            %--------------------------------------------------------------
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
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            % wighting lf and hf
            %--------------------------------------------------------------
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
                %----------------------------------------------------------
            end
            %--------------------------------------------------------------
            % Inverse Fourier Transform
            %--------------------------------------------------------------
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
            %--------------------------------------------------------------
%             %--------------------------------------------------------------
%             % Filtering
%             %--------------------------------------------------------------
%             eval(sprintf(['[hbs.syn{i_}.tha.%s,hbs.syn{i_}.thv.%s,',...
%                 'hbs.syn{i_}.thd.%s] = band_pass_filter(hbs.mon.dtm(i_),',...
%                 'hbs.syn{i_}.tha.%s);'],hbs.mon.cp{j_},hbs.mon.cp{j_},...
%                 hbs.mon.cp{j_},hbs.mon.cp{j_}));
%             %--------------------------------------------------------------
%             % FINAL FOURIER TRANSFORM
%             %--------------------------------------------------------------
%             eval(sprintf(['[~,hbs.syn{i_}.fsa.%s,~,~,~,~] = ',...
%                 'super_fft(hbs.mon.dtm(i_),hbs.syn{i_}.tha.%s,0,hbs.mon.ntm);'],...
%                 hbs.mon.cp{j_},hbs.mon.cp{j_}));
%             %--------------------------------------------------------------
        end
    end
    varargout{1} = hbs;
    return
end