%% *Compute response spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _compute_spectra_: function to compute response spectra for a set of records
%% INPUT:
%%
% * _slf (structure of numerical simulations)_
% * _shf (structure of sp synthetics)_
%% OUTPUT:
% * _slf (structure of numerical simulations)_
% * _shf (structure of sp synthetics)_
% * _hyb (hybrid synthetics structure)_
function [varargout] = lfhf_mashup(varargin)
    %% SET-UP
    slf = varargin{1};
    shf = varargin{2};
    hyb.mon = slf.mon;
    hyb.syn = cell(slf.mon.na,1);
    %% HYBRIDIZATION
    [slf,shf,hyb] = lfhf_hybridator(slf,shf,hyb);
    %% OUTPUT
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hyb;
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
% * _hyb(hybrid records)_
%% OUTPUT:
% * _slf(low frequency records)_
% * _shf(high frequency records)_
% * _hyb(hybrid records)_
function [varargout] = lfhf_hybridator(varargin)
    %% SET-UP
    slf = varargin{1};
    shf = varargin{2};
    hyb = varargin{3};
    %% BROAD-BAND SIGNALS
    for i_ = 1:slf.mon.na
        %         nfr = 2^nextpow2(hyb.mon.ntm(i_));
        %         dfr = 1/hyb.mon.dtm(i_)/(nfr-1);
        %         vfr = dfr*(0:nfr-1)';
        %         fr_max = 0.5/hyb.mon.dtm(i_);
        nfa = round(hyb.mon.fa/hyb.mon.dfr(i_));
        nfb = round(hyb.mon.fb/hyb.mon.dfr(i_));
        fac = pi./(hyb.mon.dfr(i_)*(nfb-nfa-1));
        %         hyb.mon.nfr(i_) = nfr;
        %         hyb.mon.dfr(i_) = dfr;
        %         hyb.mon.vfr(i_) = {vfr};
        %         hyb.mon.nNy(i_) = floor(fr_max/dfr)+1;
        for j_ = 1:slf.mon.nc
            %% *FOURIER TRANSFORM*
            slf.syn{i_}.fsa.(slf.mon.cp{j_}) = ...
                fft(slf.syn{i_}.tha.(slf.mon.cp{j_}),hyb.mon.nfr(i_)).*...
                hyb.mon.dtm(i_);
            shf.syn{i_}.fsa.(shf.mon.cp{j_}) = ...
                fft(shf.syn{i_}.tha.(shf.mon.cp{j_}),hyb.mon.nfr(i_)).*...
                hyb.mon.dtm(i_);
            %% *FILTERING*
            WLF = zeros(hyb.mon.nfr(i_),1);
            WLF(1:nfa) = 1.0;
            WLF(nfa+1:nfb) = 0.5*(1+cos(fac*hyb.mon.dfr(i_)*(0:nfb-nfa-1)));
            WHF = 1-WLF;
            %
            slf.syn{i_}.fsa.(hyb.mon.cp{j_}) = ...
                WLF.*slf.syn{i_}.fsa.(hyb.mon.cp{j_});
            %
            shf.syn{i_}.fsa.(hyb.mon.cp{j_}) = ...
                WHF.*shf.syn{i_}.fsa.(hyb.mon.cp{j_});
            %
            hyb.syn{i_}.fsa.(hyb.mon.cp{j_}) = ...
                slf.syn{i_}.fsa.(hyb.mon.cp{j_}) +...
                shf.syn{i_}.fsa.(hyb.mon.cp{j_});
            %             % _conjugate_
            %             slf.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+1) = 0.0;
            %             slf.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+2:hyb.mon.nfr(i_)) = ...
            %                 conj(flip(slf.syn{i_}.fsa.(hyb.mon.cp{j_})(1:hyb.mon.nfr(i_)/2-1)));
            %             %
            %             shf.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+1) = 0.0;
            %             shf.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+2:hyb.mon.nfr(i_)) = ...
            %                 conj(flip(shf.syn{i_}.fsa.(hyb.mon.cp{j_})(1:hyb.mon.nfr(i_)/2-1)));
            %             %
            %             hyb.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+1) = 0.0;
            %             hyb.syn{i_}.fsa.(hyb.mon.cp{j_})(hyb.mon.nfr(i_)/2+2:hyb.mon.nfr(i_)) = ...
            %                 conj(flip(hyb.syn{i_}.fsa.(hyb.mon.cp{j_})(1:hyb.mon.nfr(i_)/2-1)));
            %             % _Inverse Fourier Transform_
            %             slf.syn{i_}.tha.(slf.mon.cp{j_}) = ...
            %                 real(ifft(slf.syn{i_}.fsa.(slf.mon.cp{j_})))./...
            %                 slf.mon.dtm(i_);
            %             slf.syn{i_}.tha.(slf.mon.cp{j_}) = ...
            %                 slf.syn{i_}.tha.(slf.mon.cp{j_})(1:slf.mon.ntm(i_));
            %             %
            %             shf.syn{i_}.tha.(shf.mon.cp{j_}) = ...
            %                 real(ifft(shf.syn{i_}.fsa.(shf.mon.cp{j_})))./...
            %                 shf.mon.dtm(i_);
            %             shf.syn{i_}.tha.(shf.mon.cp{j_}) = ...
            %                 shf.syn{i_}.tha.(shf.mon.cp{j_})(1:shf.mon.ntm(i_));
            %             %
            %             hyb.syn{i_}.tha.(hyb.mon.cp{j_}) = ...
            %                 real(ifft(hyb.syn{i_}.fsa.(hyb.mon.cp{j_})))./...
            %                 hyb.mon.dtm(i_);
            %             hyb.syn{i_}.tha.(hyb.mon.cp{j_}) = ...
            %                 hyb.syn{i_}.tha.(hyb.mon.cp{j_})(1:hyb.mon.ntm(i_));
            %% *TIME-HISTORIES*
            slf.syn{i_}.tha.(slf.mon.cp{j_}) = super_ifft(slf.mon.dtm(i_),...
                slf.mon.ntm(i_),slf.syn{i_}.fsa.(slf.mon.cp{j_}));
            %
            shf.syn{i_}.tha.(shf.mon.cp{j_}) = super_ifft(shf.mon.dtm(i_),...
                shf.mon.ntm(i_),shf.syn{i_}.fsa.(shf.mon.cp{j_}));
            %
            hyb.syn{i_}.tha.(hyb.mon.cp{j_}) = super_ifft(hyb.mon.dtm(i_),...
                hyb.mon.ntm(i_),hyb.syn{i_}.fsa.(hyb.mon.cp{j_}));
            %% *VELOCITY-DISPLACEMENTS*
            [~,slf.syn{i_}.thv.(slf.mon.cp{j_}),slf.syn{i_}.thd.(slf.mon.cp{j_})] = ...
                band_pass_filter(slf.mon.dtm(i_),slf.syn{i_}.tha.(slf.mon.cp{j_}),[],[]);
            %
            [~,shf.syn{i_}.thv.(shf.mon.cp{j_}),shf.syn{i_}.thd.(shf.mon.cp{j_})] = ...
                band_pass_filter(shf.mon.dtm(i_),shf.syn{i_}.tha.(shf.mon.cp{j_}),[],[]);
            %
            [~,hyb.syn{i_}.thv.(hyb.mon.cp{j_}),hyb.syn{i_}.thd.(hyb.mon.cp{j_})] = ...
                band_pass_filter(hyb.mon.dtm(i_),hyb.syn{i_}.tha.(hyb.mon.cp{j_}),[],[]);            
%             slf.syn{i_}.thd.(slf.mon.cp{j_}) = ...
%                 slf.mon.dtm(i_)*cumtrapz(slf.syn{i_}.thv.(slf.mon.cp{j_}));
%             %
%             shf.syn{i_}.thv.(shf.mon.cp{j_}) = ...
%                 shf.mon.dtm(i_)*cumtrapz(shf.syn{i_}.tha.(shf.mon.cp{j_}));
%             shf.syn{i_}.thd.(shf.mon.cp{j_}) = ...
%                 shf.mon.dtm(i_)*cumtrapz(shf.syn{i_}.thv.(shf.mon.cp{j_}));
%             %
%             hyb.syn{i_}.thv.(hyb.mon.cp{j_}) = ...
%                 hyb.mon.dtm(i_)*cumtrapz(hyb.syn{i_}.tha.(hyb.mon.cp{j_}));
%             hyb.syn{i_}.thd.(hyb.mon.cp{j_}) = ...
%                 hyb.mon.dtm(i_)*cumtrapz(hyb.syn{i_}.thv.(hyb.mon.cp{j_}));
%             % _PGA-PGV-PGD_
%             [slf.syn{i_}.pga.(slf.mon.cp{j_})(1),slf.syn{i_}.pga.(slf.mon.cp{j_})(2),...
%                 slf.syn{i_}.pgv.(slf.mon.cp{j_})(1),slf.syn{i_}.pgv.(slf.mon.cp{j_})(2),...
%                 slf.syn{i_}.pgd.(slf.mon.cp{j_})(1),slf.syn{i_}.pgd.(slf.mon.cp{j_})(2)] = ...
%                 PGAVD_eval(slf.mon.dtm(i_),slf.syn{i_}.tha.(slf.mon.cp{j_}));
%             %
%             [shf.syn{i_}.pga.(shf.mon.cp{j_})(1),shf.syn{i_}.pga.(shf.mon.cp{j_})(2),...
%                 shf.syn{i_}.pgv.(shf.mon.cp{j_})(1),shf.syn{i_}.pgv.(shf.mon.cp{j_})(2),...
%                 shf.syn{i_}.pgd.(shf.mon.cp{j_})(1),shf.syn{i_}.pgd.(shf.mon.cp{j_})(2)] = ...
%                 PGAVD_eval(shf.mon.dtm(i_),shf.syn{i_}.tha.(shf.mon.cp{j_}));
%             %
%             [hyb.syn{i_}.pga.(hyb.mon.cp{j_})(1),hyb.syn{i_}.pga.(hyb.mon.cp{j_})(2),...
%                 hyb.syn{i_}.pgv.(hyb.mon.cp{j_})(1),hyb.syn{i_}.pgv.(hyb.mon.cp{j_})(2),...
%                 hyb.syn{i_}.pgd.(hyb.mon.cp{j_})(1),hyb.syn{i_}.pgd.(hyb.mon.cp{j_})(2)] = ...
%                 PGAVD_eval(hyb.mon.dtm(i_),hyb.syn{i_}.tha.(hyb.mon.cp{j_}));
%             % _Resonse Spectrum_
%             [slf.syn{i_}.rsd.(slf.mon.cp{j_}),~,~,slf.syn{i_}.psa.(slf.mon.cp{j_}),~] = ...
%                 newmark_sd(slf.syn{i_}.tha.(slf.mon.cp{j_}),slf.mon.dtm(i_),slf.mon.vTn,...
%                 slf.mon.zeta);
%             %
%             [shf.syn{i_}.rsd.(shf.mon.cp{j_}),~,~,shf.syn{i_}.psa.(shf.mon.cp{j_}),~] = ...
%                 newmark_sd(shf.syn{i_}.tha.(shf.mon.cp{j_}),shf.mon.dtm(i_),shf.mon.vTn,...
%                 shf.mon.zeta);
%             %
%             [hyb.syn{i_}.rsd.(hyb.mon.cp{j_}),~,~,hyb.syn{i_}.psa.(hyb.mon.cp{j_}),~] = ...
%                 newmark_sd(hyb.syn{i_}.tha.(hyb.mon.cp{j_}),hyb.mon.dtm(i_),hyb.mon.vTn,...
%                 hyb.mon.zeta);
        end
    end
    %% OUTPUT
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hyb;
    return
end