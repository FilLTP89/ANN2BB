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
    %% *SET-UP*
    slf = varargin{1};
    shf = varargin{2};
    
    hyb.mon = slf.mon;
    hyb.mtd = shf.mtd;
    hyb.syn = cell(slf.mon.na,1);

    %% *HYBRIDIZATION*
    [slf,shf,hyb] = lfhf_hybridator(slf,shf,hyb);
    %% *OUTPUT*
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hyb;
    return
end
%% *Hybridator*
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
    %% *SET-UP*
    slf = varargin{1};
    shf = varargin{2};
    hyb = varargin{3};
    
    %% *BROAD-BAND SIGNALS*
    switch hyb.mon.hyb
        case 'original'
            %
            % _LOOP OVER THE MONITORS_
            %
            for i_ = 1:slf.mon.na
                %
                % _LOOP OVER THE DIRECTIONS_
                %
                for j_ = 1:slf.mon.nc
                    % _current direction_
                    cpp = slf.mon.cp{j_};
                    % _find corner frequencies_
                    nfa = round(hyb.mtd.(cpp)(i_,1)/hyb.mon.dfr(i_));
                    nfb = round(hyb.mtd.(cpp)(i_,2)/hyb.mon.dfr(i_));
                    fac = pi./(hyb.mon.dfr(i_)*(nfb-nfa-1));
                    
                    %% *FOURIER TRANSFORM*
                    slf.syn{i_}.tha.(cpp)(slf.mon.ntm+1:hyb.mon.nfr(i_))=0;
                    shf.syn{i_}.tha.(cpp)(shf.mon.ntm+1:hyb.mon.nfr(i_))=0;
                    
                    slf.syn{i_}.fsa.(cpp) = fft(slf.syn{i_}.tha.(cpp),hyb.mon.nfr(i_)).*...
                        slf.mon.dtm(i_);
                    shf.syn{i_}.fsa.(cpp) = fft(shf.syn{i_}.tha.(cpp),hyb.mon.nfr(i_)).*...
                        shf.mon.dtm(i_);
                    %% *FILTERING*
                    WLF = zeros(hyb.mon.nfr(i_),1);
                    WLF(1:nfa) = 1.0;
                    WLF(nfa+1:nfb) = 0.5*(1+cos(fac*hyb.mon.dfr(i_)*(0:nfb-nfa-1)));
                    WHF = 1-WLF;
                    %
                    slf.syn{i_}.fsa.(cpp) = WLF.*slf.syn{i_}.fsa.(cpp);
                    %
                    shf.syn{i_}.fsa.(cpp) = WHF.*shf.syn{i_}.fsa.(cpp);
                    %
                    hyb.syn{i_}.fsa.(cpp) = slf.syn{i_}.fsa.(cpp) + shf.syn{i_}.fsa.(cpp);
                    
                    %% *TIME-HISTORIES*
                    
                    slf.syn{i_}.tha.(cpp) = super_ifft(slf.mon.dtm(i_),slf.syn{i_}.fsa.(cpp));
                    %
                    shf.syn{i_}.tha.(cpp) = super_ifft(shf.mon.dtm(i_),shf.syn{i_}.fsa.(cpp));
                    %
                    hyb.syn{i_}.tha.(cpp) = super_ifft(hyb.mon.dtm(i_),hyb.syn{i_}.fsa.(cpp));
                    
                    %% *VELOCITY-DISPLACEMENTS*
                    [~,slf.syn{i_}.thv.(cpp),slf.syn{i_}.thd.(cpp)] = ...
                        idc_tha(slf.mon.dtm(i_),slf.syn{i_}.tha.(cpp));
                    %
                    [~,shf.syn{i_}.thv.(cpp),shf.syn{i_}.thd.(cpp)] = ...
                        idc_tha(shf.mon.dtm(i_),shf.syn{i_}.tha.(cpp));
                    %
                    [~,hyb.syn{i_}.thv.(cpp),hyb.syn{i_}.thd.(cpp)] = ...
                        idc_tha(hyb.mon.dtm(i_),hyb.syn{i_}.tha.(cpp));
                    %
                end
                slf.mon.ntm(i_) = numel(slf.syn{i_}.tha.(cpp));
                slf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:slf.mon.ntm(i_)-1)'};
                %
                shf.mon.ntm(i_) = numel(shf.syn{i_}.tha.(cpp));
                shf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:shf.mon.ntm(i_)-1)'};
                %
                hyb.mon.ntm(i_) = numel(hyb.syn{i_}.tha.(cpp));
                hyb.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:hyb.mon.ntm(i_)-1)'};
                
            end
        case 'butter'
            %
            % _LOOP OVER THE MONITORS_
            %
            for i_ = 1:slf.mon.na
                %
                % _LOOP OVER THE DIRECTIONS_
                %
                for j_ = 1:slf.mon.nc
                    cpp = slf.mon.cp{j_};
                    
                    % _find corner frequencies_
                    nfa = round(hyb.mtd.(cpp)(i_,1)/hyb.mon.dfr(i_));
                    nfb = round(hyb.mtd.(cpp)(i_,2)/hyb.mon.dfr(i_));
                    
                    fac = pi./(hyb.mon.dfr(i_)*(nfb-nfa-1));
                    
                    slf.syn{i_}.tha.(cpp)(slf.mon.ntm+1:hyb.mon.nfr(i_))=0;
                    shf.syn{i_}.tha.(cpp)(shf.mon.ntm+1:hyb.mon.nfr(i_))=0;
                    
                    % [vfr_before,fsd_before] = super_fft(shf.mon.dtm(i_),shf.syn{i_}.thv.(cpp),0,[1,2]);
                    %% *CREATE BUTTERWORTH FILTER*
                    %
                    % _LF FILTER (LOW-PASS)_
                    %
                    [bfb.slf,bfa.slf,~] = create_butter_filter(2,[],hyb.mtd.(cpp)(i_,1),1/2/hyb.mon.dtm(i_));
                    %
                    % _HF FILTER (HIGH-PASS)_
                    %
                    [bfb.shf,bfa.shf,~] = create_butter_filter(2,hyb.mtd.(cpp)(i_,2),[],1/2/hyb.mon.dtm(i_));
                    
                    %% *TIME-HISTORIES*
                    slf.syn{i_}.thd.(cpp) = filtfilt(bfb.slf,bfa.slf,slf.syn{i_}.thd.(cpp));
%                     slf.syn{i_}.thv.(cpp) = avd_diff(slf.mon.dtm(i_),slf.syn{i_}.thd.(cpp));
%                     slf.syn{i_}.tha.(cpp) = avd_diff(slf.mon.dtm(i_),slf.syn{i_}.thv.(cpp));
                    slf.syn{i_}.thv.(cpp) = freq_differentiate(slf.syn{i_}.thd.(cpp),1./slf.mon.dtm(i_));
                    slf.syn{i_}.tha.(cpp) = freq_differentiate(slf.syn{i_}.thv.(cpp),1./slf.mon.dtm(i_));
                    
                    %
                    shf.syn{i_}.thd.(cpp) = filtfilt(bfb.shf,bfa.shf,shf.syn{i_}.thd.(cpp));
                    % shf.syn{i_}.thv.(cpp) = avd_diff(shf.mon.dtm(i_),shf.syn{i_}.thd.(cpp));
                    % shf.syn{i_}.tha.(cpp) = avd_diff(shf.mon.dtm(i_),shf.syn{i_}.thv.(cpp));
                    shf.syn{i_}.thv.(cpp) = vfr_diff(shf.syn{i_}.thd.(cpp),1./shf.mon.dtm(i_));
                    shf.syn{i_}.tha.(cpp) = vfr_diff(shf.syn{i_}.thv.(cpp),1./shf.mon.dtm(i_));
                    
                    keyboard
                    
                    hyb.syn{i_}.tha.(cpp) = slf.syn{i_}.tha.(cpp)+shf.syn{i_}.tha.(cpp);
                    
                    %% *HYBRID VELOCITY-DISPLACEMENTS*
                    %[~,hyb.syn{i_}.thv.(cpp),hyb.syn{i_}.thd.(cpp)] = ...
                    %    idc_tha(hyb.mon.dtm(i_),hyb.syn{i_}.tha.(cpp));
                    %
                    
                    hyb.syn{i_}.thv.(cpp) = freq_integrate(hyb.syn{i_}.tha.(cpp),1./hyb.mon.dtm(i_));
                    hyb.syn{i_}.thd.(cpp) = freq_integrate(hyb.syn{i_}.thv.(cpp),1./hyb.mon.dtm(i_));
                end
                slf.mon.ntm(i_) = numel(slf.syn{i_}.tha.(cpp));
                slf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:slf.mon.ntm(i_)-1)'};
                %
                shf.mon.ntm(i_) = numel(shf.syn{i_}.tha.(cpp));
                shf.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:shf.mon.ntm(i_)-1)'};
                %
                hyb.mon.ntm(i_) = numel(hyb.syn{i_}.tha.(cpp));
                hyb.mon.vtm(i_) = {slf.mon.dtm(i_)*(0:hyb.mon.ntm(i_)-1)'};
                
            end
    end
    %% OUTPUT
    varargout{1} = slf;
    varargout{2} = shf;
    varargout{3} = hyb;
    return
end
