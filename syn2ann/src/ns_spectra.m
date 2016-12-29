%% *Compute response/fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _ns_spectra_: function to compute response/fourier spectra for a set
% of synthetics records outcomes of numerical simulations
%% INPUT:
% * _nss (structure of numerical simulations)_
%% OUTPUT:
% * _nss (structure of numerical simulations)_
%% N.B.
% Need for _newmark_sd.m_
function [varargout] = ns_spectra(varargin)
    %% SET-UP
    % _numerical simulations structure_
    nss = varargin{1};
    % _response spectra parameters_
    % minimum natural period
    Tn_min   = 0;
    % maximum natural period
    Tn_max   = 5;
    % natural period step
    dTn      = 0.005;
    nss.mon.vTn  = (Tn_min:dTn:Tn_max)';
    nss.mon.nT   = numel(nss.mon.vTn);
    nss.mon.zeta = 0.05;
    % _fourier spectra parameters_
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
    %% SD and PSA SPECTRA
    for i_ = 1:nss.mon.na
        for j_ = 1:nss.mon.nc
            [nss.syn{i_}.rsd.(nss.mon.cp{j_}),~,~,nss.syn{i_}.psa.(nss.mon.cp{j_}),~] = ...
                newmark_sd(nss.syn{i_}.tha.(nss.mon.cp{j_}),nss.mon.dtm(i_),...
                nss.mon.vTn,nss.mon.zeta);
        end
    end
    %% FOURIER SPECTRUM
    for i_ = 1:nss.mon.na
        for j_ = 1:nss.mon.nc
            [~,~,~,nss.syn{i_}.fsa.(nss.mon.cp{j_}),~,~] = ...
                super_fft(nss.mon.dtm(i_),nss.syn{i_}.tha.(nss.mon.cp{j_}),0);
        end
    end
    %% OUTPUT
    varargout{1} = nss;
    return
end
