%% *Compute response/fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_spp_: function to compute response/fourier spectra for a set
% of sasorded sasords
%% INPUT:
% * _sas (structure of sasorded signals)_
%% OUTPUT:
% * _sas (structure of sasorded signals))_
%% N.B.
% Need for _newmark_sd.m_
function [varargout] = syn2ann_spp(varargin)
    %% *SET-UP*
    % _numerical simulations structure_
    sas = varargin{1};
    % _response spectra parameters_
    % minimum natural period
    Tn_min   = 0;
    % maximum natural period
    Tn_max   = 5;
    % natural period step
    dTn      = 0.005;
    sas.mon.vTn  = (Tn_min:dTn:Tn_max)';
    sas.mon.nT   = numel(sas.mon.vTn);
    sas.mon.zeta = 0.05;
    % _fourier spectra parameters_
    for i_ = 1:sas.mon.na
        % Nyquist frequency
        fr_max      = 1/2/sas.mon.dtm(i_);
        % fft points
        sas.mon.nfr(i_) = 2^nextpow2(sas.mon.ntm(i_));
        % frequency period step
        sas.mon.dfr(i_) = 1/sas.mon.dtm(i_)/(sas.mon.nfr(i_)-1);
        % Nyquist frequency index
        sas.mon.nNy(i_) = floor(fr_max/sas.mon.dfr(i_))+1;
        % frequency vector
        sas.mon.vfr{i_} = sas.mon.dfr(i_)*(0:sas.mon.nfr(i_)-1)';
    end
    %% *SD and PSA SPECTRA*
    for i_ = 1:sas.mon.na
        for j_ = 1:sas.mon.nc
            [sas.syn{i_}.rsd.(sas.mon.cp{j_}),...
                ~,~,sas.syn{i_}.psa.(sas.mon.cp{j_}),~] = ...
                newmark_sd(sas.syn{i_}.tha.(sas.mon.cp{j_}),...
                sas.mon.dtm(i_),sas.mon.vTn,...
                sas.mon.zeta);
        end
    end
    %% *FOURIER SPECTRUM*
    for i_ = 1:sas.mon.na
        for j_ = 1:sas.mon.nc
            sas.syn{i_}.fsa.(sas.mon.cp{j_}) = ...
                sas.mon.dtm(i_)*fft(sas.syn{i_}.tha.(sas.mon.cp{j_}),...
                sas.mon.nfr(i_));
        end
    end
    %% OUTPUT
    varargout{1} = sas;
    return
end
