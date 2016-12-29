%% *Compute response/fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _sp_spectra_: function to compute response/fourier spectra for a set
% of synthetics records generated with Sabetta and Pugliese method. 
%% INPUT:  
%% 
% * _sps    = structure of sp synthetics_
%% OUTPUT: 
% * _sps    = structure of sp synthetics_
%% N.B. 
% Need for _newmark_sd.m_
function [varargout] = sp_spectra(varargin)
    %% SET-UP
    sps = varargin{1};
    % _response spectra parameters_
    % minimum natural period
    Tn_min   = 0;
    % maximum natural period
    Tn_max   = 5;
    % natural period step
    dTn      = 0.005;
    sps.mon.vTn  = (Tn_min:dTn:Tn_max)';
    sps.mon.nT   = numel(sps.mon.vTn);
    sps.mon.zeta = 0.05;
    % _fourier spectra parameters_
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
    %% SD and PSA SPECTRA
    for i_ = 1:sps.mon.na
        for j_ = 1:sps.mon.nc
            [sps.syn{i_}.rsd.(sps.mon.cp{j_}),~,~,sps.syn{i_}.psa.(sps.mon.cp{j_}),~] = ...
                newmark_sd(sps.syn{i_}.tha.(sps.mon.cp{j_}),sps.mon.dtm(i_),...
                sps.mon.vTn,sps.mon.zeta);
        end
    end
    %% FOURIER SPECTRUM
    for i_ = 1:sps.mon.na
        for j_ = 1:sps.mon.nc
            [~,~,~,sps.syn{i_}.fsa.(sps.mon.cp{j_}),~,~] = ...
                super_fft(sps.mon.dtm(i_),sps.syn{i_}.tha.(sps.mon.cp{j_}),0);
        end
    end
    %% OUTPUT
    varargout{1} = sps;
    return
end
