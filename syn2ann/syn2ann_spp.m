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
    out_sel=1:2;
    if nargin>1
        out_sel=varargin{2};
    end
    %% *SD and PSA SPECTRA*
    if any(out_sel==1)
        % minimum natural period
        Tn_min   = 0;
        % maximum natural period
        Tn_max   = 5;
        % natural period step
        dTn      = 0.05;
        sas.mon.vTn  = (Tn_min:dTn:Tn_max)';
        sas.mon.nT   = numel(sas.mon.vTn);
        sas.mon.zeta = 0.05;
        for i_ = 1:sas.mon.na
            for j_ = 1:sas.mon.nc
                [sas.syn{i_}.psa.(sas.mon.cp{j_}),sas.syn{i_}.rsd.(sas.mon.cp{j_})] = ...
                    SDOF_response(sas.syn{i_}.tha.(sas.mon.cp{j_}),...
                    sas.mon.dtm(i_),sas.mon.vTn,...
                    sas.mon.zeta,[1,2]);
            end
        end
    end
    %% *FOURIER SPECTRUM*
    if any(out_sel==2)
        
        sas.mon.vfr = cell(sas.mon.na,1);
 
        for i_ = 1:sas.mon.na
            % frequency vector
            sas.mon.vfr{i_} = super_fft(sas.mon.dtm(i_),...
                sas.syn{i_}.tha.(sas.mon.cp{1}),0,1);
            % number of frequency-points
            sas.mon.nfr(i_) = numel(sas.mon.vfr{i_});
            % frequency step
            sas.mon.dfr(i_) = mean(diff(sas.mon.vfr{i_}));
            % Nyquist frequency index
            sas.mon.nNy(i_) = floor(0.5/sas.mon.dtm(i_)/sas.mon.dfr(i_)+0.5);
            for j_ = 1:sas.mon.nc
                sas.syn{i_}.fsa.(sas.mon.cp{j_}) = super_fft(sas.mon.dtm(i_),...
                    sas.syn{i_}.tha.(sas.mon.cp{j_}),0,4);
            end
        end
    end
    
    %% OUTPUT
    varargout{1} = sas;
    return
end
