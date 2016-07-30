%% *Generate synthetics - Sabetta & Pugliese approach*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_sp96_generator_: function to produce non-stationary accelerograms
% according to the approach by Sabetta & Pugliese, 1996.
%% INPUT:
% * _mon (monitor structure)_
% * _str (string: extra metadata file)
%% OUTPUT:
% * _sps (structure of sp synthetics)_
%% N.B.:
% it has been verified that with ssc = 2 the method of sabetta
% gives as output acceleration THs with an average response spectrum
% compliant with the (median+isig/2) spectrum given by the GMPE
% Need for _sabetta.m_
%% REFERENCES
% @Article{Sabetta_Pugliese_1996,
%   author =  {Sabetta, F. and Pugliese, A.},
%   title =   {{Estimation of Response Spectra and Simulation of Nonstationary Earthquake Ground Motions}},
%   journal = {Bulletin of the Seismological Society of America},
%   year =    {1996},
%   volume =  {86},
%   number =  {2},
%   pages =   {337--352},
% }
function [varargout] = syn2ann_sp96_generator(varargin)
    %% *SET-UP*
    sps.mon = varargin{1};
    sps.mtd = varargin{2};
    %% *SABETTA & PUGLIESE SYNTHETICS*
    
    for i_ = 1:sps.mon.na
        for j_ = 1:sps.mon.nc
            cpp = sps.mon.cp{j_};
            [vtm,sps.syn{i_}.tha.(cpp)] = sabetta(sps.mtd.mw(i_),sps.mon.dep(i_),...
                sps.mtd.scc(i_),sps.mtd.sst(i_),sps.mtd.dtm_sp96(i_),sps.mtd.scl(i_));
            %
            [sps.syn{i_}.tha.(cpp),sps.syn{i_}.thv.(cpp),sps.syn{i_}.thd.(cpp)] = ...
                band_pass_filter(sps.mtd.dtm_sp96(i_),sps.syn{i_}.tha.(cpp),...
                sps.mon.lfr,sps.mon.hfr);
            
        end
        sps.mon.dtm(i_) = sps.mtd.dtm_sp96(i_);
        sps.mon.vtm{i_} = vtm;
        sps.mon.ntm(i_) = numel(sps.mon.vtm{i_});
    end
    %% *OUTPUT*
    varargout{1} = sps;
    return
end
