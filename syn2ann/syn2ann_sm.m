%% *Spectral Matching*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_sm_: function to match hybrid synthetic
% spectra with trained ANN results
%% INPUT:
% * _hbs (structure of hybrid synthetics)_
% * _trs (trained/simulated ann structure)_
%% OUTPUT:
% * _spm (spectral matched structure)_
%% N.B.
function [varargout] = syn2ann_sm(varargin)
    
    %% SET-UP
    hbs = varargin{1};
    trs = varargin{2};
    smr.mon.cp = hbs.mon.cp;
    %% SPECTRAL SCALING
    for i_ = 1:hbs.mon.na
        for j_=1:hbs.mon.nc
            smr.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp(j_);
            smr.(hbs.mon.cp{j_}).mon.nc = 1;
            smr.(hbs.mon.cp{j_}).mon.na = hbs.mon.na;
            % _target period - indexes_
            vTn_idx = 1:numel(trs.e.tid);
            %
            % _spectral scaling_
            %
            [smr.(hbs.mon.cp{j_}).mon.dtm(i_),...
                smr.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}),...
                smr.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_}),...
                smr.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_}),...
                smr.(hbs.mon.cp{j_}).mon.vfr{1},smr.(hbs.mon.cp{j_}).mon.vTn,...
                smr.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}),...
                smr.(hbs.mon.cp{j_}).syn{i_}.psd.(hbs.mon.cp{j_})] = ...
                spectral_scaling(hbs.mon.dtm(i_),...
                hbs.syn{i_}.tha.(hbs.mon.cp{j_}),...
                trs.(trs.mon.cp{j_}).syn{i_}.psa.(trs.mon.cp{j_}),...
                trs.(trs.mon.cp{j_}).mon.vTn,vTn_idx);
            %
            smr.(hbs.mon.cp{j_}).mon.ntm(i_) = numel(smr.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
            %
            smr.(hbs.mon.cp{j_}).mon.vtm{i_} = smr.(hbs.mon.cp{j_}).mon.dtm(i_)*...
                (0:smr.(hbs.mon.cp{j_}).mon.ntm(i_)-1);
        end
    end
    %% OUTPUT
    varargout{1} = smr;
    return
end
