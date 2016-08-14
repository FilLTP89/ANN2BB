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
    spm.mon.cp = hbs.mon.cp;
    %% SPECTRAL SCALING
    for i_ = 1:hbs.mon.na
        for j_=1:hbs.mon.nc
            spm.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp(j_);
            spm.(hbs.mon.cp{j_}).mon.nc = 1;
            spm.(hbs.mon.cp{j_}).mon.na = hbs.mon.na;
            %
            % _spectral scaling_
            %
%             vTn_idx = 1:numel(trs.e.tid);
            
            [spm.(hbs.mon.cp{j_}).mon.dtm(i_),...
                spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).mon.vTn,...
                spm.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_})] = ...
                spectral_scaling_original(hbs.mon.dtm(i_),...
                hbs.syn{i_}.tha.(hbs.mon.cp{j_}),...
                trs.(trs.mon.cp{j_}).mon.vTn,...
                trs.(trs.mon.cp{j_}).syn{i_}.psa.(trs.mon.cp{j_}));
            
            [spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_})] = ...
                integr_diff_avd(spm.(hbs.mon.cp{j_}).mon.dtm(i_),...
                spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
            
            spm.(hbs.mon.cp{j_}).mon.ntm(i_) = numel(spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
            %
            spm.(hbs.mon.cp{j_}).mon.vtm{i_} = spm.(hbs.mon.cp{j_}).mon.dtm(i_)*...
                (0:spm.(hbs.mon.cp{j_}).mon.ntm(i_)-1);
            %
            if strcmpi(hbs.mon.cp{j_},'z')
                figure;
                hold all;
                plot(trs.z.mon.vTn,trs.z.syn{1}.psa.z*100,'md-');
                plot(spm.z.mon.vTn,spm.z.syn{1}.psa.z*100,'b+-');
            end
        end
    end
    %% *OUTPUT*
    varargout{1} = spm;
    return
end