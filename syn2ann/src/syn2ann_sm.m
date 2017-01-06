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
    
    %% *SETUP*
    hbs = varargin{1};
    trs = varargin{2};
    spm.mon.cp = hbs.mon.cp;
    
    %% *SPECTRAL SCALING*
    for i_ = 1:hbs.mon.na
        for j_=1:hbs.mon.nc
            spm.(hbs.mon.cp{j_}).mon.cp = hbs.mon.cp(j_);
            spm.(hbs.mon.cp{j_}).mon.nc = 1;
            spm.(hbs.mon.cp{j_}).mon.na = hbs.mon.na;
            %
            % _spectral scaling_
            %
            
%             [out_t,out_acc,out_vel,out_dis,out_T,out_Se,out_freq,out_FAS,n_] = ...
%        SpectralMatching(target_Se,hybrid_acc,niter,tol_upp,tol_low,tol_upp_pga,tol_low_pga,i);
   
            [spm.(hbs.mon.cp{j_}).mon.vtm{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).mon.vTn{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).mon.vfr{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.fsa.(hbs.mon.cp{j_}),NC] = ...
                spectral_matching_new(...
                [trs.(trs.mon.cp{j_}).mon.vTn(:),...
                trs.(trs.mon.cp{j_}).syn{i_}.psa.(trs.mon.cp{j_})(:)*100],...
                [hbs.mon.vtm{i_}(:),hbs.syn{i_}.tha.(hbs.mon.cp{j_})(:)],...
                trs.(trs.mon.cp{j_}).mon.nit,trs.(trs.mon.cp{j_}).mon.tol.hgh.psa,...
                trs.(trs.mon.cp{j_}).mon.tol.low.psa,...
                trs.(trs.mon.cp{j_}).mon.tol.hgh.pga,...
                trs.(trs.mon.cp{j_}).mon.tol.low.pga,trs.mon.cp{j_});
            
            spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}) = ...
                spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_})./100;
            
            spm.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_}) = ...
                spm.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_})./100;
            
            spm.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_}) = ...
                spm.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_})./100;
            
            spm.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}) = ...
                spm.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_})./100;
            
            spm.(hbs.mon.cp{j_}).mon.dtm(i_) = mean(diff(spm.(hbs.mon.cp{j_}).mon.vtm{i_}));
            %
            spm.(hbs.mon.cp{j_}).mon.ntm(i_) = numel(spm.(hbs.mon.cp{j_}).mon.vtm{i_});
            %
%             spm.(hbs.mon.cp{j_}).mon.vtm{i_} = spm.(hbs.mon.cp{j_}).mon.dtm(i_)*...
%                 (0:spm.(hbs.mon.cp{j_}).mon.ntm(i_)-1)'+hbs.mon.vtm{i_}(1);
        end
    end
    %% *OUTPUT*
    varargout{1} = spm;
    return
end
