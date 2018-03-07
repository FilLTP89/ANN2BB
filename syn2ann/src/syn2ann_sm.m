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
            [spm.(hbs.mon.cp{j_}).mon.vtm{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thv.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).syn{i_}.thd.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).mon.vTn{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_}),...
                spm.(hbs.mon.cp{j_}).mon.vfr{i_},...
                spm.(hbs.mon.cp{j_}).syn{i_}.fsa.(hbs.mon.cp{j_}),NC] = ...
                spectral_matching_review(...
                [trs.(trs.mon.cp{j_}).mon.vTn(:),...
                trs.(trs.mon.cp{j_}).syn{i_}.psa.(trs.mon.cp{j_})(:)*100],...
                [hbs.mon.vtm{i_}(:),hbs.syn{i_}.tha.(hbs.mon.cp{j_})(:)],...
                trs.mon.cp{j_},trs.(hbs.mon.cp{j_}).mon.TnC);
            %
            % _post process_
            %            
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
%             if NC<27
%                 disp(' ')
%                 disp('Info: Match is not found to be spectrally satisfactory.')
%                 disp('      You may want to observe the agreement, increase the iteration number and re-run.')
%                 
%                 %%% COMMENT IS OPTIONAL, in this version program does not need any
%                 %%% input from the user.
%                 %        disp('    Press any key to continue...')
%                 %        pause
%                 %%% COMMENT IS OPTIONAL
%                 
%                 %MAIN
%             else
%                 disp('Info: Match is acceptable.')
%             end
        end
    end
    %% *OUTPUT*
    varargout{1} = spm;
    return
end
