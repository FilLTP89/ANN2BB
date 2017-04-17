%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_write_maps_: function to write contour maps on csv file and on
% different procs
%% *N.B.*
% Need for:
% _syn2ann_run.m,rotationComponent2D.m,super_csvwrite.m_

function [varargout] = syn2ann_write_maps(varargin)
    NJB = varargin{1};
    fprintf('JOB: %u',NJB);
    status = 1;
    load(sprintf('/mssmat2/home/gattif/Documents/PHD_passing_through_polimi/syn2ann/database/maps/syn2ann_input_maps_%u.mat',NJB));
%     load(sprintf('/home/gattif/PHD_passing_through_polimi/syn2ann/src/syn2ann_input_maps_%u.mat',NJB));
    syn2ann_run;
    %% *SPM MAP*
    for i_ = 1:mon.na
        % _SPM_
        map_spm.X{i_,1} = mon.eutm(i_);
        map_spm.Y{i_,1} = mon.nutm(i_);
        % _HYB_
        map_hyb.X{i_,1} = mon.eutm(i_);
        map_hyb.Y{i_,1} = mon.nutm(i_);
        % _PBS_
        map_pbs.X{i_,1} = mon.eutm(i_);
        map_pbs.Y{i_,1} = mon.nutm(i_);
    end
    
    %% *WRITE PGA MAPS*
    if any(strcmpi(mon.map.typ,'pga'))
        for i_=1:mon.na
            for j_=1:mon.nc
                % _SPM_
                map_spm.(strcat('pga_',mon.cp{j_})){i_,1} = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.pga.(mon.cp{j_})(2));
                % _HYB_
                map_hyb.(strcat('pga_',mon.cp{j_})){i_,1} = ...
                    abs(hbs.bst.syn{i_}.pga.(mon.cp{j_})(2));
                % _PBS_
                map_pbs.(strcat('pga_',mon.cp{j_})){i_,1} = ...
                    abs(pbs.org.syn{i_}.pga.(mon.cp{j_})(2));
            end
            % _SPM - FP/FN DIRECTIONS_
            [map_spm.pga_pl{i_,1},map_spm.pga_pp{i_,1}] = ...
                rotationComponent2D(spm.sps.ew.syn{i_}.pga.ew(2),...
                spm.sps.ns.syn{i_}.pga.ns(2),mon.map.stk);
            map_spm.pga_pl{i_,1} = abs(map_spm.pga_pl{i_,1});
            map_spm.pga_pp{i_,1} = abs(map_spm.pga_pp{i_,1});
            map_spm.pga_gh{i_,1} = ...
                sqrt(abs(spm.sps.ew.syn{i_}.pga.ew(2)*spm.sps.ns.syn{i_}.pga.ns(2)));
            % _HYB - FP/FN DIRECTIONS_
            [map_hyb.pga_pl{i_,1},map_hyb.pga_pp{i_,1}] = ...
                rotationComponent2D(hbs.bst.syn{i_}.pga.ew(2),...
                hbs.bst.syn{i_}.pga.ns(2),mon.map.stk);
            map_hyb.pga_pl{i_,1} = abs(map_hyb.pga_pl{i_,1});
            map_hyb.pga_pp{i_,1} = abs(map_hyb.pga_pp{i_,1});
            map_hyb.pga_gh{i_,1} = ...
                sqrt(abs(hbs.bst.syn{i_}.pga.ew(2)*hbs.bst.syn{i_}.pga.ns(2)));
            % _PBS - FP/FN DIRECTIONS_
            [map_pbs.pga_pl{i_,1},map_pbs.pga_pp{i_,1}] = ...
                rotationComponent2D(pbs.org.syn{i_}.pga.ew(2),...
                pbs.org.syn{i_}.pga.ns(2),mon.map.stk);
            map_pbs.pga_pl{i_,1} = abs(map_pbs.pga_pl{i_,1});
            map_pbs.pga_pp{i_,1} = abs(map_pbs.pga_pp{i_,1});
            map_pbs.pga_gh{i_,1} = ...
                sqrt(abs(pbs.org.syn{i_}.pga.ew(2)*pbs.org.syn{i_}.pga.ns(2)));
        end
        % _SPM_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_spm_pga_%u.csv',mon.map.tag,NJB)),...
            [map_spm.X,map_spm.Y,map_spm.pga_ew,map_spm.pga_ns,...
            map_spm.pga_gh,map_spm.pga_pl,map_spm.pga_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
        % _HYB_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_hyb_pga_%u.csv',mon.map.tag,NJB)),...
            [map_hyb.X,map_hyb.Y,map_hyb.pga_ew,map_hyb.pga_ns,...
            map_hyb.pga_gh,map_hyb.pga_pl,map_hyb.pga_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
        % _PBS_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_pbs_pga_%u.csv',mon.map.tag,NJB)),...
            [map_pbs.X,map_pbs.Y,map_pbs.pga_ew,map_pbs.pga_ns,...
            map_pbs.pga_gh,map_pbs.pga_pl,map_pbs.pga_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
    end
    
    %% *WRITE PGV MAPS*
    if any(strcmpi(mon.map.typ,'pgv'))
        for i_=1:mon.na
            for j_=1:mon.nc
                % _SPM_
                map_spm.(strcat('pgv_',mon.cp{j_})){i_,1} = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.pgv.(mon.cp{j_})(2));
                % _HYB_
                map_hyb.(strcat('pgv_',mon.cp{j_})){i_,1} = ...
                    abs(hbs.bst.syn{i_}.pgv.(mon.cp{j_})(2));
                % _PBS_
                map_pbs.(strcat('pgv_',mon.cp{j_})){i_,1} = ...
                    abs(pbs.org.syn{i_}.pgv.(mon.cp{j_})(2));
            end
            % _SPM - FP/FN DIRECTIONS_
            [map_spm.pgv_pl{i_,1},map_spm.pgv_pp{i_,1}] = ...
                rotationComponent2D(spm.sps.ew.syn{i_}.pgv.ew(2),...
                spm.sps.ns.syn{i_}.pgv.ns(2),mon.map.stk);
            map_spm.pgv_pl{i_,1} = abs(map_spm.pgv_pl{i_,1});
            map_spm.pgv_pp{i_,1} = abs(map_spm.pgv_pp{i_,1});
            map_spm.pgv_gh{i_,1} = ...
                sqrt(abs(spm.sps.ew.syn{i_}.pgv.ew(2)*spm.sps.ns.syn{i_}.pgv.ns(2)));
            % _HYB - FP/FN DIRECTIONS_
            [map_hyb.pgv_pl{i_,1},map_hyb.pgv_pp{i_,1}] = ...
                rotationComponent2D(hbs.bst.syn{i_}.pgv.ew(2),...
                hbs.bst.syn{i_}.pgv.ns(2),mon.map.stk);
            map_hyb.pgv_pl{i_,1} = abs(map_hyb.pgv_pl{i_,1});
            map_hyb.pgv_pp{i_,1} = abs(map_hyb.pgv_pp{i_,1});
            map_hyb.pgv_gh{i_,1} = ...
                sqrt(abs(hbs.bst.syn{i_}.pgv.ew(2)*hbs.bst.syn{i_}.pgv.ns(2)));
            % _PBS - FP/FN DIRECTIONS_
            [map_pbs.pgv_pl{i_,1},map_pbs.pgv_pp{i_,1}] = ...
                rotationComponent2D(pbs.org.syn{i_}.pgv.ew(2),...
                pbs.org.syn{i_}.pgv.ns(2),mon.map.stk);
            map_pbs.pgv_pl{i_,1} = abs(map_pbs.pgv_pl{i_,1});
            map_pbs.pgv_pp{i_,1} = abs(map_pbs.pgv_pp{i_,1});
            map_pbs.pgv_gh{i_,1} = ...
                sqrt(abs(pbs.org.syn{i_}.pgv.ew(2)*pbs.org.syn{i_}.pgv.ns(2)));
        end
        % _SPM_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_spm_pgv_%u.csv',mon.map.tag,NJB)),...
            [map_spm.X,map_spm.Y,map_spm.pgv_ew,map_spm.pgv_ns,...
            map_spm.pgv_gh,map_spm.pgv_pl,map_spm.pgv_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
        % _HYB_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_hyb_pgv_%u.csv',mon.map.tag,NJB)),...
            [map_hyb.X,map_hyb.Y,map_hyb.pgv_ew,map_hyb.pgv_ns,...
            map_hyb.pgv_gh,map_hyb.pgv_pl,map_hyb.pgv_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
        % _PBS_
        super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_pbs_pgv_%u.csv',mon.map.tag,NJB)),...
            [map_pbs.X,map_pbs.Y,map_pbs.pgv_ew,map_pbs.pgv_ns,...
            map_pbs.pgv_gh,map_pbs.pgv_pl,map_pbs.pgv_pp],...
            '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
            {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
    end
    
    
    %% *WRITE PSA MAPS*
    if any(strcmpi(mon.map.typ,'psa'))
        for k_=1:numel(mon.map.vTn.psa)
            idx = find(abs(spm.sps.ew.mon.vTn-mon.map.vTn.psa(k_))<1e-2);
            if isempty(idx)
                keyboard
            end
            vTn = round(mon.map.vTn.psa(k_)*100);
            for i_=1:mon.na
                for j_=1:mon.nc
                    % _SPM_
                    map_spm.(strcat('psa_',mon.cp{j_},'_',num2str(vTn))){i_,1} = ...
                        abs(spm.sps.(mon.cp{j_}).syn{i_}.psa.(mon.cp{j_})(idx));
                    % _HYB_
                    map_hyb.(strcat('psa_',mon.cp{j_},'_',num2str(vTn))){i_,1} = ...
                        abs(hbs.bst.syn{i_}.psa.(mon.cp{j_})(idx));
                    % _PBS_
                    map_pbs.(strcat('psa_',mon.cp{j_},'_',num2str(vTn))){i_,1} = ...
                        abs(pbs.org.syn{i_}.psa.(mon.cp{j_})(idx));
                end
                % _SPM_
                map_spm.(strcat('psa_gh_',num2str(vTn))){i_,1} =...
                    sqrt(abs(spm.sps.ew.syn{i_}.psa.ew(idx)*spm.sps.ns.syn{i_}.psa.ns(idx)));
                % _HYB_
                map_hyb.(strcat('psa_gh_',num2str(vTn))){i_,1} =...
                    sqrt(abs(hbs.bst.syn{i_}.psa.ew(idx)*hbs.bst.syn{i_}.psa.ns(idx)));
                % _PBS_
                map_pbs.(strcat('psa_gh_',num2str(vTn))){i_,1} =...
                    sqrt(abs(pbs.org.syn{i_}.psa.ew(idx)*pbs.org.syn{i_}.psa.ns(idx)));
            end
            % _SPM_
            super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_spm_psa_%u_%u.csv',...
                mon.map.tag,round(mon.map.vTn.psa(k_)*100),NJB)),...
                [map_spm.X,map_spm.Y,map_spm.(strcat('psa_ew_',num2str(vTn))),...
                map_spm.(strcat('psa_ns_',num2str(vTn))),...
                map_spm.(strcat('psa_gh_',num2str(vTn)))],...
                '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
                {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]'});
            % _HYB_
            super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_hyb_psa_%u_%u.csv',...
                mon.map.tag,round(mon.map.vTn.psa(k_)*100),NJB)),...
                [map_hyb.X,map_hyb.Y,map_hyb.(strcat('psa_ew_',num2str(vTn))),...
                map_hyb.(strcat('psa_ns_',num2str(vTn))),...
                map_hyb.(strcat('psa_gh_',num2str(vTn)))],...
                '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
                {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]'});
            % _PBS_
            super_csvwrite(fullfile(mon.map.fnm,sprintf('%s_pbs_psa_%u_%u.csv',...
                mon.map.tag,round(mon.map.vTn.psa(k_)*100),NJB)),...
                [map_pbs.X,map_pbs.Y,map_pbs.(strcat('psa_ew_',num2str(vTn))),...
                map_pbs.(strcat('psa_ns_',num2str(vTn))),...
                map_pbs.(strcat('psa_gh_',num2str(vTn)))],...
                '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
                {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]'});
        end
    end
    varargout{1} = status;
end
