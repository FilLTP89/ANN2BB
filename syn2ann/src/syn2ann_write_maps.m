%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2amon.na_write_map_spms_: function to write contour map_spms
%% *N.B.*
% Need for:
% __

function [varargout] = syn2ann_write_maps(varargin)
    NJB = varargin{1};
    fprintf('JOB: %u',NJB); 
    status = 1;
        load(sprintf('syn2ann_input_maps_%u.mat',NJB));
        syn2ann_run;
        %% *SPM MAP*
    %     [map_spm(1:mon.na).Geometry] = deal('Point');
    %     [map_hyb(1:mon.na).Geometry] = deal('Point');

    %     for i_ = 1:mon.na
    %         [map_spm(i_).X] = mon.eutm(i_);
    %         [map_spm(i_).Y] = mon.nutm(i_);
    % %         [map_hyb(i_).X] = mon.eutm(i_);
    % %         [map_hyb(i_).Y] = mon.nutm(i_);
    %     end
        for i_ = 1:mon.na
            map_spm.X{i_,1} = mon.eutm(i_);
            map_spm.Y{i_,1} = mon.nutm(i_);
        end
        
        %% *WRITE PGA MAPS*
        if any(strcmpi(mon.map.typ,'pga'))
            for i_=1:mon.na
                for j_=1:mon.nc
                    map_spm.(strcat('R',mon.cp{j_})){i_,1} = ...
                        abs(spm.sps.(mon.cp{j_}).syn{i_}.pga.(mon.cp{j_})(2));
    %                 map_hyb(i_).(strcat('R',mon.cp{j_})) = ...
    %                     abs(hbs.bst.syn{i_}.pga.(mon.cp{j_})(2));
                end
                [map_spm.Rp{i_,1},map_spm.Ro{i_,1}] = ...
                    rotationComponent2D(spm.sps.ew.syn{i_}.pga.ew(2),...
                    spm.sps.ns.syn{i_}.pga.ns(2),mon.map.stk);
                map_spm.Rp{i_,1} = abs(map_spm.Rp{i_,1});
                map_spm.Ro{i_,1} = abs(map_spm.Ro{i_,1});
                map_spm.Rgh{i_,1} = ...
                    sqrt(abs(spm.sps.ew.syn{i_}.pga.ew(2)*spm.sps.ns.syn{i_}.pga.ns(2)));
                %
    %             [map_hyb(i_).Rp,map_hyb(i_).Ro] = ...
    %                 rotationComponent2D(hbs.bst.syn{i_}.pga.ew(2),...
    %                 hbs.bst.syn{i_}.pga.ns(2),mon.map.stk);
    %             map_hyb(i_).Rp = abs(map_hyb(i_).Rp);
    %             map_hyb(i_).Ro = abs(map_hyb(i_).Ro);
    %             map_hyb(i_).Rgh = ...
    %                 sqrt(abs(hbs.bst.syn{i_}.pga.ew(2)*hbs.bst.syn{i_}.pga.ns(2)));
                
            end
            super_csvwrite(sprintf('%s_PGA_%u.csv',mon.map.fnm,NJB),...
                [map_spm.X,map_spm.Y,map_spm.Rew,map_spm.Rns,...
                    map_spm.Rgh,map_spm.Rp,map_spm.Ro],...
                    '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
                {'E-UTM [m]','N-UTM [m]','Rew [m/s2]','Rns [m/s2]','Rgh [m/s2]','Rp [m/s2]','Ro [m/s2]'});
    %         shapewrite(map_spm,strcat(mon.map.fnm,'_spm_pga'));
    %         shapewrite(map_hyb,strcat(mon.map.fnm,'_hyb_pga'));
        end
        
    %     %% *WRITE PGV MAPS*
    %     if any(strcmpi(mon.map.typ,'pgv'))
    %         for i_=1:mon.na
    %             for j_=1:mon.nc
    %                 map_spm(i_).(strcat('R',mon.cp{j_})) = ...
    %                     abs(spm.sps.(mon.cp{j_}).syn{i_}.pgv.(mon.cp{j_})(2));
    %             end
    %             [map_spm(i_).Rp,map_spm(i_).Ro] = ...
    %                 rotationComponent2D(spm.sps.ew.syn{i_}.pgv.ew(2),...
    %                 spm.sps.ns.syn{i_}.pgv.ns(2),mon.map.stk);
    %             map_spm(i_).Rp = abs(map_spm(i_).Rp);
    %             map_spm(i_).Ro = abs(map_spm(i_).Ro);
    %             map_spm(i_).Rgh = ...
    %                 sqrt(abs(spm.sps.ew.syn{i_}.pgv.ew(2)*spm.sps.ns.syn{i_}.pgv.ns(2)));
    %         end
    %         shapewrite(map_spm,strcat(mon.map.fnm,'_pgv'));
    %     end
    %     
    %     %% *WRITE PGD MAPS*
    %     if any(strcmpi(mon.map.typ,'pgd'))
    %         for i_=1:mon.na
    %             for j_=1:mon.nc
    %                 map_spm(i_).(strcat('R',mon.cp{j_})) = ...
    %                     abs(spm.sps.(mon.cp{j_}).syn{i_}.pgd.(mon.cp{j_})(2));
    %             end
    %             [map_spm(i_).Rp,map_spm(i_).Ro] = ...
    %                 rotationComponent2D(spm.sps.ew.syn{i_}.pgd.ew(2),...
    %                 spm.sps.ns.syn{i_}.pgd.ns(2),mon.map.stk);
    %             map_spm(i_).Rp = abs(map_spm(i_).Rp);
    %             map_spm(i_).Ro = abs(map_spm(i_).Ro);
    %             map_spm(i_).Rgh = ...
    %                 sqrt(abs(spm.sps.ew.syn{i_}.pgd.ew(2)*spm.sps.ns.syn{i_}.pgd.ns(2)));
    %         end
    %         shapewrite(map_spm,strcat(mon.map.fnm,'_pgd'));
    %     end
    %     
    %     %% *WRITE RSD MAPS*
    %     if any(strcmpi(mon.map.typ,'rsd'))
    %         idx = find(spm.sps.ew.mon.vTn==mon.map.vTn.rsd);
    %         for i_=1:mon.na
    %             for j_=1:mon.nc
    %                 map_spm(i_).(strcat('R',mon.cp{j_})) = spm.sps.(mon.cp{j_}).syn{i_}.rsd.(mon.cp{j_})(idx);
    %             end
    %             map_spm(i_).Rgh = sqrt(spm.sps.ew.syn{i_}.rsd.ew(idx)*spm.sps.ns.syn{i_}.rsd.ns(idx));
    %         end
    %         shapewrite(map_spm,strcat(mon.map.fnm,'_rsd'));
    %     end
    %     
    %     %% *WRITE PSA MAPS*
    %     if any(strcmpi(mon.map.typ,'psa'))
    %         idx = find(spm.sps.ew.mon.vTn==mon.map.vTn.psa);
    %         for i_=1:mon.na
    %             for j_=1:mon.nc
    %                 map_spm(i_).(strcat('R',mon.cp{j_})) = ...
    %                     abs(spm.sps.(mon.cp{j_}).syn{i_}.psa.(mon.cp{j_})(idx));
    %             end
    %             map_spm(i_).Rgh =...
    %                 sqrt(abs(spm.sps.ew.syn{i_}.psa.ew(idx)*spm.sps.ns.syn{i_}.psa.ns(idx)));
    %         end
    %         shapewrite(map_spm,strcat(mon.map.fnm,'_psa'));
    %     end
    varargout{1} = status;
end
